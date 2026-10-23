//! Port of `x_orphan_transduction.py` (`XOrphanTransduction`): sibling / orphan transductions.

use crate::alignments::BamInfo;
use crate::annotation::process_chrm_name;
use crate::bam::{AlignmentFile, PyAln, Record, BAM_CSOFT_CLIP};
use crate::black_list::XBlackList;
use crate::chrom::is_decoy_contig_chrms;
use crate::clip_disc_filter::{py_head, py_tail, right_clip_pos};
use crate::globals as gv;
use crate::intermediate_sites::SiteMap;
use crate::num::PyNum;
use crate::polya;
use crate::post_filter::XTEARsltParser;
use crate::transduction::{TdAnnotation, TdSite, XTransduction};
use crate::utils::{self, split_ws};
use anyhow::{anyhow, Result};
use indexmap::IndexMap;
use std::io::Write;

/// CPython `a is b` for two equal strings coming from different places (split of a file line
/// vs. a value returned by a worker process): only the cached singletons (empty string and
/// one-character ASCII strings) are the same object.
fn py_str_is(a: &str, b: &str) -> bool {
    a == b && a.len() <= 1
}

/// Per-site statistics gathered by the read loops of the orphan/sibling functions.
#[derive(Default, Debug)]
struct SiteScan {
    m_ldisc: IndexMap<String, Vec<i64>>,
    m_rdisc: IndexMap<String, Vec<i64>>,
    n_l_rc: i64,
    n_l_nrc: i64,
    n_r_rc: i64,
    n_r_nrc: i64,
    n_lclip: i64,
    n_rclip: i64,
    n_lpolya: i64,
    n_lpolyt: i64,
    n_rpolya: i64,
    n_rpolyt: i64,
    n_total_lreads: i64,
    n_total_rreads: i64,
    n_bkgrnd_low_mapq: i64,
    l_lclip_pos: Vec<i64>,
    l_rclip_pos: Vec<i64>,
    l_local_rc_disc_pos: Vec<i64>,
    l_local_nrc_disc_pos: Vec<i64>,
}

/// The read loop shared by `parse_td_sibling_by_site_novel_sites`,
/// `filter_ori_td_by_clip_disc_position` and `check_features_for_one_site`.
/// `buggy_rclip` reproduces the use of the *first* cigar operation for the right clipped part
/// in `check_features_for_one_site`.
#[allow(clippy::too_many_arguments)]
fn scan_site(samfile: &mut AlignmentFile, names: &[String], chrm: &str, chrm_in_bam: &str, insertion_pos: i64, start_pos: i64, end_pos: i64, buggy_rclip: bool) -> Result<SiteScan> {
    let mut s = SiteScan::default();
    let n_refs = names.len() as i32;
    let i_offset = crate::py::fmt::py_int(gv::read_length() / 4.0);
    let dft_is = gv::dft_is();
    let min_pair_dist = gv::MIN_ORPHAN_DISC_PAIR_INSERT_SIZE;
    samfile.fetch(chrm_in_bam, start_pos, end_pos, |a: &Record| {
        if a.is_unmap() {
            return;
        }
        let l_cigar = a.cigar_tuples();
        if l_cigar.is_empty() {
            return;
        }
        let map_pos = a.reference_start();
        if (map_pos - insertion_pos).abs() < gv::LOCAL_COV_WIN {
            if map_pos >= insertion_pos - i_offset {
                s.n_total_rreads += 1;
            } else {
                s.n_total_lreads += 1;
            }
        }
        if a.mapping_quality() < gv::MAX_BKGRND_LOW_MAPQ && (map_pos - insertion_pos).abs() < gv::LOCAL_COV_WIN {
            s.n_bkgrnd_low_mapq += 1;
        }
        if a.mapping_quality() < gv::minimum_disc_mapq() {
            return;
        }
        if a.is_sec() || a.is_supp() {
            return;
        }
        if a.is_dup() {
            return;
        }
        let b_rc = a.is_rev();
        let query_seq = a.query_sequence().unwrap_or_default();
        let first = l_cigar[0];
        let last = l_cigar[l_cigar.len() - 1];
        let n_ck = gv::CK_POLYA_SEQ_MAX;
        if first.0 == BAM_CSOFT_CLIP && (map_pos - insertion_pos).abs() < gv::NEARBY_CLIP {
            s.n_lclip += 1;
            let clipped_seq = py_head(&query_seq, first.1 as i64);
            let s_chk = py_tail(clipped_seq, n_ck);
            let (b_at, b_a) = polya::is_consecutive_polya_t_with_ori(s_chk);
            if b_at {
                if b_a {
                    s.n_lpolya += 1;
                } else {
                    s.n_lpolyt += 1;
                }
            }
            s.l_lclip_pos.push(map_pos);
        }
        if last.0 == BAM_CSOFT_CLIP {
            let clip_pos = right_clip_pos(&l_cigar, map_pos);
            if (clip_pos - insertion_pos).abs() < gv::NEARBY_CLIP {
                s.n_rclip += 1;
                let n = if buggy_rclip { first.1 } else { last.1 } as i64;
                let clipped_seq = py_tail(&query_seq, n);
                let s_chk = py_head(clipped_seq, n_ck);
                let (b_at, b_a) = polya::is_consecutive_polya_t_with_ori(s_chk);
                if b_at {
                    if b_a {
                        s.n_rpolya += 1;
                    } else {
                        s.n_rpolyt += 1;
                    }
                }
                s.l_rclip_pos.push(clip_pos);
            }
        }
        let ntid = a.next_reference_id();
        if ntid >= 0 && ntid < n_refs && !a.mate_is_unmapped() {
            let mate_chrm = &names[ntid as usize];
            if is_decoy_contig_chrms(mate_chrm) {
                return;
            }
            let mate_pos = a.next_reference_start();
            if chrm != mate_chrm || (mate_pos - insertion_pos).abs() > min_pair_dist {
                let near = ((map_pos - insertion_pos).abs() as f64) < dft_is;
                if map_pos < insertion_pos - i_offset {
                    s.m_ldisc.entry(mate_chrm.clone()).or_default().push(mate_pos);
                    if near {
                        if b_rc {
                            s.n_l_rc += 1;
                            s.l_local_rc_disc_pos.push(map_pos);
                        } else {
                            s.n_l_nrc += 1;
                            s.l_local_nrc_disc_pos.push(map_pos);
                        }
                    }
                } else {
                    s.m_rdisc.entry(mate_chrm.clone()).or_default().push(mate_pos);
                    if near {
                        if b_rc {
                            s.n_r_rc += 1;
                            s.l_local_rc_disc_pos.push(map_pos);
                        } else {
                            s.n_r_nrc += 1;
                            s.l_local_nrc_disc_pos.push(map_pos);
                        }
                    }
                }
            }
        }
    })?;
    Ok(s)
}

/// `_form_unique_cluster` -> (is_cluster, dominant chrm, mid position, # reads of the chrm)
fn form_unique_cluster(m_disc: &IndexMap<String, Vec<i64>>, n_half_disc: f64, chrm_ratio: f64, pos_ratio: f64) -> (bool, String, i64, i64) {
    let mut n_pairs = 0i64;
    let mut n_max_chrm = 0i64;
    let mut s_max_chrm = String::new();
    let mid_pos = -1i64;
    for (mate_chrm, l) in m_disc.iter() {
        let n_temp = l.len() as i64;
        if n_temp > n_max_chrm {
            n_max_chrm = n_temp;
            s_max_chrm = mate_chrm.clone();
        }
        n_pairs += n_temp;
    }
    if (n_pairs as f64) < n_half_disc || n_pairs <= 0 {
        return (false, s_max_chrm, mid_pos, n_max_chrm);
    }
    if (n_max_chrm as f64) > n_half_disc && (n_max_chrm as f64) / (n_pairs as f64) > chrm_ratio {
        let (b, mid) = position_of_same_chrm_form_cluster(&m_disc[&s_max_chrm], pos_ratio);
        return (b, s_max_chrm, mid, n_max_chrm);
    }
    (false, s_max_chrm, mid_pos, n_max_chrm)
}

/// `_position_of_same_chrm_form_cluster`
fn position_of_same_chrm_form_cluster(l_pos: &[i64], ratio: f64) -> (bool, i64) {
    let mut l = l_pos.to_vec();
    l.sort();
    let mid_pos = l[l.len() / 2];
    let n_in = l.iter().filter(|&&p| (p - mid_pos).abs() < gv::MAX_NORMAL_INSERT_SIZE).count();
    ((n_in as f64) / (l.len() as f64) > ratio, mid_pos)
}

/// `is_disc_cluster_consist_with_clip_position`
fn is_disc_cluster_consist_with_clip_position(l_nrc_pos: &mut [i64], i_lclip_pos: i64, l_rc_pos: &mut [i64], i_rclip_pos: i64) -> bool {
    l_nrc_pos.sort();
    let mut i_mid_nrc_pos = -1;
    if !l_nrc_pos.is_empty() {
        i_mid_nrc_pos = l_nrc_pos[l_nrc_pos.len() / 4];
    }
    l_rc_pos.sort();
    let mut i_mid_rc_pos = -1;
    if !l_rc_pos.is_empty() {
        i_mid_rc_pos = l_rc_pos[crate::py::fmt::py_int(l_rc_pos.len() as f64 * 3.0 / 4.0) as usize];
    }
    if i_lclip_pos < i_mid_nrc_pos {
        return false;
    }
    !(i_rclip_pos > i_mid_rc_pos && i_mid_rc_pos > 0)
}

/// `_check_form_repeat_or_same_cluster`
fn check_form_repeat_or_same_cluster(xannotation: &TdAnnotation, src_chrm: &str, src_cluster_pos: i64, m_disc: &IndexMap<String, Vec<i64>>, n_half_disc: f64, f_rep_ratio: f64, f_cluster_ratio: f64) -> (bool, i64) {
    let mut n_in_rep = 0i64;
    let mut n_same = 0i64;
    let mut n_total = 0i64;
    let mut l_mate_cluster = Vec::new();
    for (mate_chrm, l) in m_disc.iter() {
        for &mate_pos in l {
            let (b_in_rep, _) = xannotation.is_within_repeat_region_interval_tree(mate_chrm, mate_pos);
            n_total += 1;
            if b_in_rep {
                n_in_rep += 1;
            } else if src_chrm == mate_chrm && (mate_pos - src_cluster_pos).abs() < gv::MAX_SIBLING_FLANK_LEN {
                l_mate_cluster.push(mate_pos);
                n_same += 1;
            }
        }
    }
    if (n_in_rep as f64) > n_half_disc && (n_in_rep as f64) / (n_total as f64) > f_rep_ratio {
        return (true, -1);
    }
    if (n_same as f64) > n_half_disc && (n_same as f64) / (n_total as f64) > f_cluster_ratio {
        l_mate_cluster.sort();
        let n = l_mate_cluster.len();
        return (true, l_mate_cluster[crate::py::fmt::py_int(n as f64 / 2.0) as usize]);
    }
    (false, -1)
}

/// `_is_focal_cov_normal`
fn is_focal_cov_normal(n_lread: i64, n_rreads: i64, i_window: i64, i_times: i64) -> (bool, f64, f64) {
    if i_window <= 0 {
        return (false, -1.0, -1.0);
    }
    let rl = gv::read_length();
    let n_ldepth = (n_lread as f64 * rl) / i_window as f64;
    let n_rdepth = (n_rreads as f64 * rl) / i_window as f64;
    let i_max_cov = i_times as f64 * gv::average_cov();
    (!(n_ldepth > i_max_cov || n_rdepth > i_max_cov), n_ldepth, n_rdepth)
}

/// A candidate of `parse_td_sibling_by_site_novel_sites`.
#[derive(Clone, Debug)]
struct NovelRcd {
    chrm: String,
    pos: i64,
    s_rmchrm: String,
    i_lpos: i64,
    i_rpos: i64,
    n_lclip: i64,
    n_rclip: i64,
    n_lpolya: i64,
    n_rpolyt: i64,
    nlspt: i64,
    nrspt: i64,
    n_l_rc: i64,
    n_l_nrc: i64,
    n_r_rc: i64,
    n_r_nrc: i64,
    f_lcov: f64,
    f_rcov: f64,
}

/// Merged record of `call_novel_sibling_TD_from_raw_list`.
#[derive(Clone, Debug)]
struct NovelMerged {
    src_chrm: String,
    src_lpos: i64,
    src_rpos: i64,
    n_lclip: i64,
    n_rclip: i64,
    n_lpolya: i64,
    n_rpolya: i64,
    nlspot: i64,
    nrspot: i64,
    n_l_rc: i64,
    n_l_nrc: i64,
    n_r_rc: i64,
    n_r_nrc: i64,
    f_lcov: f64,
    f_rcov: f64,
}

pub struct XOrphanTransduction {
    pub td: XTransduction,
}

impl XOrphanTransduction {
    pub fn new(working_folder: &str, n_jobs: usize, sf_reference: &str) -> Self {
        XOrphanTransduction { td: XTransduction::new(working_folder, n_jobs, sf_reference) }
    }

    fn pool(&self) -> rayon::ThreadPool {
        rayon::ThreadPoolBuilder::new().num_threads(self.td.n_jobs.max(1)).build().unwrap()
    }

    fn bams(sf_bam_list: &str) -> Result<Vec<String>> {
        let mut v = Vec::new();
        for line in utils::read_lines(sf_bam_list)? {
            let f = split_ws(&line);
            if let Some(b) = f.first() {
                v.push(b.to_string());
            }
        }
        Ok(v)
    }

    fn open_site_bam(&self, sf_bam: &str, chrm: &str) -> Result<(AlignmentFile, Vec<String>, String)> {
        let bam_info = BamInfo::new(sf_bam, &self.td.sf_reference);
        let b_with_chr = bam_info.is_chrm_contain_chr()?;
        let chrm_in_bam = process_chrm_name(chrm, b_with_chr);
        let samfile = AlignmentFile::open(sf_bam, &self.td.sf_reference)?;
        let names = samfile.references().to_vec();
        Ok((samfile, names, chrm_in_bam))
    }

    /// `parse_td_sibling_by_site` -> (chrm, pos, bl, br, src chrm, src pos, n support)
    fn parse_td_sibling_by_site(&self, chrm: &str, insertion_pos: i64, extnd: i64, sf_bam: &str, n_half_disc_cutoff: f64) -> Result<(String, i64, bool, bool, String, i64, i64)> {
        let mut start_pos = insertion_pos - extnd;
        if start_pos <= 0 {
            start_pos = 1;
        }
        let end_pos = insertion_pos + extnd;
        let (mut samfile, names, chrm_in_bam) = self.open_site_bam(sf_bam, chrm)?;
        let n_refs = names.len() as i32;
        let mut m_ldisc: IndexMap<String, Vec<i64>> = IndexMap::new();
        let mut m_rdisc: IndexMap<String, Vec<i64>> = IndexMap::new();
        let mut l_nrc: Vec<i64> = Vec::new();
        let mut l_rc: Vec<i64> = Vec::new();
        let max_is = gv::MAX_NORMAL_INSERT_SIZE;
        samfile.fetch(&chrm_in_bam, start_pos, end_pos, |a: &Record| {
            if a.is_sec() || a.is_supp() || a.is_dup() || a.is_unmap() {
                return;
            }
            if a.cigar_tuples().is_empty() {
                return;
            }
            if a.mapping_quality() < gv::minimum_disc_mapq() {
                return;
            }
            let map_pos = a.reference_start();
            let is_rc = a.is_rev();
            let ntid = a.next_reference_id();
            if ntid >= 0 && ntid < n_refs && !a.mate_is_unmapped() {
                let mate_chrm = &names[ntid as usize];
                if is_decoy_contig_chrms(mate_chrm) {
                    return;
                }
                let mate_pos = a.next_reference_start();
                if chrm != mate_chrm || (mate_pos - insertion_pos).abs() > max_is {
                    if map_pos < insertion_pos {
                        m_ldisc.entry(mate_chrm.clone()).or_default().push(mate_pos);
                    } else {
                        m_rdisc.entry(mate_chrm.clone()).or_default().push(mate_pos);
                    }
                    if !is_rc {
                        l_nrc.push(map_pos);
                    } else {
                        l_rc.push(map_pos);
                    }
                }
            }
        })?;
        let r = gv::TRANSDCT_MULTI_SOURCE_MIN_RATIO;
        let (b_l, s_lm, i_lpos, nlspt) = form_unique_cluster(&m_ldisc, n_half_disc_cutoff, r, r);
        let (b_r, s_rm, i_rpos, nrspt) = form_unique_cluster(&m_rdisc, n_half_disc_cutoff, r, r);
        let consist = is_disc_cluster_consist_with_clip_position(&mut l_nrc, insertion_pos, &mut l_rc, insertion_pos);
        if b_l && !b_r && consist {
            return Ok((chrm.to_string(), insertion_pos, true, false, s_lm, i_lpos, nlspt));
        } else if !b_l && b_r && consist {
            return Ok((chrm.to_string(), insertion_pos, false, true, s_rm, i_rpos, nrspt));
        }
        Ok((chrm.to_string(), insertion_pos, false, false, String::new(), -1, 0))
    }

    /// `call_sibling_TD_from_existing_list`
    #[allow(clippy::too_many_arguments)]
    pub fn call_sibling_td_from_existing_list(
        &self,
        sf_cns_sites: &str,
        sf_bam_list: &str,
        extnd: i64,
        n_half_disc_cutoff: f64,
        i_search_win: i64,
        xannotation: &TdAnnotation,
        i_rep_type: i64,
        i_max_cov: f64,
        sf_updated_cns: &str,
        sf_sibling_td: &str,
    ) -> Result<()> {
        let mut m_merged: SiteMap<TdSite> = IndexMap::new();
        for sf_bam in Self::bams(sf_bam_list)? {
            let mut l_sites: Vec<(String, i64)> = Vec::new();
            for line in utils::read_lines(sf_cns_sites)? {
                let f = split_ws(utils::rstrip(&line));
                if f.len() <= 1 {
                    continue;
                }
                l_sites.push((f[0].to_string(), utils::pint(f[1])?));
            }
            let rslts: Vec<Result<(String, i64, bool, bool, String, i64, i64)>> =
                self.pool().install(|| utils::par_map(l_sites.len(), |i| self.parse_td_sibling_by_site(&l_sites[i].0, l_sites[i].1, extnd, &sf_bam, 0.0)));
            for r in rslts {
                let (tmp_chrm, tmp_pos, bl, br, s_src_chrm, i_src_pos, nspot) = r?;
                if !bl && !br {
                    continue;
                }
                if bl && br {
                    continue;
                }
                let e = m_merged.entry(tmp_chrm).or_default();
                match e.get(&tmp_pos).cloned() {
                    None => {
                        e.insert(tmp_pos, TdSite::Sibling { bl, br, src_chrm: s_src_chrm, src_pos: i_src_pos, nspot });
                    }
                    Some(TdSite::Sibling { bl: obl, br: obr, src_chrm: osrc, src_pos: opos, nspot: ons }) => {
                        if bl != obl || br != obr || s_src_chrm != osrc {
                            continue;
                        }
                        if (i_src_pos - opos).abs() > gv::MAX_NORMAL_INSERT_SIZE {
                            continue;
                        }
                        let new_pos = if ons > nspot { opos } else { i_src_pos };
                        e.insert(tmp_pos, TdSite::Sibling { bl, br, src_chrm: s_src_chrm, src_pos: new_pos, nspot: nspot + ons });
                    }
                    Some(_) => unreachable!(),
                }
            }
        }
        // the "one-side" related candidates
        let mut m_one_side: SiteMap<f64> = IndexMap::new();
        for rcd in XTEARsltParser::load_in_xtea_rslt_ori(sf_cns_sites)? {
            if XTEARsltParser::is_sibling_transduction_related_ins(&rcd)? {
                let p = utils::pint(&rcd[1])?;
                m_one_side.entry(rcd[0].clone()).or_default().insert(p, XTEARsltParser::get_two_side_coverage_diff(&rcd)?);
            }
        }
        let empty: SiteMap<[i64; 6]> = IndexMap::new();
        let (m_new_ins, _m_depth) = self.td.filter_transduction(&m_merged, xannotation, i_rep_type, sf_bam_list, i_max_cov, &empty, 0.0)?;
        self.td.update_existing_results(sf_cns_sites, &m_new_ins, self.td.same_sites_extnd, sf_updated_cns, None)?;
        let mut fout = utils::create(sf_sibling_td)?;
        for (tmp_chrm, m) in m_new_ins.iter() {
            for (&tmp_pos, rcd) in m.iter() {
                let (src_chrm, src_pos, n_sprt) = match rcd {
                    TdSite::Sibling { src_chrm, src_pos, nspot, .. } => (src_chrm.clone(), *src_pos, *nspot),
                    _ => unreachable!(),
                };
                if (n_sprt as f64) < n_half_disc_cutoff {
                    continue;
                }
                // `_is_sibling_within_candidates` (NOTE: `is not` comparison of chromosome names)
                let mut hit: Option<(String, i64, f64)> = None;
                'outer: for (ins_chrm, mp) in m_one_side.iter() {
                    if !py_str_is(ins_chrm, &src_chrm) {
                        continue;
                    }
                    for (&ins_pos, &f_cov_diff) in mp.iter() {
                        if (ins_pos - src_pos).abs() < i_search_win {
                            hit = Some((ins_chrm.clone(), ins_pos, f_cov_diff));
                            break 'outer;
                        }
                    }
                }
                let sinfo = match hit {
                    None => format!("{}\t{}\t{}\t{}\n", tmp_chrm, tmp_pos, src_chrm, src_pos),
                    Some((hc, hp, hcov)) => {
                        let mut s = format!("{}\t{}\t{}\t{}\n", tmp_chrm, tmp_pos, hc, hp);
                        if let Some(&own) = m_one_side.get(tmp_chrm).and_then(|x| x.get(&tmp_pos)) {
                            if own > hcov {
                                s = format!("{}\t{}\t{}\t{}\n", hc, hp, tmp_chrm, tmp_pos);
                            }
                        }
                        s
                    }
                };
                fout.write_all(sinfo.as_bytes())?;
            }
        }
        fout.flush()?;
        Ok(())
    }

    /// `parse_td_sibling_by_site_novel_sites` for one chromosome
    #[allow(clippy::too_many_arguments)]
    fn parse_td_sibling_by_site_novel_sites(&self, chrm: &str, extnd: i64, l_positions: &[i64], sf_bam: &str, xannotation: &TdAnnotation, n_polya_cutoff: f64, n_half_disc_cutoff: f64) -> Result<Vec<NovelRcd>> {
        let (mut samfile, names, chrm_in_bam) = self.open_site_bam(sf_bam, chrm)?;
        let mut l_rslts = Vec::new();
        for &insertion_pos in l_positions {
            let mut start_pos = insertion_pos - extnd;
            if start_pos <= 0 {
                start_pos = 1;
            }
            let end_pos = insertion_pos + extnd;
            let mut s = scan_site(&mut samfile, &names, chrm, &chrm_in_bam, insertion_pos, start_pos, end_pos, false)?;
            let r_multi = gv::TRANSDCT_MULTI_SOURCE_MIN_RATIO;
            let r_region = gv::TRANSDCT_REGION_CLUSTER_MIN_RATIO;
            let (b_lclstr, s_lmchrm, mut i_lpos, nlspt) = form_unique_cluster(&s.m_ldisc, n_half_disc_cutoff, r_multi, r_region);
            let (b_rclstr, s_rmchrm, mut i_rpos, nrspt) = form_unique_cluster(&s.m_rdisc, n_half_disc_cutoff, r_multi, r_region);
            let mut b_valid_polya = false;
            if (s.n_lpolya as f64) >= n_polya_cutoff || (s.n_rpolyt as f64) >= n_polya_cutoff {
                if (s.n_lclip > 0 && (s.n_lpolya as f64) / (s.n_lclip as f64) > gv::MIN_POLYA_CLIP_RATIO)
                    || (s.n_rclip > 0 && (s.n_rpolyt as f64) / (s.n_rclip as f64) > gv::MIN_POLYA_CLIP_RATIO)
                {
                    b_valid_polya = true;
                }
            }
            let b_disc_clip_pos_consist =
                is_disc_cluster_consist_with_clip_position(&mut s.l_local_nrc_disc_pos, insertion_pos, &mut s.l_local_rc_disc_pos, insertion_pos);
            let n_total_focal_ldisc = s.n_l_rc + s.n_l_nrc;
            let b_l_dom = n_total_focal_ldisc > 0 && (s.n_l_nrc as f64) / (n_total_focal_ldisc as f64) > gv::MIN_CLUSTER_RC_RATIO;
            let n_total_focal_rdisc = s.n_r_rc + s.n_r_nrc;
            let b_r_dom = n_total_focal_rdisc > 0 && (s.n_r_rc as f64) / (n_total_focal_rdisc as f64) > gv::MIN_CLUSTER_RC_RATIO;
            let mut b_consistent = false;
            if b_lclstr && b_rclstr {
                if s_lmchrm == s_rmchrm && (i_lpos - i_rpos).abs() < gv::MAX_SIBLING_FLANK_LEN && b_l_dom && b_r_dom {
                    b_consistent = true;
                }
            } else if b_lclstr && !b_rclstr {
                if b_l_dom && b_r_dom {
                    let (b_rcluster, rp) = check_form_repeat_or_same_cluster(
                        xannotation,
                        &s_lmchrm,
                        i_lpos,
                        &s.m_rdisc,
                        n_half_disc_cutoff,
                        gv::TRANSDCT_MULTI_SOURCE_MIN_RATIO,
                        gv::TRANSDCT_OTHER_SIDE_CLUSTER_MIN_RATIO,
                    );
                    i_rpos = rp;
                    if b_rcluster {
                        b_consistent = true;
                    }
                }
            } else if !b_lclstr && b_rclstr && b_l_dom && b_r_dom {
                let (b_lcluster, lp) = check_form_repeat_or_same_cluster(
                    xannotation,
                    &s_rmchrm,
                    i_rpos,
                    &s.m_ldisc,
                    n_half_disc_cutoff,
                    gv::TRANSDCT_MULTI_SOURCE_MIN_RATIO,
                    gv::TRANSDCT_OTHER_SIDE_CLUSTER_MIN_RATIO,
                );
                i_lpos = lp;
                if b_lcluster {
                    b_consistent = true;
                }
            }
            let (b_cov_normal, f_lcov, f_rcov) = is_focal_cov_normal(s.n_total_lreads, s.n_total_rreads, gv::LOCAL_COV_WIN, gv::MAX_COV_TIMES);
            let n_tot = s.n_total_lreads + s.n_total_rreads;
            let b_pass_bkgrnd_chk = !(n_tot <= 0 || (s.n_bkgrnd_low_mapq as f64) / (n_tot as f64) > gv::MAX_BKGRND_LOW_MAPQ_RATIO);
            if b_consistent && b_valid_polya && b_disc_clip_pos_consist && b_cov_normal && b_pass_bkgrnd_chk {
                l_rslts.push(NovelRcd {
                    chrm: chrm.to_string(),
                    pos: insertion_pos,
                    s_rmchrm: s_rmchrm.clone(),
                    i_lpos,
                    i_rpos,
                    n_lclip: s.n_lclip,
                    n_rclip: s.n_rclip,
                    n_lpolya: s.n_lpolya,
                    n_rpolyt: s.n_rpolyt,
                    nlspt,
                    nrspt,
                    n_l_rc: s.n_l_rc,
                    n_l_nrc: s.n_l_nrc,
                    n_r_rc: s.n_r_rc,
                    n_r_nrc: s.n_r_nrc,
                    f_lcov,
                    f_rcov,
                });
            }
        }
        Ok(l_rslts)
    }

    /// `_cnt_low_mapq_reads`
    fn cnt_low_mapq_reads(&self, sf_bam: &str, chrm: &str, start_pos: i64, end_pos: i64) -> Result<(i64, i64)> {
        let (mut samfile, _names, chrm_in_bam) = self.open_site_bam(sf_bam, chrm)?;
        let (mut n_total, mut n_low) = (0i64, 0i64);
        samfile.fetch(&chrm_in_bam, start_pos, end_pos, |a: &Record| {
            if a.is_unmap() || a.cigar_tuples().is_empty() {
                return;
            }
            n_total += 1;
            if a.mapping_quality() < gv::MAX_BKGRND_LOW_MAPQ {
                n_low += 1;
            }
        })?;
        Ok((n_total, n_low))
    }

    /// `is_qualified_source_region`
    fn is_qualified_source_region(&self, rcd: &NovelRcd, sf_bam: &str) -> Result<Option<NovelRcd>> {
        let (mut start_pos, mut end_pos) = (rcd.i_lpos, rcd.i_rpos);
        if rcd.i_lpos > rcd.i_rpos {
            start_pos = rcd.i_rpos;
            end_pos = rcd.i_lpos;
        }
        if start_pos <= 0 && end_pos <= 0 {
            return Ok(None);
        }
        if start_pos < 0 {
            start_pos = end_pos - 500;
        }
        let (n_total, n_low) = self.cnt_low_mapq_reads(sf_bam, &rcd.s_rmchrm, start_pos, end_pos)?;
        if n_total > 0 && (n_low as f64) / (n_total as f64) < gv::MAX_BKGRND_LOW_MAPQ_RATIO {
            return Ok(Some(rcd.clone()));
        }
        Ok(None)
    }

    /// Write the black-list filtered sites of `sf_sites` to `sf_tmp` -> chromosomes in order
    fn filter_sites_by_blacklist(sf_sites: &str, sf_tmp: &str, x_blklist: &XBlackList) -> Result<Vec<String>> {
        let mut m_chrms: IndexMap<String, ()> = IndexMap::new();
        let mut fout = utils::create(sf_tmp)?;
        for line in utils::read_lines(sf_sites)? {
            let f = split_ws(utils::rstrip(&line));
            if f.len() <= 1 {
                continue;
            }
            let chrm = f[0];
            let pos = utils::pint(f[1])?;
            if x_blklist.fall_in_region(chrm, pos).0 {
                println!("{}:{} fall in black list region, filtered out!", chrm, pos);
                continue;
            }
            m_chrms.insert(chrm.to_string(), ());
            writeln!(fout, "{}\t{}", chrm, pos)?;
        }
        fout.flush()?;
        Ok(m_chrms.keys().cloned().collect())
    }

    fn positions_of_chrm(sf_sites: &str, chrm: &str) -> Result<Vec<i64>> {
        let mut v = Vec::new();
        for line in utils::read_lines(sf_sites)? {
            let f = split_ws(&line);
            if f.len() < 2 {
                return Err(anyhow!("IndexError: malformed line in {}", sf_sites));
            }
            let p = utils::pint(f[1])?;
            if f[0] == chrm {
                v.push(p);
            }
        }
        Ok(v)
    }

    /// `call_novel_sibling_TD_from_raw_list`
    #[allow(clippy::too_many_arguments)]
    pub fn call_novel_sibling_td_from_raw_list(&self, sf_raw_sites: &str, sf_bam_list: &str, extnd: i64, n_clip_cutoff: i64, n_disc_cutoff: i64, sf_black_list: &str, sf_rmsk: &str, sf_sibling_td: &str) -> Result<()> {
        let mut x_blklist = XBlackList::new();
        x_blklist.load_index_regions(sf_black_list)?;
        // every worker builds the same annotation in the original; build it once
        let xannotation = self.td.prep_annotation_interval_tree(sf_rmsk, 225)?;
        let mut m_merged: SiteMap<NovelMerged> = IndexMap::new();
        for sf_bam in Self::bams(sf_bam_list)? {
            let sf_tmp_file = format!("{}.tmp", sf_raw_sites);
            let chrms = Self::filter_sites_by_blacklist(sf_raw_sites, &sf_tmp_file, &x_blklist)?;
            let n_tmp_polya_cutoff = n_clip_cutoff as f64 / 2.0;
            let n_half_disc_cutoff = n_disc_cutoff as f64 / 2.0;
            let rslts: Vec<Result<Vec<NovelRcd>>> = self.pool().install(|| {
                utils::par_map(chrms.len(), |i| {
                    let l_pos = Self::positions_of_chrm(&sf_tmp_file, &chrms[i])?;
                    self.parse_td_sibling_by_site_novel_sites(&chrms[i], extnd, &l_pos, &sf_bam, &xannotation, n_tmp_polya_cutoff, n_half_disc_cutoff)
                })
            });
            let mut l_tmp: Vec<NovelRcd> = Vec::new();
            for r in rslts {
                l_tmp.extend(r?);
            }
            let checked: Vec<Result<Option<NovelRcd>>> = self.pool().install(|| utils::par_map(l_tmp.len(), |i| self.is_qualified_source_region(&l_tmp[i], &sf_bam)));
            for r in checked {
                let rcd = match r? {
                    None => continue,
                    Some(x) => x,
                };
                let e = m_merged.entry(rcd.chrm.clone()).or_default();
                match e.get(&rcd.pos).cloned() {
                    None => {
                        e.insert(
                            rcd.pos,
                            NovelMerged {
                                src_chrm: rcd.s_rmchrm.clone(),
                                src_lpos: rcd.i_lpos,
                                src_rpos: rcd.i_rpos,
                                n_lclip: rcd.n_lclip,
                                n_rclip: rcd.n_rclip,
                                n_lpolya: rcd.n_lpolya,
                                n_rpolya: rcd.n_rpolyt,
                                nlspot: rcd.nlspt,
                                nrspot: rcd.nrspt,
                                n_l_rc: rcd.n_l_rc,
                                n_l_nrc: rcd.n_l_nrc,
                                n_r_rc: rcd.n_r_rc,
                                n_r_nrc: rcd.n_r_nrc,
                                f_lcov: rcd.f_lcov,
                                f_rcov: rcd.f_rcov,
                            },
                        );
                    }
                    Some(t) => {
                        if t.src_chrm != rcd.s_rmchrm {
                            continue;
                        }
                        if (rcd.i_lpos - t.src_lpos).abs() > gv::MAX_NORMAL_INSERT_SIZE {
                            continue;
                        }
                        if (rcd.i_rpos - t.src_rpos).abs() > gv::MAX_NORMAL_INSERT_SIZE {
                            continue;
                        }
                        e.insert(
                            rcd.pos,
                            NovelMerged {
                                src_chrm: rcd.s_rmchrm.clone(),
                                src_lpos: rcd.i_lpos,
                                src_rpos: rcd.i_rpos,
                                n_lclip: rcd.n_lclip + t.n_lclip,
                                n_rclip: rcd.n_rclip + t.n_rclip,
                                n_lpolya: rcd.n_lpolya + t.n_lpolya,
                                n_rpolya: rcd.n_rpolyt + t.n_rpolya,
                                nlspot: rcd.nlspt + t.nlspot,
                                nrspot: rcd.nrspt + t.nrspot,
                                n_l_rc: rcd.n_l_rc + t.n_l_rc,
                                n_l_nrc: rcd.n_l_nrc + t.n_l_nrc,
                                n_r_rc: rcd.n_r_rc + t.n_r_rc,
                                n_r_nrc: rcd.n_r_nrc + t.n_r_nrc,
                                f_lcov: rcd.f_lcov + t.f_lcov,
                                f_rcov: rcd.f_rcov + t.f_rcov,
                            },
                        );
                    }
                }
            }
        }
        let mut m_site_info: SiteMap<TdSite> = IndexMap::new();
        let sf_new_sites = format!("{}.tmp", sf_sibling_td);
        {
            let mut fout = utils::create(&sf_new_sites)?;
            let half_clip = n_clip_cutoff as f64 / 2.0;
            let half_disc = n_disc_cutoff as f64 / 2.0;
            for (tmp_chrm, m) in m_merged.iter() {
                for (&tmp_pos, t) in m.iter() {
                    if (t.n_lclip as f64) < half_clip && (t.n_rclip as f64) < half_clip {
                        println!("{}:{} is filtered out, as no enough left or right clipped reads!", tmp_chrm, tmp_pos);
                        continue;
                    }
                    if (t.n_lpolya as f64) < half_clip && (t.n_rpolya as f64) < half_clip {
                        println!("{}:{} is filtered out, as no enough polyA(T) reads!", tmp_chrm, tmp_pos);
                        continue;
                    }
                    if (t.nlspot as f64) < half_disc || (t.nrspot as f64) < half_disc {
                        println!("{}:{} is filtered out, as no enough left ({}) and right ({}) disc reads!", tmp_chrm, tmp_pos, t.nlspot, t.nrspot);
                        continue;
                    }
                    let (mut i_start, mut i_end) = (t.src_lpos, t.src_rpos);
                    if t.src_lpos > t.src_rpos {
                        i_start = t.src_rpos;
                        i_end = t.src_lpos;
                    }
                    writeln!(
                        fout,
                        "{}\t{}\t{}:{}-{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        tmp_chrm, tmp_pos, t.src_chrm, i_start, i_end, t.n_lclip, t.n_rclip, t.n_lpolya, t.n_rpolya, t.nlspot, t.nrspot, t.n_l_rc, t.n_l_nrc, t.n_r_rc, t.n_r_nrc
                    )?;
                    let s_src = format!("{}:{}-{}", t.src_chrm, i_start, i_end);
                    m_site_info.entry(tmp_chrm.clone()).or_default().insert(
                        tmp_pos,
                        TdSite::Info {
                            n_l_disc: t.nlspot,
                            lcns_start: PyNum::Int(i_start),
                            lcns_end: PyNum::Int(i_start),
                            n_r_disc: t.nrspot,
                            rcns_start: PyNum::Int(i_end),
                            rcns_end: PyNum::Int(i_end),
                            n_cns_clip: t.n_lclip + t.n_rclip,
                            n_polya_clip: t.n_lpolya + t.n_rpolya,
                            lcov: PyNum::Float(t.f_lcov),
                            rcov: PyNum::Float(t.f_rcov),
                            s_src,
                        },
                    );
                }
            }
            fout.flush()?;
        }
        let i_max_cov = gv::MAX_COV_TIMES as f64 * gv::average_cov();
        let m_gntp = self.td.collect_gntp_features(&sf_new_sites, sf_bam_list, i_max_cov)?;
        self.td.add_new_rcd_with_cutoff(&m_site_info, &m_gntp, gv::average_cov(), sf_sibling_td, gv::ORPHAN_TRANSDUCTION)
    }

    /// `filter_ori_td_by_clip_disc_position` for one chromosome -> failed sites
    fn filter_ori_td_by_clip_disc_position(&self, chrm: &str, extnd: i64, l_positions: &[i64], sf_bam: &str) -> Result<Vec<(String, i64)>> {
        let (mut samfile, names, chrm_in_bam) = self.open_site_bam(sf_bam, chrm)?;
        let mut l_failed = Vec::new();
        for &insertion_pos in l_positions {
            let mut start_pos = insertion_pos - extnd;
            if start_pos <= 0 {
                start_pos = 1;
            }
            let end_pos = insertion_pos + extnd;
            let mut s = scan_site(&mut samfile, &names, chrm, &chrm_in_bam, insertion_pos, start_pos, end_pos, false)?;
            if !is_disc_cluster_consist_with_clip_position(&mut s.l_local_nrc_disc_pos, insertion_pos, &mut s.l_local_rc_disc_pos, insertion_pos) {
                l_failed.push((chrm.to_string(), insertion_pos));
            }
        }
        Ok(l_failed)
    }

    /// `distinguish_source_from_insertion_for_td` -> failed sites
    #[allow(clippy::too_many_arguments)]
    pub fn distinguish_source_from_insertion_for_td(&self, sf_ori_td_sites: &str, sf_bam_list: &str, extnd: i64, sf_black_list: &str, sf_rmsk: &str) -> Result<SiteMap<()>> {
        let mut x_blklist = XBlackList::new();
        x_blklist.load_index_regions(sf_black_list)?;
        // the original builds (but does not use) the annotation in every worker; parse it once
        // so that malformed files fail the same way
        let _ = self.td.prep_annotation_interval_tree(sf_rmsk, 225)?;
        let mut m_failed: SiteMap<()> = IndexMap::new();
        for sf_bam in Self::bams(sf_bam_list)? {
            let sf_tmp_file = format!("{}.tmp", sf_ori_td_sites);
            let chrms = Self::filter_sites_by_blacklist(sf_ori_td_sites, &sf_tmp_file, &x_blklist)?;
            let rslts: Vec<Result<Vec<(String, i64)>>> = self.pool().install(|| {
                utils::par_map(chrms.len(), |i| {
                    let l_pos = Self::positions_of_chrm(&sf_tmp_file, &chrms[i])?;
                    self.filter_ori_td_by_clip_disc_position(&chrms[i], extnd, &l_pos, &sf_bam)
                })
            });
            for r in rslts {
                for (c, p) in r? {
                    m_failed.entry(c).or_default().insert(p, ());
                }
            }
        }
        Ok(m_failed)
    }

    /// `update_existing_list_only`
    pub fn update_existing_list_only(sf_sites: &str, m_failed_td: &SiteMap<()>) -> Result<()> {
        let m_existing = Self::load_existing(sf_sites, m_failed_td)?;
        let mut fout = utils::create(sf_sites)?;
        for (_, m) in m_existing.iter() {
            for (_, l) in m.iter() {
                writeln!(fout, "{}", l)?;
            }
        }
        fout.flush()?;
        Ok(())
    }

    fn load_existing(sf: &str, m_failed_td: &SiteMap<()>) -> Result<SiteMap<String>> {
        let mut m: SiteMap<String> = IndexMap::new();
        for line in utils::read_raw_lines(sf)? {
            let f = split_ws(&line);
            if f.len() < 2 {
                return Err(anyhow!("IndexError: malformed line in {}", sf));
            }
            let chrm = f[0];
            let pos = utils::pint(f[1])?;
            if m_failed_td.get(chrm).map_or(false, |x| x.contains_key(&pos)) {
                continue;
            }
            m.entry(chrm.to_string()).or_default().insert(pos, utils::rstrip(&line).to_string());
        }
        Ok(m)
    }

    /// `append_to_existing_list`
    pub fn append_to_existing_list(sf_new_sites: &str, sf_existing: &str, m_failed_td: &SiteMap<()>) -> Result<()> {
        let m_existing = Self::load_existing(sf_existing, m_failed_td)?;
        {
            let mut fout = utils::create(sf_existing)?;
            for (_, m) in m_existing.iter() {
                for (_, l) in m.iter() {
                    writeln!(fout, "{}", l)?;
                }
            }
            fout.flush()?;
        }
        let mut fapp = utils::append(sf_existing)?;
        for line in utils::read_raw_lines(sf_new_sites)? {
            let f = split_ws(&line);
            if f.len() < 2 {
                return Err(anyhow!("IndexError: malformed line in {}", sf_new_sites));
            }
            let pos = utils::pint(f[1])?;
            if m_existing.get(f[0]).map_or(false, |x| x.contains_key(&pos)) {
                continue;
            }
            writeln!(fapp, "{}", utils::rstrip(&line))?;
        }
        fapp.flush()?;
        Ok(())
    }

    /// `check_features_for_one_site` -> (chrm, pos, passed)
    fn check_features_for_one_site(&self, chrm: &str, pos: i64, sf_bam: &str, n_clip_cutoff: i64, n_disc_cutoff: i64, s_src: &str) -> Result<(String, i64, bool)> {
        let fail = || Ok((chrm.to_string(), pos, false));
        let mut region_chrm = String::new();
        let mut region_start: i64 = 0;
        let mut region_end: i64 = 0;
        if s_src.contains(gv::SIBLING_LABEL) {
            if s_src.contains(':') {
                let f: Vec<&str> = s_src.split(':').collect();
                if f.len() < 2 {
                    return fail();
                }
                region_chrm = f[0].to_string();
                let r: Vec<&str> = f[1].split('~').collect();
                let v = utils::pint(r[0])?;
                region_start = v - 100;
                region_end = v + 100;
            } else if s_src.contains('_') {
                let f: Vec<&str> = s_src.split('_').collect();
                if f.len() < 2 {
                    return fail();
                }
                region_chrm = f[0].to_string();
                let v = utils::pint(f[1])?;
                region_start = v - 100;
                region_end = v + 100;
            }
        } else {
            let f: Vec<&str> = s_src.split(':').collect();
            if f.len() < 2 {
                return fail();
            }
            let mut s_region = f[1];
            if s_region.starts_with('-') {
                s_region = &f[1][1..];
            }
            let r: Vec<&str> = s_region.split('-').collect();
            if r.len() < 2 {
                return fail();
            }
            region_chrm = f[0].to_string();
            region_start = utils::pint(r[0])?;
            if r[1].contains('~') {
                region_end = utils::pint(r[1].split('~').next().unwrap())?;
            } else {
                region_end = utils::pint(r[r.len() - 1])?;
            }
            if region_start.abs() == 1 {
                region_start = region_end - 100;
            }
            if region_end.abs() == 1 {
                region_end = region_start + 100;
            }
        }
        let (mut samfile, names, chrm_in_bam) = self.open_site_bam(sf_bam, chrm)?;
        let start_pos = pos - gv::MAX_NORMAL_INSERT_SIZE;
        let end_pos = pos + gv::MAX_NORMAL_INSERT_SIZE;
        let s = scan_site(&mut samfile, &names, chrm, &chrm_in_bam, pos, start_pos, end_pos, true)?;
        region_start -= gv::MAX_NORMAL_INSERT_SIZE;
        region_end += gv::MAX_NORMAL_INSERT_SIZE;
        let cnt = |m: &IndexMap<String, Vec<i64>>| -> i64 {
            match m.get(&region_chrm) {
                None => 0,
                Some(l) => l.iter().filter(|&&p| p >= region_start && p <= region_end).count() as i64,
            }
        };
        if (cnt(&s.m_ldisc) as f64) >= n_disc_cutoff as f64 / 2.0 {
            return fail();
        }
        if (cnt(&s.m_rdisc) as f64) >= n_disc_cutoff as f64 / 2.0 {
            return fail();
        }
        if (s.n_lpolya as f64) > n_clip_cutoff as f64 / 2.0 || (s.n_rpolyt as f64) > n_clip_cutoff as f64 / 2.0 {
            return fail();
        }
        Ok((chrm.to_string(), pos, true))
    }

    /// `check_features_for_given_sites` -> sites to filter out
    pub fn check_features_for_given_sites(&self, sf_raw_sites: &str, sf_bam_list: &str, n_clip_cutoff: i64, n_disc_cutoff: i64) -> Result<SiteMap<()>> {
        let mut m_out: SiteMap<()> = IndexMap::new();
        for sf_bam in Self::bams(sf_bam_list)? {
            let mut l_rcd: Vec<(String, i64, String)> = Vec::new();
            for line in utils::read_lines(sf_raw_sites)? {
                let f: Vec<String> = split_ws(utils::rstrip(&line)).into_iter().map(|s| s.to_string()).collect();
                if f.len() <= 1 {
                    continue;
                }
                let p = utils::pint(&f[1])?;
                let s_src = XTEARsltParser::get_transduction_source(&f)?;
                l_rcd.push((f[0].clone(), p, s_src));
            }
            let rslts: Vec<Result<(String, i64, bool)>> = self
                .pool()
                .install(|| utils::par_map(l_rcd.len(), |i| self.check_features_for_one_site(&l_rcd[i].0, l_rcd[i].1, &sf_bam, n_clip_cutoff, n_disc_cutoff, &l_rcd[i].2)));
            for r in rslts {
                let (c, p, ok) = r?;
                if !ok {
                    m_out.entry(c).or_default().insert(p, ());
                }
            }
        }
        Ok(m_out)
    }
}
