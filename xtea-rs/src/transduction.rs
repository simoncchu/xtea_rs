//! Port of `x_transduction.py` (`XTransduction`): calling 3' transductions by re-aligning the
//! clipped and discordant reads to the flanking regions of full-length reference (and
//! polymorphic) copies.

use crate::bam::{AlignmentFile, PyAln, Record};
use crate::bwa_align::BwAlign;
use crate::clip_disc_filter::{is_clipped_part_qualified_algnmt, py_head, py_tail, DiscFlt, XClipDiscFilter};
use crate::cluster::is_disc_cluster;
use crate::coverage::ReadDepth;
use crate::fasta::FastaFile;
use crate::genotype_feature::{GntpFeat, XGenotyper};
use crate::globals as gv;
use crate::intermediate_sites::{self as xsites, SiteMap};
use crate::intervaltree::IntervalTree;
use crate::num::PyNum;
use crate::post_filter::XTEARsltParser;
use crate::utils::{self, split_ws};
use anyhow::{anyhow, Result};
use indexmap::IndexMap;
use std::collections::HashMap;
use std::io::Write;

// ------------------------------------------------------------------------------------------------
// Annotation used by the transduction modules (`prep_annotation_interval_tree`)

/// One record of `load_rmsk_annotation_with_extnd_div_with_lenth_cutoff`.
#[derive(Clone, Debug)]
struct TdRec {
    end: f64,
    div: f64,
    sub_type: String,
    family: String,
}

/// The RepeatMasker annotation loaded with `load_rmsk_annotation_with_extnd_div_with_lenth_cutoff`
/// and indexed with interval trees.  The extension may be a float (the mean insert size in
/// the "sibling" step), in which case the interval boundaries are floats, exactly as in Python.
pub struct TdAnnotation {
    b_with_chr: bool,
    m_rmsk: IndexMap<String, IndexMap<u64, (f64, Vec<TdRec>)>>,
    m_tree: HashMap<String, IntervalTree>,
}

#[inline]
fn fk(v: f64) -> u64 {
    if v == 0.0 {
        0
    } else {
        v.to_bits()
    }
}

/// Result of `get_div_subfamily_with_min_div`: (div_rate, sub_family, family, pos, pos_end);
/// `div_rate` is the int -1 and the family `None` when absent.
#[derive(Clone, Debug)]
pub struct DivInfo {
    pub div_rate: PyNum,
    pub sub_family: String,
    pub family: Option<String>,
    pub pos: f64,
    pub pos_end: i64,
}

impl TdAnnotation {
    /// `prep_annotation_interval_tree`
    pub fn prep(sf_rmsk: &str, i_extnd: f64, i_min_copy_len: i64, b_with_chr: bool) -> Result<TdAnnotation> {
        let mut a = TdAnnotation { b_with_chr, m_rmsk: IndexMap::new(), m_tree: HashMap::new() };
        for line in utils::read_lines(sf_rmsk)? {
            if utils::rstrip(&line).is_empty() {
                continue;
            }
            let fields = split_ws(&line);
            if fields.len() < 15 {
                continue;
            }
            if fields[0] == "SW" || fields[0] == "score" {
                continue;
            }
            let mismatch_rate = crate::py::fmt::parse_py_float(fields[1]).ok_or_else(|| anyhow!("ValueError: {}", fields[1]))?;
            let chrm = crate::annotation::process_chrm_name(fields[4], b_with_chr);
            let start_pos = utils::pint(fields[5])?;
            let end_pos = utils::pint(fields[6])?;
            let b_rc = fields[8] == "C";
            let sub_type = fields[9].to_string();
            let family = fields[10].to_string();
            let _ = if b_rc { (utils::pint(fields[13])?, utils::pint(fields[12])?) } else { (utils::pint(fields[11])?, utils::pint(fields[12])?) };
            if (end_pos - start_pos).abs() < i_min_copy_len {
                continue;
            }
            let m = a.m_rmsk.entry(chrm.clone()).or_default();
            if m.contains_key(&fk(start_pos as f64)) {
                println!("Position {}:{} has more than  1 annotation!!!!", chrm, start_pos);
            }
            let extd_start = start_pos as f64 - i_extnd;
            let extd_end = end_pos as f64 + i_extnd;
            m.entry(fk(extd_start)).or_insert_with(|| (extd_start, Vec::new())).1.push(TdRec { end: extd_end, div: mismatch_rate, sub_type, family });
        }
        for (chrm, m) in a.m_rmsk.iter() {
            let mut t = IntervalTree::new();
            for (_, (pos, recs)) in m.iter() {
                t.addf(*pos, recs[0].end);
            }
            a.m_tree.insert(chrm.clone(), t);
        }
        Ok(a)
    }

    fn chrm(&self, c: &str) -> String {
        crate::annotation::process_chrm_name(c, self.b_with_chr)
    }

    /// `is_within_repeat_region_interval_tree` -> (hit, start of the first hit)
    pub fn is_within_repeat_region_interval_tree(&self, chrm1: &str, pos: i64) -> (bool, f64) {
        match self.m_tree.get(&self.chrm(chrm1)) {
            None => (false, -1.0),
            Some(t) => match t.at_f(pos as f64).first() {
                Some((b, _)) => (true, *b),
                None => (false, -1.0),
            },
        }
    }

    /// `is_within_repeat_region_interval_tree2`
    pub fn is_within_repeat_region_interval_tree2(&self, chrm1: &str, pos: i64) -> Vec<(bool, f64)> {
        match self.m_tree.get(&self.chrm(chrm1)) {
            None => Vec::new(),
            Some(t) => t.at_f(pos as f64).into_iter().map(|(b, _)| (true, b)).collect(),
        }
    }

    /// `get_div_subfamily_with_min_div`
    pub fn get_div_subfamily_with_min_div(&self, chrm1: &str, pos: f64) -> DivInfo {
        let none = DivInfo { div_rate: PyNum::Int(-1), sub_family: String::new(), family: None, pos: -1.0, pos_end: -1 };
        let chrm = self.chrm(chrm1);
        let recs = match self.m_rmsk.get(&chrm).and_then(|m| m.get(&fk(pos))) {
            None => return none,
            Some((_, r)) => r,
        };
        let mut best = &recs[0];
        for r in recs.iter().skip(1) {
            if r.div < best.div {
                best = r;
            }
        }
        DivInfo { div_rate: PyNum::Float(best.div), sub_family: best.sub_type.clone(), family: Some(best.family.clone()), pos, pos_end: best.end.trunc() as i64 }
    }

    /// `get_div_subfamily_with_min_div2`
    pub fn get_div_subfamily_with_min_div2(&self, chrm1: &str, pos: f64) -> Vec<DivInfo> {
        let chrm = self.chrm(chrm1);
        match self.m_rmsk.get(&chrm).and_then(|m| m.get(&fk(pos))) {
            None => Vec::new(),
            Some((_, recs)) => recs
                .iter()
                .map(|r| DivInfo { div_rate: PyNum::Float(r.div), sub_family: r.sub_type.clone(), family: Some(r.family.clone()), pos, pos_end: r.end.trunc() as i64 })
                .collect(),
        }
    }
}

// ------------------------------------------------------------------------------------------------
// record types

/// A transduction candidate of `call_out_candidate_transduct`:
/// (b_l_trsdct, b_r_trsdct, n_disc_trsdct, trsdct_info)
#[derive(Clone, Debug)]
pub struct TdCand {
    pub b_l: bool,
    pub b_r: bool,
    pub n_disc: i64,
    pub info: String,
}

/// Site information used to update/add records (`collect_ins_info` and the sibling variants).
#[derive(Clone, Debug)]
pub enum TdSite {
    /// (n_l_disc, lcns_start, lcns_end, n_r_disc, rcns_start, rcns_end, n_cns_clip, n_polyA_clip,
    ///  i_lcov2, i_rcov2, s_src)
    Info {
        n_l_disc: i64,
        lcns_start: PyNum,
        lcns_end: PyNum,
        n_r_disc: i64,
        rcns_start: PyNum,
        rcns_end: PyNum,
        n_cns_clip: i64,
        n_polya_clip: i64,
        lcov: PyNum,
        rcov: PyNum,
        s_src: String,
    },
    /// (blcluster, brcluster, s_src_chrm, i_src_pos, nspot)
    Sibling { bl: bool, br: bool, src_chrm: String, src_pos: i64, nspot: i64 },
}

/// A disc read aligned to a flank: (hit_ref_pos, anchor_map_pos, dist_from_rep)
pub type TdHit = (i64, i64, i64);

/// `m_trsd_picked` record: (max_source, l_pos, max_cnt)
#[derive(Clone, Debug)]
pub struct TdPicked {
    pub source: String,
    pub l_pos: Vec<TdHit>,
    pub max_cnt: i64,
}

/// `m_orphan` record: (tmp_source, n_left, n_right)
pub type TdOrphan = (String, i64, i64);

/// `_gnrt_reverse_complementary`
pub fn gnrt_reverse_complementary(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            _ => 'N',
        })
        .collect()
}

#[inline]
fn is_consecutive_polya(s: &str) -> bool {
    s.contains("AAAAA") || s.contains("AATAA")
}
#[inline]
fn is_consecutive_polyt(s: &str) -> bool {
    s.contains("TTTTT") || s.contains("TTATT")
}

/// `parse_s_cigar` -> (n_map, n_all)
fn parse_s_cigar(s_cigar: &str) -> Result<(i64, i64)> {
    let mut s_len = String::new();
    let mut n_map = 0;
    let mut n_all = 0;
    for ch in s_cigar.chars() {
        if ch.is_ascii_digit() {
            s_len.push(ch);
        } else {
            let tmp: i64 = s_len.parse().map_err(|_| anyhow!("ValueError: invalid literal for int(): '{}'", s_len))?;
            if ch == 'M' {
                n_map += tmp;
            }
            n_all += tmp;
            s_len.clear();
        }
    }
    Ok((n_map, n_all))
}

/// `check_algnmt_tag_decoy` -> (hit, n_Alu, n_L1, n_SVA)
fn check_algnmt_tag_decoy(a: &Record, f_mapped_cutoff: f64) -> Result<(bool, i64, i64, i64)> {
    let (mut n_alu, mut n_l1, mut n_sva) = (0, 0, 0);
    let mut b_hit = false;
    if let Some(s_alt) = a.get_tag_string("XA") {
        for s_rcd in s_alt.split(';') {
            let f: Vec<&str> = s_rcd.split(',').collect();
            if f.len() <= 1 {
                continue;
            }
            let s_cigar = *f.get(2).ok_or_else(|| anyhow!("IndexError: XA field"))?;
            let ratio = |c: &str| -> Result<f64> {
                let (n_map, n_all) = parse_s_cigar(c)?;
                Ok(n_map as f64 / n_all as f64)
            };
            if f[0].contains(gv::TD_DECOY_SVA) {
                if ratio(s_cigar)? >= f_mapped_cutoff {
                    b_hit = true;
                    n_sva += 1;
                }
            } else if f[0].contains(gv::TD_DECOY_LINE) {
                if ratio(s_cigar)? >= f_mapped_cutoff {
                    b_hit = true;
                    n_l1 += 1;
                }
            } else if f[0].contains(gv::TD_DECOY_ALU) && ratio(s_cigar)? >= f_mapped_cutoff {
                b_hit = true;
                n_alu += 1;
            }
        }
    }
    Ok((b_hit, n_alu, n_l1, n_sva))
}

/// The polyA/T check part of the clipped read (same windows as in x_clip_disc_filter).
fn clip_seq_ck(clipped_seq: &str, b_left: bool, b_rc: bool) -> String {
    let n = gv::CK_POLYA_SEQ_MAX;
    if !b_rc {
        if !b_left {
            py_head(clipped_seq, n).to_string()
        } else {
            py_tail(clipped_seq, n).to_string()
        }
    } else {
        let s = if b_left { py_head(clipped_seq, n) } else { py_tail(clipped_seq, n) };
        gnrt_reverse_complementary(s)
    }
}

/// `is_matched_rep_type`
pub fn is_matched_rep_type(i_rep_type: i64, s_sub_family: &str) -> bool {
    if (i_rep_type & 1 == 1) && (s_sub_family.contains("L1") || s_sub_family.contains("LINE/L1") || s_sub_family.contains("LINE1")) {
        return true;
    }
    i_rep_type & 4 == 4 && s_sub_family.contains("SVA")
}

#[inline]
fn site_get<'a, T>(m: &'a SiteMap<T>, chrm: &str, pos: i64) -> Option<&'a T> {
    m.get(chrm).and_then(|x| x.get(&pos))
}

// ------------------------------------------------------------------------------------------------
// XTransduction

pub struct XTransduction {
    pub working_folder: String,
    pub n_jobs: usize,
    pub sf_reference: String,
    pub abnormal_insert_size: i64,
    pub b_rslt_with_chr: bool,
    pub boundary_extnd: f64,
    pub same_sites_extnd: i64,
    pub max_indel_ratio: f64,
    pub f_min_uniq_ratio: f64,
}

impl XTransduction {
    pub fn new(working_folder: &str, n_jobs: usize, sf_reference: &str) -> Self {
        XTransduction {
            working_folder: working_folder.to_string(),
            n_jobs,
            sf_reference: sf_reference.to_string(),
            abnormal_insert_size: 3000,
            b_rslt_with_chr: true,
            boundary_extnd: 450.0,
            same_sites_extnd: 75,
            max_indel_ratio: 0.2,
            f_min_uniq_ratio: 0.5,
        }
    }

    pub fn set_boundary_extend(&mut self, extnd: f64) {
        self.boundary_extnd = extnd;
    }

    pub fn set_abnormal_insert_size(&mut self, i_is: i64) {
        if i_is > self.abnormal_insert_size {
            self.abnormal_insert_size = i_is;
        }
    }

    /// `prep_annotation_interval_tree`
    pub fn prep_annotation_interval_tree(&self, sf_rmsk: &str, i_min_copy_len: i64) -> Result<TdAnnotation> {
        TdAnnotation::prep(sf_rmsk, self.boundary_extnd, i_min_copy_len, self.b_rslt_with_chr)
    }

    /// `re_slct_with_clip_raw_disc_sites`
    #[allow(clippy::too_many_arguments)]
    pub fn re_slct_with_clip_raw_disc_sites(
        &self,
        sf_raw_disc: &str,
        sf_candidates: &str,
        n_disc_cutoff: i64,
        xannotation: &TdAnnotation,
        i_rep_type: i64,
        b_tumor: bool,
        sf_out: &str,
        b_with_candidate: bool,
    ) -> Result<()> {
        let m_te_candidates = xsites::load_in_candidate_list_str_version(sf_candidates)?;
        let mut fout = utils::create(sf_out)?;
        if b_with_candidate {
            for (ins_chrm, m) in m_te_candidates.iter() {
                for ins_pos in m.keys() {
                    writeln!(fout, "{}\t{}", ins_chrm, ins_pos)?;
                }
            }
        }
        for line in utils::read_raw_lines(sf_raw_disc)? {
            let fields = split_ws(&line);
            if fields.len() < 2 {
                return Err(anyhow!("IndexError: malformed line in {}: {}", sf_raw_disc, line.trim_end()));
            }
            let ins_chrm = fields[0];
            let ins_pos = utils::pint(fields[1])?;
            if let Some(m) = m_te_candidates.get(ins_chrm) {
                if m.keys().any(|&tmp_pos| xsites::are_sites_close(ins_pos, tmp_pos, gv::TWO_SITES_CLOSE_DIST)) {
                    continue;
                }
            }
            if fields.len() < 4 {
                return Err(anyhow!("IndexError: malformed line in {}: {}", sf_raw_disc, line.trim_end()));
            }
            let n_ldisc = utils::pint(fields[2])?;
            let n_rdisc = utils::pint(fields[3])?;
            let mut n_half_cutoff = crate::py::fmt::py_int(n_disc_cutoff as f64 / 2.0);
            if !b_tumor {
                n_half_cutoff = crate::py::fmt::py_int(n_disc_cutoff as f64 * 3.0 / 4.0);
            }
            if n_ldisc < n_half_cutoff && n_rdisc < n_half_cutoff {
                continue;
            }
            let (b_in_rep, i_pos) = xannotation.is_within_repeat_region_interval_tree(ins_chrm, ins_pos);
            let d = xannotation.get_div_subfamily_with_min_div(ins_chrm, i_pos);
            let b_match_rep = is_matched_rep_type(i_rep_type, &d.sub_family);
            if b_in_rep && b_match_rep && d.div_rate.f() < gv::TD_REP_DIVERGENT_CUTOFF as f64 {
                println!("[re-select step]:Filtered out: {}:{} fall in repetitive region.", ins_chrm, ins_pos);
                continue;
            }
            writeln!(fout, "{}", utils::rstrip(&line))?;
        }
        fout.flush()?;
        Ok(())
    }

    /// `pick_full_length_insertion`
    pub fn pick_full_length_insertion(&self, sf_candidates: &str, ileftmost_cns: i64) -> Result<Vec<(String, i64)>> {
        let m = xsites::load_in_candidate_list_str_version(sf_candidates)?;
        let mut l_picked = Vec::new();
        let get = |rcd: &Vec<String>, i: usize| -> Result<String> { rcd.get(i).cloned().ok_or_else(|| anyhow!("IndexError: list index out of range")) };
        for (ins_chrm, mp) in m.iter() {
            for (&ins_pos, rcd) in mp.iter() {
                if Self::reach_5_prime_end(&get(rcd, 17)?, &get(rcd, 18)?, ileftmost_cns)? {
                    l_picked.push((ins_chrm.clone(), ins_pos));
                    continue;
                }
                if Self::reach_5_prime_end(&get(rcd, 19)?, &get(rcd, 20)?, ileftmost_cns)? {
                    l_picked.push((ins_chrm.clone(), ins_pos));
                }
            }
        }
        Ok(l_picked)
    }

    fn reach_5_prime_end(s_lcluster: &str, s_rcluster: &str, i_leftmost: i64) -> Result<bool> {
        let l: Vec<&str> = s_lcluster.split(':').collect();
        let r: Vec<&str> = s_rcluster.split(':').collect();
        if l.len() < 2 || r.len() < 2 {
            return Err(anyhow!("IndexError: malformed cluster"));
        }
        let i_lstart = crate::post_filter::int_float(l[0])?;
        let _ = crate::post_filter::int_float(l[1])?;
        let i_rstart = crate::post_filter::int_float(r[0])?;
        let _ = crate::post_filter::int_float(r[1])?;
        Ok((i_lstart > 0 && i_lstart < i_leftmost) || (i_rstart > 0 && i_rstart < i_leftmost))
    }

    /// `construct_novel_flanks` (with `gnrt_flank_regions_of_polymerphic_insertions`)
    pub fn construct_novel_flanks(&self, l_fl_polymorphic: &[(String, i64)], sf_flank: &str, i_flank_lenth: i64, sf_new_flank: &str) -> Result<()> {
        let sf_polymerphic_flanks = format!("{}.polymerphic_only.fa", sf_new_flank);
        gnrt_flank_regions_of_polymerphic_insertions(l_fl_polymorphic, i_flank_lenth, &self.sf_reference, &sf_polymerphic_flanks)?;
        utils::run_cmd_to_file(&format!("cat {} {}", sf_flank, sf_polymerphic_flanks), sf_new_flank);
        Ok(())
    }

    /// `parse_transduction_clip_algnmt` -> (left clip transductions, right, polyA counters)
    #[allow(clippy::type_complexity)]
    pub fn parse_transduction_clip_algnmt(&self, sf_clip_alignmt: &str, flank_lth: i64) -> Result<(SiteMap<IndexMap<i64, Vec<String>>>, SiteMap<IndexMap<i64, Vec<String>>>, SiteMap<[i64; 6]>)> {
        let mut m_l: SiteMap<IndexMap<i64, Vec<String>>> = IndexMap::new();
        let mut m_r: SiteMap<IndexMap<i64, Vec<String>>> = IndexMap::new();
        let mut m_polya: SiteMap<[i64; 6]> = IndexMap::new();
        let mut samfile = AlignmentFile::open(sf_clip_alignmt, &self.sf_reference)?;
        let names: Vec<String> = samfile.references().to_vec();
        let mut err: Option<anyhow::Error> = None;
        samfile.fetch_all(|a: &Record| {
            if err.is_some() {
                return;
            }
            let r = (|| -> Result<()> {
                let read_info = a.query_name();
                let f: Vec<&str> = read_info.split(gv::SEPERATOR).collect();
                if f.len() < 5 {
                    return Err(anyhow!("IndexError: unexpected read name {}", read_info));
                }
                let ori_chrm = f[0].to_string();
                let ref_mpos = utils::pint(f[1])?;
                let ori_clip_flag = f[2];
                let ori_insertion_pos = utils::pint(f[4])?;
                if (ref_mpos - ori_insertion_pos).abs() > gv::NEARBY_CLIP {
                    return Ok(());
                }
                if a.is_unmap() {
                    return Ok(());
                }
                let hit_flank_id = names.get(a.tid() as usize).cloned().unwrap_or_default();
                let b_left = ori_clip_flag != gv::FLAG_RIGHT_CLIP;
                let (b_ok, _) = is_clipped_part_qualified_algnmt(&a.cigar_tuples(), gv::TD_CLIP_QLFD_RATIO);
                if !b_ok {
                    return Ok(());
                }
                let mut source_start: i64 = -1;
                let mut source_end: i64 = -1;
                let mut sourc_rc = "0".to_string();
                let mut chrm_fl_l1 = String::new();
                let flank_id_fields: Vec<&str> = hit_flank_id.split(gv::SEPERATOR).collect();
                if flank_id_fields.len() > 1 {
                    chrm_fl_l1 = flank_id_fields[0].to_string();
                    source_start = utils::pint(flank_id_fields[1])?;
                    source_end = utils::pint(*flank_id_fields.get(2).ok_or_else(|| anyhow!("IndexError"))?)?;
                    sourc_rc = flank_id_fields[flank_id_fields.len() - 1].chars().next().map(|c| c.to_string()).ok_or_else(|| anyhow!("IndexError"))?;
                    if ori_chrm == chrm_fl_l1
                        && ((ref_mpos - source_start).abs() < gv::MIN_POLYMORPHIC_SOURCE_DIST || (ref_mpos - source_end).abs() < gv::MIN_POLYMORPHIC_SOURCE_DIST)
                    {
                        return Ok(());
                    }
                }
                let clipped_seq = a.query_sequence().unwrap_or_default();
                let s_ck = clip_seq_ck(&clipped_seq, b_left, a.is_rev());
                let b_polya = is_consecutive_polya(&s_ck);
                let b_polyt = is_consecutive_polyt(&s_ck);
                if b_polya || b_polyt {
                    let p = m_polya.entry(ori_chrm.clone()).or_default().entry(ori_insertion_pos).or_insert([0; 6]);
                    let k = if b_polya && !b_polyt {
                        0
                    } else if b_polyt && !b_polya {
                        2
                    } else {
                        4
                    };
                    p[if b_left { k } else { k + 1 }] += 1;
                }
                if chrm_fl_l1.len() > 3 && chrm_fl_l1.starts_with("chr") && (ori_chrm.len() <= 3 || !ori_chrm.starts_with("chr")) {
                    chrm_fl_l1 = chrm_fl_l1[3..].to_string();
                }
                let map_pos = a.reference_start();
                let mut hit_ref_pos = ori_insertion_pos;
                if flank_id_fields.len() > 1 {
                    hit_ref_pos = utils::pint(flank_id_fields[2])? + map_pos;
                    if flank_id_fields[flank_id_fields.len() - 1].ends_with('L') {
                        hit_ref_pos = utils::pint(flank_id_fields[1])? - flank_lth + map_pos;
                    }
                }
                let s_transduct = format!("{}:{}-{}~{}~{}", chrm_fl_l1, source_start, source_end, sourc_rc, hit_ref_pos);
                let m = if b_left { &mut m_l } else { &mut m_r };
                m.entry(ori_chrm).or_default().entry(ori_insertion_pos).or_default().entry(ref_mpos).or_default().push(s_transduct);
                Ok(())
            })();
            if let Err(e) = r {
                err = Some(e);
            }
        })?;
        if let Some(e) = err {
            return Err(e);
        }
        Ok((m_l, m_r, m_polya))
    }

    /// `filter_by_clip_algnmt` -> {chrm: {pos: (n_cnt, n_polyA)}}
    pub fn filter_by_clip_algnmt(
        m_ltsdct: &SiteMap<IndexMap<i64, Vec<String>>>,
        m_rtsdct: &SiteMap<IndexMap<i64, Vec<String>>>,
        m_polya: &SiteMap<[i64; 6]>,
        n_clip_cutoff: i64,
        n_polya_cutoff: f64,
    ) -> SiteMap<(i64, i64)> {
        let mut m_clip_cnt: SiteMap<i64> = IndexMap::new();
        for m_t in [m_ltsdct, m_rtsdct] {
            for (ins_chrm, m) in m_t.iter() {
                let e = m_clip_cnt.entry(ins_chrm.clone()).or_default();
                for (&ins_pos, mp) in m.iter() {
                    let n: i64 = mp.values().map(|v| v.len() as i64).sum();
                    *e.entry(ins_pos).or_insert(0) += n;
                }
            }
        }
        let mut m_slct: SiteMap<(i64, i64)> = IndexMap::new();
        for (ins_chrm, m) in m_clip_cnt.iter() {
            for (&ins_pos, &n_cnt) in m.iter() {
                if n_cnt < n_clip_cutoff {
                    continue;
                }
                let p = match site_get(m_polya, ins_chrm, ins_pos) {
                    None => continue,
                    Some(p) => p,
                };
                let n_polya: i64 = p.iter().sum();
                if (n_polya as f64) < n_polya_cutoff {
                    continue;
                }
                m_slct.entry(ins_chrm.clone()).or_default().insert(ins_pos, (n_cnt, n_polya));
            }
        }
        m_slct
    }

    /// `re_select_disc_reads`
    pub fn re_select_disc_reads(&self, sf_disc_sam: &str, m_transdct: &SiteMap<TdPicked>, sf_out: &str) -> Result<()> {
        let mut fout = utils::create(sf_out)?;
        let mut samfile = AlignmentFile::open(sf_disc_sam, &self.sf_reference)?;
        let mut err: Option<anyhow::Error> = None;
        let mut buf = String::new();
        samfile.fetch_all(|a: &Record| {
            if err.is_some() {
                return;
            }
            let read_info = a.query_name();
            let f: Vec<&str> = read_info.split(gv::SEPERATOR).collect();
            let n = f.len();
            if n < 4 {
                err = Some(anyhow!("IndexError: unexpected read name {}", read_info));
                return;
            }
            if let Err(e) = utils::pint(f[n - 4]) {
                err = Some(e);
                return;
            }
            let ins_chrm = f[n - 3];
            let ins_pos = match utils::pint(f[n - 2]) {
                Ok(v) => v,
                Err(e) => {
                    err = Some(e);
                    return;
                }
            };
            if site_get(m_transdct, ins_chrm, ins_pos).is_some() {
                buf.push('>');
                buf.push_str(&read_info);
                buf.push('\n');
                buf.push_str(&a.query_sequence().unwrap_or_default());
                buf.push('\n');
            }
        })?;
        if let Some(e) = err {
            return Err(e);
        }
        fout.write_all(buf.as_bytes())?;
        fout.flush()?;
        Ok(())
    }

    /// `parse_transduction_disc_algnmt` -> (picked transductions, orphan transductions)
    pub fn parse_transduction_disc_algnmt(&self, sf_disc_alignmt: &str, flank_lth: i64, n_half_cutoff: f64, i_rep_type: i64) -> Result<(SiteMap<TdPicked>, SiteMap<TdOrphan>)> {
        let mut m_transduction: SiteMap<IndexMap<String, Vec<TdHit>>> = IndexMap::new();
        let mut m_raw_td_src: SiteMap<IndexMap<String, i64>> = IndexMap::new();
        let mut m_alu: SiteMap<i64> = IndexMap::new();
        let mut m_l1: SiteMap<i64> = IndexMap::new();
        let mut m_sva: SiteMap<i64> = IndexMap::new();
        let mut samfile = AlignmentFile::open(sf_disc_alignmt, &self.sf_reference)?;
        let names: Vec<String> = samfile.references().to_vec();
        let mut err: Option<anyhow::Error> = None;
        samfile.fetch_all(|a: &Record| {
            if err.is_some() {
                return;
            }
            let r = (|| -> Result<()> {
                let read_info = a.query_name();
                let f: Vec<&str> = read_info.split(gv::SEPERATOR).collect();
                let n = f.len();
                if n < 4 {
                    return Err(anyhow!("IndexError: unexpected read name {}", read_info));
                }
                let anchor_map_pos = utils::pint(f[n - 4])?;
                let ins_chrm = f[n - 3].to_string();
                let ins_pos = utils::pint(f[n - 2])?;
                if a.is_unmap() {
                    return Ok(());
                }
                let hit_flank_id = names.get(a.tid() as usize).cloned().unwrap_or_default();
                let f_mapped_cutoff = gv::F_MIN_TRSDCT_DISC_MAP_RATION;
                let (b_ok, _) = is_clipped_part_qualified_algnmt(&a.cigar_tuples(), f_mapped_cutoff);
                if !b_ok {
                    return Ok(());
                }
                let inc = |m: &mut SiteMap<i64>| {
                    *m.entry(ins_chrm.clone()).or_default().entry(ins_pos).or_insert(0) += 1;
                };
                if hit_flank_id.contains(gv::TD_DECOY_LINE)
                    || hit_flank_id.contains(gv::TD_DECOY_SVA)
                    || hit_flank_id.contains(gv::TD_DECOY_ALU)
                    || hit_flank_id.contains(gv::TD_DECOY_HERV)
                {
                    if hit_flank_id.contains(gv::TD_DECOY_ALU) {
                        inc(&mut m_alu);
                    }
                    if hit_flank_id.contains(gv::TD_DECOY_SVA) {
                        inc(&mut m_sva);
                    }
                    if hit_flank_id.contains(gv::TD_DECOY_LINE) {
                        inc(&mut m_l1);
                    }
                    return Ok(());
                }
                let (b_xa_hit, n_a, n_l, n_s) = check_algnmt_tag_decoy(a, f_mapped_cutoff)?;
                if b_xa_hit {
                    if n_a > 0 {
                        inc(&mut m_alu);
                    }
                    if n_l > 0 {
                        inc(&mut m_l1);
                    }
                    if n_s > 0 {
                        inc(&mut m_sva);
                    }
                }
                let map_pos = a.reference_start();
                let flank_id_fields: Vec<&str> = hit_flank_id.split(gv::SEPERATOR).collect();
                if flank_id_fields.len() < 3 {
                    println!("[Error]: Wrong flanking/decoy id: {}", hit_flank_id);
                    return Ok(());
                }
                let mut chrm_fl_l1 = flank_id_fields[0].to_string();
                let source_start = utils::pint(flank_id_fields[1])?;
                let source_end = utils::pint(flank_id_fields[2])?;
                let last = flank_id_fields[flank_id_fields.len() - 1];
                let sourc_rc = last.chars().next().map(|c| c.to_string()).ok_or_else(|| anyhow!("IndexError"))?;
                if chrm_fl_l1.len() > 3 && chrm_fl_l1.starts_with("chr") && (ins_chrm.len() <= 3 || !ins_chrm.starts_with("chr")) {
                    chrm_fl_l1 = chrm_fl_l1[3..].to_string();
                }
                let mut hit_ref_pos = utils::pint(flank_id_fields[2])? + map_pos;
                let mut dist_from_rep = map_pos;
                if last.ends_with('L') {
                    hit_ref_pos = utils::pint(flank_id_fields[1])? - flank_lth + map_pos;
                    dist_from_rep = flank_lth - map_pos;
                }
                if ins_chrm == chrm_fl_l1
                    && ((ins_pos - source_start).abs() < gv::MIN_POLYMORPHIC_SOURCE_DIST || (ins_pos - source_end).abs() < gv::MIN_POLYMORPHIC_SOURCE_DIST)
                {
                    println!("[DISC-TD-STEP:] Filter out {}:{}, aligned to itself flanking regions", ins_chrm, ins_pos);
                    return Ok(());
                }
                let s_source = format!("{}:{}-{}~{}", chrm_fl_l1, source_start, source_end, sourc_rc);
                if a.mapping_quality() <= gv::TRANSDCT_UNIQ_MAPQ {
                    *m_raw_td_src.entry(ins_chrm).or_default().entry(ins_pos).or_default().entry(s_source).or_insert(0) += 1;
                    return Ok(());
                }
                m_transduction.entry(ins_chrm).or_default().entry(ins_pos).or_default().entry(s_source).or_default().push((
                    hit_ref_pos,
                    anchor_map_pos,
                    dist_from_rep,
                ));
                Ok(())
            })();
            if let Err(e) = r {
                err = Some(e);
            }
        })?;
        if let Some(e) = err {
            return Err(e);
        }
        let m_decoy_slct = Self::chk_decoy_cns_cnt(&m_transduction, &m_alu, &m_l1, &m_sva, n_half_cutoff * 2.0, i_rep_type);
        let m_trsd_picked = Self::pick_trsdct_rescued(&m_decoy_slct, n_half_cutoff, &m_raw_td_src, self.f_min_uniq_ratio);
        let m_orphan = Self::slct_orphan_transduction(&m_decoy_slct, &m_trsd_picked, n_half_cutoff);
        Ok((m_trsd_picked, m_orphan))
    }

    /// `_chk_decoy_cns_cnt`
    fn chk_decoy_cns_cnt(
        m_trsdct: &SiteMap<IndexMap<String, Vec<TdHit>>>,
        m_alu: &SiteMap<i64>,
        m_l1: &SiteMap<i64>,
        m_sva: &SiteMap<i64>,
        n_cutoff: f64,
        i_rep_type: i64,
    ) -> SiteMap<IndexMap<String, Vec<TdHit>>> {
        if i_rep_type < 0 {
            return m_trsdct.clone();
        }
        let mut m_slct: SiteMap<IndexMap<String, Vec<TdHit>>> = IndexMap::new();
        for (ins_chrm, m) in m_trsdct.iter() {
            for (&ins_pos, v) in m.iter() {
                let n_alu = *site_get(m_alu, ins_chrm, ins_pos).unwrap_or(&0) as f64;
                let n_l1 = *site_get(m_l1, ins_chrm, ins_pos).unwrap_or(&0) as f64;
                let n_sva = *site_get(m_sva, ins_chrm, ins_pos).unwrap_or(&0) as f64;
                let mut b_slct = true;
                if i_rep_type & 1 != 0 {
                    if n_alu >= n_cutoff || n_sva >= n_cutoff {
                        b_slct = false;
                    }
                } else if i_rep_type & 4 != 0 && (n_alu >= n_cutoff || n_l1 >= n_cutoff) {
                    b_slct = false;
                }
                if b_slct {
                    m_slct.entry(ins_chrm.clone()).or_default().insert(ins_pos, v.clone());
                }
            }
        }
        m_slct
    }

    /// `_pick_trsdct_rescued` (the second, unused return value is dropped)
    fn pick_trsdct_rescued(m_transduction: &SiteMap<IndexMap<String, Vec<TdHit>>>, n_disc_cutoff: f64, m_raw_td_src: &SiteMap<IndexMap<String, i64>>, f_min_ratio: f64) -> SiteMap<TdPicked> {
        let mut m_picked: SiteMap<TdPicked> = IndexMap::new();
        for (ins_chrm, m) in m_transduction.iter() {
            for (&ins_pos, m_src) in m.iter() {
                let mut max_cnt = 0i64;
                let mut max_source = String::new();
                let mut n_all = 0i64;
                for (src, v) in m_src.iter() {
                    let c = v.len() as i64;
                    n_all += c;
                    if c > max_cnt {
                        max_cnt = c;
                        max_source = src.clone();
                    }
                }
                if (n_all as f64) < n_disc_cutoff {
                    println!("[DISC-TD-STEP:] Filter out {}:{}, no enough disc support!", ins_chrm, ins_pos);
                    continue;
                }
                let f_tmp_ratio = max_cnt as f64 / n_all as f64;
                let n_non_unique = m_raw_td_src.get(ins_chrm).and_then(|x| x.get(&ins_pos)).and_then(|x| x.get(&max_source)).cloned().unwrap_or(0);
                if (max_cnt as f64) >= n_disc_cutoff
                    && f_tmp_ratio >= gv::TRANSDCT_MULTI_SOURCE_MIN_RATIO
                    && (max_cnt as f64) / ((max_cnt + n_non_unique) as f64) > f_min_ratio
                {
                    let l_pos = m_src[&max_source].clone();
                    m_picked.entry(ins_chrm.clone()).or_default().insert(ins_pos, TdPicked { source: max_source, l_pos, max_cnt });
                } else {
                    println!("[DISC-TD-STEP:] Filter out {}:{}, ratio is low!", ins_chrm, ins_pos);
                }
            }
        }
        m_picked
    }

    /// `_slct_orphan_transduction`
    fn slct_orphan_transduction(m_transduction: &SiteMap<IndexMap<String, Vec<TdHit>>>, m_slcted: &SiteMap<TdPicked>, n_half_cutoff: f64) -> SiteMap<TdOrphan> {
        let mut m_orphan: SiteMap<TdOrphan> = IndexMap::new();
        for (ins_chrm, m) in m_transduction.iter() {
            for (&ins_pos, m_src) in m.iter() {
                let picked = match site_get(m_slcted, ins_chrm, ins_pos) {
                    None => continue,
                    Some(p) => p,
                };
                let (mut n_left, mut n_right) = (0i64, 0i64);
                for &(_, anchor_pos, _) in &m_src[&picked.source] {
                    if anchor_pos + 5 < ins_pos {
                        n_left += 1;
                    } else {
                        n_right += 1;
                    }
                }
                if (n_left as f64) >= n_half_cutoff && (n_right as f64) >= n_half_cutoff {
                    m_orphan.entry(ins_chrm.clone()).or_default().insert(ins_pos, (picked.source.clone(), n_left, n_right));
                }
            }
        }
        m_orphan
    }

    /// `prep_representative_pos`
    pub fn prep_representative_pos<T>(m_disc: &SiteMap<T>) -> SiteMap<i64> {
        let mut m: SiteMap<i64> = IndexMap::new();
        for (c, mp) in m_disc.iter() {
            let e = m.entry(c.clone()).or_default();
            for &p in mp.keys() {
                e.insert(p, p);
            }
        }
        m
    }

    /// `call_out_candidate_transduct`
    #[allow(clippy::too_many_arguments)]
    pub fn call_out_candidate_transduct(
        m_repsnt_pos: &SiteMap<i64>,
        m_disc_transdct: &SiteMap<TdPicked>,
        icorcord: i64,
        ratio: f64,
        m_clip_slct: &SiteMap<(i64, i64)>,
        m_disc_polya: &SiteMap<i64>,
        m_orphan: &SiteMap<TdOrphan>,
    ) -> SiteMap<TdCand> {
        let mut m_transduction: SiteMap<TdCand> = IndexMap::new();
        for (ins_chrm, m) in m_repsnt_pos.iter() {
            for (&ins_pos, &clip_pos) in m.iter() {
                if site_get(m_clip_slct, ins_chrm, ins_pos).is_none() {
                    continue;
                }
                if site_get(m_disc_polya, ins_chrm, ins_pos).is_none() {
                    continue;
                }
                let mut b_l = false;
                let mut b_r = false;
                let mut trsdct_l = -1i64;
                let mut trsdct_r = -1i64;
                let mut n_disc_trsdct = 0i64;
                if let Some(p) = site_get(m_disc_transdct, ins_chrm, ins_pos) {
                    let mut l_l: Vec<i64> = Vec::new();
                    let mut l_r: Vec<i64> = Vec::new();
                    for &(hit_ref_pos, anchor_pos, _) in &p.l_pos {
                        if anchor_pos < clip_pos {
                            l_l.push(hit_ref_pos);
                        } else {
                            l_r.push(hit_ref_pos);
                        }
                        // NOTE: inside the loop in the original
                        if l_l.len() > l_r.len() {
                            let (b, s, e) = is_disc_cluster(&mut l_l, icorcord, ratio);
                            b_l = b;
                            trsdct_l = s;
                            trsdct_r = e;
                            n_disc_trsdct = l_l.len() as i64;
                        } else {
                            let (b, s, e) = is_disc_cluster(&mut l_r, icorcord, ratio);
                            b_r = b;
                            trsdct_l = s;
                            trsdct_r = e;
                            n_disc_trsdct = l_r.len() as i64;
                        }
                    }
                    if b_l || b_r {
                        let info = format!("{}~{}-{}", p.source, trsdct_l, trsdct_r);
                        m_transduction.entry(ins_chrm.clone()).or_default().insert(ins_pos, TdCand { b_l, b_r, n_disc: n_disc_trsdct, info });
                    }
                }
            }
        }
        for (ins_chrm, m) in m_clip_slct.iter() {
            for &ins_pos in m.keys() {
                if let Some((s_source, n_left, n_right)) = site_get(m_orphan, ins_chrm, ins_pos) {
                    let info = format!("{}~orphan", s_source);
                    m_transduction.entry(ins_chrm.clone()).or_default().insert(ins_pos, TdCand { b_l: false, b_r: true, n_disc: n_left + n_right, info });
                }
            }
        }
        m_transduction
    }

    /// `fall_in_low_div_same_type_rep`
    pub fn fall_in_low_div_same_type_rep(xannotation: &TdAnnotation, ins_chrm: &str, ins_pos: i64, i_rep_type: i64) -> bool {
        for (b_in_rep, i_pos) in xannotation.is_within_repeat_region_interval_tree2(ins_chrm, ins_pos) {
            for d in xannotation.get_div_subfamily_with_min_div2(ins_chrm, i_pos) {
                if b_in_rep && is_matched_rep_type(i_rep_type, &d.sub_family) && d.div_rate.f() < gv::TD_REP_DIVERGENT_CUTOFF as f64 {
                    return true;
                }
            }
        }
        false
    }

    /// `filter_transduction` -> (kept sites, read depth of all the input sites)
    #[allow(clippy::too_many_arguments)]
    pub fn filter_transduction<T: Clone>(
        &self,
        m_ori_ins: &SiteMap<T>,
        xannotation: &TdAnnotation,
        i_rep_type: i64,
        sf_bam_list: &str,
        i_max_cov: f64,
        m_clip_polya: &SiteMap<[i64; 6]>,
        n_half_clip: f64,
    ) -> Result<(SiteMap<T>, SiteMap<[PyNum; 4]>)> {
        let rd = ReadDepth::new(&self.working_folder, self.n_jobs, &self.sf_reference);
        let l_sites: Vec<(String, i64)> = m_ori_ins.iter().flat_map(|(c, m)| m.keys().map(move |p| (c.clone(), *p))).collect();
        let m_read_depth = rd.calc_coverage_of_two_regions(&l_sites, sf_bam_list, gv::COV_SEARCH_WINDOW, gv::LOCAL_COV_WIN, gv::COV_ISD_CHK_WIN)?;
        let mut m_new: SiteMap<T> = IndexMap::new();
        for (ins_chrm, m) in m_ori_ins.iter() {
            for (&ins_pos, v) in m.iter() {
                if Self::fall_in_low_div_same_type_rep(xannotation, ins_chrm, ins_pos, i_rep_type) {
                    println!("[Transduction filter]:Filtered out: {}:{} fall in repetitive region.", ins_chrm, ins_pos);
                    continue;
                }
                let d = site_get(&m_read_depth, ins_chrm, ins_pos).ok_or_else(|| anyhow!("KeyError: {}:{}", ins_chrm, ins_pos))?;
                if d.iter().any(|x| x.f() > i_max_cov) {
                    println!("[Transduction filter]:Filtered out: {}:{} fall in alignment island.", ins_chrm, ins_pos);
                    println!(
                        "{}:{}, Details {} {} {} {} {}",
                        ins_chrm,
                        ins_pos,
                        d[0],
                        d[1],
                        d[2],
                        d[3],
                        crate::py::fmt::float_repr(i_max_cov)
                    );
                    continue;
                }
                if let Some(p) = site_get(m_clip_polya, ins_chrm, ins_pos) {
                    if (p[0] as f64) > n_half_clip && (p[1] as f64) > n_half_clip {
                        println!("[Transduction filter]:Filtered out: {}:{} has two sides polyA reads.", ins_chrm, ins_pos);
                        continue;
                    }
                    if (p[2] as f64) > n_half_clip && (p[3] as f64) > n_half_clip {
                        println!("[Transduction filter]:Filtered out: {}:{} has two sides polyT reads.", ins_chrm, ins_pos);
                        continue;
                    }
                }
                m_new.entry(ins_chrm.clone()).or_default().insert(ins_pos, v.clone());
            }
        }
        Ok((m_new, m_read_depth))
    }

    /// `collect_ins_info`
    pub fn collect_ins_info(
        m_td_ins: &SiteMap<TdCand>,
        m_disc_cns: &SiteMap<DiscFlt>,
        m_orphan: &SiteMap<TdOrphan>,
        m_clip_cns: &SiteMap<(i64, i64)>,
        m_read_depth: &SiteMap<[PyNum; 4]>,
    ) -> Result<SiteMap<TdSite>> {
        let mut m_site_info: SiteMap<TdSite> = IndexMap::new();
        for (ins_chrm, m) in m_td_ins.iter() {
            let e = m_site_info.entry(ins_chrm.clone()).or_default();
            for (&ins_pos, td) in m.iter() {
                let (mut n_l_disc, mut lcns_start, mut lcns_end, mut n_r_disc, mut rcns_start, mut rcns_end) =
                    (0i64, PyNum::Int(-1), PyNum::Int(-1), 0i64, PyNum::Int(-1), PyNum::Int(-1));
                if let Some(d) = site_get(m_disc_cns, ins_chrm, ins_pos) {
                    n_l_disc = d.n_l_disc;
                    lcns_start = d.lcns_start;
                    lcns_end = d.lcns_end;
                    n_r_disc = d.n_r_disc;
                    rcns_start = d.rcns_start;
                    rcns_end = d.rcns_end;
                } else {
                    let o = site_get(m_orphan, ins_chrm, ins_pos).ok_or_else(|| anyhow!("KeyError: {}:{}", ins_chrm, ins_pos))?;
                    n_l_disc = o.1;
                    n_r_disc = o.2;
                }
                if td.b_l && !td.b_r {
                    n_l_disc = td.n_disc;
                    lcns_start = PyNum::Int(-1);
                    lcns_end = PyNum::Int(-1);
                } else if !td.b_l && td.b_r {
                    n_r_disc = td.n_disc;
                    rcns_start = PyNum::Int(-1);
                    rcns_end = PyNum::Int(-1);
                }
                let c = site_get(m_clip_cns, ins_chrm, ins_pos).ok_or_else(|| anyhow!("KeyError: {}:{}", ins_chrm, ins_pos))?;
                let dpt = site_get(m_read_depth, ins_chrm, ins_pos).ok_or_else(|| anyhow!("KeyError: {}:{}", ins_chrm, ins_pos))?;
                e.insert(
                    ins_pos,
                    TdSite::Info {
                        n_l_disc,
                        lcns_start,
                        lcns_end,
                        n_r_disc,
                        rcns_start,
                        rcns_end,
                        n_cns_clip: c.0,
                        n_polya_clip: c.1,
                        lcov: dpt[2],
                        rcov: dpt[3],
                        s_src: td.info.clone(),
                    },
                );
            }
        }
        Ok(m_site_info)
    }

    /// `_update_existing_results`
    pub fn update_existing_results(&self, sf_xtea_rslts: &str, m_td_site_info: &SiteMap<TdSite>, i_extnd: i64, sf_updt: &str, sf_new_ins: Option<&str>) -> Result<()> {
        let mut l_rslts = XTEARsltParser::load_in_xtea_rslt_ori(sf_xtea_rslts)?;
        let mut m_overlap: SiteMap<()> = IndexMap::new();
        for rcd in l_rslts.iter_mut() {
            if rcd.len() < 2 {
                return Err(anyhow!("IndexError: list index out of range"));
            }
            let ins_chrm = rcd[0].clone();
            let ins_pos = utils::pint(&rcd[1])?;
            if let Some(m) = m_td_site_info.get(&ins_chrm) {
                for i_tmp in (ins_pos - i_extnd)..(ins_pos + i_extnd) {
                    if let Some(info) = m.get(&i_tmp) {
                        m_overlap.entry(ins_chrm.clone()).or_default().insert(i_tmp, ());
                        Self::update_one_td_rcd(rcd, info)?;
                        break;
                    }
                }
            }
        }
        let mut fout = utils::create(sf_updt)?;
        for rcd in &l_rslts {
            let mut s = String::new();
            for f in rcd {
                s.push_str(f);
                s.push('\t');
            }
            s.push('\n');
            fout.write_all(s.as_bytes())?;
        }
        fout.flush()?;
        if let Some(sf) = sf_new_ins {
            let mut fout = utils::create(sf)?;
            for (ins_chrm, m) in m_td_site_info.iter() {
                for ins_pos in m.keys() {
                    if m_overlap.get(ins_chrm).map_or(false, |x| x.contains_key(ins_pos)) {
                        continue;
                    }
                    writeln!(fout, "{}\t{}", ins_chrm, ins_pos)?;
                }
            }
            fout.flush()?;
        }
        Ok(())
    }

    /// `_update_one_td_rcd`
    fn update_one_td_rcd(rcd: &mut [String], info: &TdSite) -> Result<()> {
        match info {
            TdSite::Info { n_l_disc, n_r_disc, n_cns_clip, n_polya_clip, s_src, .. } => {
                XTEARsltParser::update_td_info(rcd, s_src, *n_cns_clip, n_l_disc + n_r_disc, *n_polya_clip)
            }
            TdSite::Sibling { src_chrm, src_pos, .. } => {
                // "{rcd_td_info[2]}:{rcd_td_info[3]}~sibling"
                XTEARsltParser::update_td_source_only(rcd, &format!("{}:{}~{}", src_chrm, src_pos, gv::SIBLING_LABEL))
            }
        }
    }

    /// `add_new_rcd_with_cutoff`
    pub fn add_new_rcd_with_cutoff(&self, m_site_info: &SiteMap<TdSite>, m_gntp: &SiteMap<GntpFeat>, ave_cov: f64, sf_out: &str, s_type: &str) -> Result<()> {
        let mut fout = utils::create(sf_out)?;
        for (ins_chrm, m) in m_gntp.iter() {
            for (&ins_pos, g) in m.iter() {
                if (g.n_disc_large_indel as f64) > ave_cov * self.max_indel_ratio {
                    continue;
                }
                let info = site_get(m_site_info, ins_chrm, ins_pos).ok_or_else(|| anyhow!("KeyError: {}:{}", ins_chrm, ins_pos))?;
                add_new_rcd(ins_chrm, ins_pos, info, g, s_type, &mut fout)?;
            }
        }
        fout.flush()?;
        Ok(())
    }

    /// `_collect_gntp_features`
    pub fn collect_gntp_features(&self, sf_sites: &str, sf_bam_list: &str, i_max_cov: f64) -> Result<SiteMap<GntpFeat>> {
        let x_gntper = XGenotyper::new(&self.sf_reference, &self.working_folder, self.n_jobs);
        let sf_gntp_feature = format!("{}.gntp.features", sf_sites);
        x_gntper.call_genotype_f(sf_bam_list, sf_sites, gv::dft_is(), &sf_gntp_feature)?;
        XGenotyper::load_in_features_from_file_with_cov_cutoff(&sf_gntp_feature, i_max_cov)
    }

    /// `update_high_confident_callset`
    pub fn update_high_confident_callset(&self, sf_ori_hc: &str, sf_all_cns: &str, sf_new_hc: &str) -> Result<()> {
        let mut m_cns: SiteMap<String> = IndexMap::new();
        let mut m_all_td: SiteMap<String> = IndexMap::new();
        for line in utils::read_raw_lines(sf_all_cns)? {
            let fields: Vec<String> = split_ws(&line).into_iter().map(|s| s.to_string()).collect();
            if fields.len() < 2 {
                continue;
            }
            let ins_chrm = fields[0].clone();
            let ins_pos = utils::pint(&fields[1])?;
            let s = utils::rstrip(&line).to_string();
            m_cns.entry(ins_chrm.clone()).or_default().insert(ins_pos, s.clone());
            if XTEARsltParser::is_transduction(&fields)? {
                m_all_td.entry(ins_chrm).or_default().insert(ins_pos, s);
            }
        }
        let mut fout = utils::create(sf_new_hc)?;
        let mut m_saved: SiteMap<()> = IndexMap::new();
        for line in utils::read_raw_lines(sf_ori_hc)? {
            let fields = split_ws(&line);
            if fields.len() < 2 {
                continue;
            }
            let ins_chrm = fields[0];
            let ins_pos = utils::pint(fields[1])?;
            let s = site_get(&m_cns, ins_chrm, ins_pos).ok_or_else(|| anyhow!("KeyError: {}:{}", ins_chrm, ins_pos))?;
            writeln!(fout, "{}", s)?;
            m_saved.entry(ins_chrm.to_string()).or_default().insert(ins_pos, ());
        }
        for (ins_chrm, m) in m_all_td.iter() {
            for (ins_pos, s) in m.iter() {
                if m_saved.get(ins_chrm).map_or(false, |x| x.contains_key(ins_pos)) {
                    continue;
                }
                writeln!(fout, "{}", s)?;
            }
        }
        fout.flush()?;
        Ok(())
    }

    /// `call_candidate_transduction_v3`
    #[allow(clippy::too_many_arguments)]
    pub fn call_candidate_transduction_v3(
        &self,
        sf_raw_slct: &str,
        sf_te_candidates: &str,
        x_filter: &XClipDiscFilter,
        sf_flank: &str,
        sf_cns: &str,
        i_flank_length: i64,
        extnd: i64,
        bin_size: i64,
        n_clip_cutoff: i64,
        n_disc_cutoff: i64,
        i_concord_dist: i64,
        f_concord_ratio: f64,
        xannotation: &TdAnnotation,
        sf_bam_list: &str,
        i_rep_type: i64,
        i_max_cov: f64,
        ave_cov: f64,
        sf_out: &str,
    ) -> Result<()> {
        let ileftmost_cns = if gv::is_call_sva() { 5000 } else { 500 };
        let l_fl_polymorpic = self.pick_full_length_insertion(sf_te_candidates, ileftmost_cns)?;
        let sf_flank_with_poly = format!("{}all_with_polymerphic_flanks.fa", self.working_folder);
        self.construct_novel_flanks(&l_fl_polymorpic, sf_flank, i_flank_length, &sf_flank_with_poly)?;
        let bwa_align = BwAlign::new(gv::BWA_PATH, gv::BWA_REALIGN_CUTOFF, self.n_jobs as i64);
        utils::run_cmd_small_output(&format!("{} index {}", gv::BWA_PATH, sf_flank_with_poly));

        let sf_clip_fq = format!("{}raw_candidate_sites_all_clip.fq", self.working_folder);
        x_filter.collect_clipped_reads(sf_raw_slct, extnd, &sf_clip_fq)?;
        let sf_clip_algnmt = format!("{}temp_transduction_clip.sam", self.working_folder);
        bwa_align.realign_clipped_read_with_polya(&sf_flank_with_poly, &sf_clip_fq, &sf_clip_algnmt)?;
        XClipDiscFilter::clean_file_by_path(&sf_clip_fq);

        let (m_lclip_trsd, m_rclip_trsd, m_clip_polya) = self.parse_transduction_clip_algnmt(&sf_clip_algnmt, i_flank_length)?;
        let mut n_polya_cutoff = n_clip_cutoff as f64 / 2.0 - 1.0;
        if n_polya_cutoff < 0.0 {
            n_polya_cutoff = 0.0;
        }
        let m_clip_slct = Self::filter_by_clip_algnmt(&m_lclip_trsd, &m_rclip_trsd, &m_clip_polya, n_clip_cutoff, n_polya_cutoff);
        let sf_clip_tmp_out = format!("{}.unique_trsdct_half_clip", sf_out);
        {
            let mut f = utils::create(&sf_clip_tmp_out)?;
            for (c, m) in m_clip_slct.iter() {
                for p in m.keys() {
                    writeln!(f, "{}\t{}", c, p)?;
                }
            }
        }
        let sf_disc_fa = format!("{}raw_candidate_sites_all_disc.fa", self.working_folder);
        x_filter.collect_disc_reads(&sf_clip_tmp_out, extnd, bin_size, &sf_disc_fa)?;
        let sf_disc_algnmt = format!("{}temp_transduction_disc.sam", self.working_folder);
        bwa_align.realign_disc_reads(&sf_flank_with_poly, &sf_disc_fa, &sf_disc_algnmt);
        XClipDiscFilter::clean_file_by_path(&sf_disc_fa);

        let mut n_half_cutoff = n_disc_cutoff as f64 / 2.0;
        if n_half_cutoff <= 1.0 {
            n_half_cutoff = n_disc_cutoff as f64;
        }
        let (m_disc_transd, m_orphan) = self.parse_transduction_disc_algnmt(&sf_disc_algnmt, i_flank_length, n_half_cutoff, i_rep_type)?;
        {
            let mut f = utils::create(&format!("{}.unique_trsdct_disc_only", sf_out))?;
            for (c, m) in m_disc_transd.iter() {
                for (p, r) in m.iter() {
                    writeln!(f, "{}\t{}\t{}", c, p, r.source)?;
                }
            }
        }
        {
            let mut f = utils::create(&format!("{}.unique_trsdct_disc_orphan", sf_out))?;
            for (c, m) in m_orphan.iter() {
                for (p, r) in m.iter() {
                    writeln!(f, "{}\t{}\t{}\t{}\t{}", c, p, r.0, r.1, r.2)?;
                }
            }
        }
        let sf_disc_fa2 = format!("{}raw_candidate_sites_all_disc_focal_sites.fa", self.working_folder);
        self.re_select_disc_reads(&sf_disc_algnmt, &m_disc_transd, &sf_disc_fa2)?;
        let sf_disc_algnmt_cns = format!("{}temp_transduction_disc_cns.sam", self.working_folder);
        bwa_align.realign_disc_reads(sf_cns, &sf_disc_fa2, &sf_disc_algnmt_cns);
        let disc = x_filter.parse_disc_algnmt_consensus(&sf_disc_algnmt_cns, gv::MIN_DISC_MAPPED_RATIO)?;
        let m_represent_pos = Self::prep_representative_pos(&m_disc_transd);
        let m_disc_filtered = XClipDiscFilter::check_seprt_disc_consistency(&m_represent_pos, &disc.m_disc_pos, i_concord_dist, f_concord_ratio, n_half_cutoff, false);
        {
            let mut f = utils::create(&format!("{}.unique_trsdct_disc_only_half_clip_half_disc", sf_out))?;
            for (c, m) in m_disc_filtered.iter() {
                for p in m.keys() {
                    writeln!(f, "{}\t{}", c, p)?;
                }
            }
        }
        let m_represent_pos2 = Self::prep_representative_pos(&m_disc_filtered);
        let m_transduct_candidates =
            Self::call_out_candidate_transduct(&m_represent_pos2, &m_disc_transd, extnd, f_concord_ratio, &m_clip_slct, &disc.m_polya, &m_orphan);
        {
            let mut f = utils::create(&format!("{}.unique_trsdct_disc_only_half_clip_half_disc_polyA", sf_out))?;
            for (c, m) in m_transduct_candidates.iter() {
                for p in m.keys() {
                    writeln!(f, "{}\t{}", c, p)?;
                }
            }
        }
        let (m_new_ins, m_depth) =
            self.filter_transduction(&m_transduct_candidates, xannotation, i_rep_type, sf_bam_list, i_max_cov, &m_clip_polya, n_clip_cutoff as f64 / 2.0)?;
        let m_site_info = Self::collect_ins_info(&m_new_ins, &m_disc_filtered, &m_orphan, &m_clip_slct, &m_depth)?;
        {
            let mut f = utils::create(&format!("{}.unique_trsdct_disc_only_half_clip_half_disc_polyA_after_filter", sf_out))?;
            for (c, m) in m_site_info.iter() {
                for p in m.keys() {
                    writeln!(f, "{}\t{}", c, p)?;
                }
            }
        }
        let sf_new_ins = format!("{}.tmp_new_sites_position_only", sf_out);
        self.update_existing_results(sf_te_candidates, &m_site_info, self.same_sites_extnd, sf_out, Some(&sf_new_ins))?;
        let m_gntp = self.collect_gntp_features(&sf_new_ins, sf_bam_list, i_max_cov)?;
        let sf_new_rcds = format!("{}{}", sf_out, gv::TD_NEW_SITES_SUFFIX);
        self.add_new_rcd_with_cutoff(&m_site_info, &m_gntp, ave_cov, &sf_new_rcds, gv::ONE_SIDE_TRSDCT)?;
        let new_lines = utils::read_raw_lines(&sf_new_rcds)?;
        let mut fout = utils::append(sf_out)?;
        for l in new_lines {
            fout.write_all(l.as_bytes())?;
        }
        fout.flush()?;
        Ok(())
    }
}

/// `XTEARsltParser.add_new_rcd`
pub fn add_new_rcd(ins_chrm: &str, ins_pos: i64, info: &TdSite, g: &GntpFeat, s_type: &str, fout: &mut dyn Write) -> Result<()> {
    let (n_l_disc, lcns_start, lcns_end, n_r_disc, rcns_start, rcns_end, n_cns_clip, n_polya_clip, lcov, rcov, s_src) = match info {
        TdSite::Info { n_l_disc, lcns_start, lcns_end, n_r_disc, rcns_start, rcns_end, n_cns_clip, n_polya_clip, lcov, rcov, s_src } => {
            (*n_l_disc, *lcns_start, *lcns_end, *n_r_disc, *rcns_start, *rcns_end, *n_cns_clip, *n_polya_clip, *lcov, *rcov, s_src.clone())
        }
        TdSite::Sibling { .. } => return Err(anyhow!("IndexError: tuple index out of range")),
    };
    let s = format!(
        "{c}\t{p}\t{p}\t-1\t-1\t{n6}\t0\t{n0}\t{n3}\t{n7}\t0\t{lcov}\t{rcov}\t0\t0\t0\t0\t0\t0\t-1:-1\t-1:-1\t{ls}:{le}\t{rs}:{re}\t{src}\t{n6}\t{nd}\t{n7}\t0\t0\t0\t0\t{inv}\t{stype}\t{oec}\t{hit}\t{g0}\t{g1}\t{g2}\t{g3}\t{g4}\t{g5}\t{g6}\t{g7}\t0\t0\t0\t0\t0\n",
        c = ins_chrm,
        p = ins_pos,
        n6 = n_cns_clip,
        n0 = n_l_disc,
        n3 = n_r_disc,
        n7 = n_polya_clip,
        lcov = lcov,
        rcov = rcov,
        ls = lcns_start,
        le = lcns_end,
        rs = rcns_start,
        re = rcns_end,
        src = s_src,
        nd = n_l_disc + n_r_disc,
        inv = gv::NOT_FIVE_PRIME_INV,
        stype = s_type,
        oec = gv::ONE_END_CONSISTNT,
        hit = gv::HIT_END_OF_CNS,
        g0 = g.n_af_clip,
        g1 = g.n_full_map,
        g2 = g.n_l_raw_clip,
        g3 = g.n_r_raw_clip,
        g4 = g.n_disc_pairs,
        g5 = g.n_concd_pairs,
        g6 = g.n_disc_large_indel,
        g7 = g.s_clip_lens,
    );
    fout.write_all(s.as_bytes())?;
    Ok(())
}

/// `XReference.gnrt_flank_regions_of_polymerphic_insertions`
pub fn gnrt_flank_regions_of_polymerphic_insertions(l_sites: &[(String, i64)], i_extend1: i64, sf_ref: &str, sf_out: &str) -> Result<()> {
    let i_extend = i_extend1.abs();
    let f_fa = FastaFile::open(sf_ref)?;
    let b_with_chr = f_fa.has_reference("chr1");
    let mut fout = utils::create(sf_out)?;
    let d = gv::S_DELIM;
    for (ins_chrm, ins_pos) in l_sites {
        let ins_pos = *ins_pos;
        if ins_pos - i_extend < 0 {
            continue;
        }
        let istart = ins_pos - i_extend;
        let iend = ins_pos + i_extend;
        let ref_chrm = crate::annotation::process_chrm_name(ins_chrm, b_with_chr);
        if !f_fa.has_reference(&ref_chrm) {
            continue;
        }
        let s_left = f_fa.fetch(&ref_chrm, istart, ins_pos)?;
        writeln!(fout, ">{}{}{}{}{}{}{}{}0L", ins_chrm, d, ins_pos, d, ins_pos, d, gv::S_POLYMORPHIC, d)?;
        writeln!(fout, "{}", s_left)?;
        let s_right = f_fa.fetch(&ref_chrm, ins_pos, iend)?;
        writeln!(fout, ">{}{}{}{}{}{}{}{}0R", ins_chrm, d, ins_pos, d, ins_pos, d, gv::S_POLYMORPHIC, d)?;
        writeln!(fout, "{}", s_right)?;
    }
    fout.flush()?;
    Ok(())
}
