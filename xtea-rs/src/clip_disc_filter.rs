//! Port of `x_clip_disc_filter.py`: `XClipDisc` (collection of clipped parts and discordant
//! reads around candidate sites) and `XClipDiscFilter` (the "-N" consensus based filter).
//!
//! Notes on fidelity:
//! * Python 3 true division is kept where it produces floats (e.g. the disc read position on the
//!   consensus `reference_start + len(seq) / 2`, the "half" cutoffs `n / 2`), and such values are
//!   written with Python's `str()` (see [`PyNum`]).
//! * dict iteration order is reproduced with `IndexMap`.

use crate::alignments::BamInfo;
use crate::annotation::process_chrm_name;
use crate::bam::{AlignmentFile, PyAln, Record, BAM_CDEL, BAM_CHARD_CLIP, BAM_CINS, BAM_CMATCH, BAM_CSOFT_CLIP};
use crate::bwa_align::BwAlign;
use crate::chrom::is_decoy_contig_chrms;
use crate::cluster::{find_first_second_peak, is_disc_cluster_f};
use crate::coverage::ReadDepth;
use crate::genotype_feature::{GntpFeat, XGenotyper};
use crate::globals as gv;
use crate::intermediate_sites::SiteMap;
use crate::num::PyNum;
use crate::polya;
use crate::utils::{self, split_ws};
use anyhow::{anyhow, Result};
use indexmap::IndexMap;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as FmtWrite;
use std::io::Write;

// ------------------------------------------------------------------------------------------------
// small helpers

/// Python `s[:n]` for n >= 0
#[inline]
pub(crate) fn py_head(s: &str, n: i64) -> &str {
    let n = n.max(0) as usize;
    if n >= s.len() {
        s
    } else {
        &s[..n]
    }
}

/// Python `s[-n:]` (note `s[-0:]` is the whole string)
#[inline]
pub(crate) fn py_tail(s: &str, n: i64) -> &str {
    let n = n.max(0) as usize;
    if n == 0 || n >= s.len() {
        s
    } else {
        &s[s.len() - n..]
    }
}

#[inline]
fn py_head_q(q: &[u8], n: i64) -> &[u8] {
    let n = n.max(0) as usize;
    if n >= q.len() {
        q
    } else {
        &q[..n]
    }
}

#[inline]
fn py_tail_q(q: &[u8], n: i64) -> &[u8] {
    let n = n.max(0) as usize;
    if n == 0 || n >= q.len() {
        q
    } else {
        &q[q.len() - n..]
    }
}

/// `is_discordant` (XClipDisc / XGenotyper)
#[inline]
pub fn is_discordant(chrm: &str, map_pos: i64, mate_chrm: &str, mate_pos: i64, is_threshold: i64) -> bool {
    if chrm != mate_chrm {
        return true;
    }
    (mate_pos - map_pos).abs() > is_threshold
}

/// `is_two_side_clipped`
#[inline]
fn is_two_side_clipped(l_cigar: &[(u32, u32)], i_min_clip: i64) -> bool {
    let first = l_cigar[0];
    let last = l_cigar[l_cigar.len() - 1];
    last.0 == BAM_CSOFT_CLIP && first.0 == BAM_CSOFT_CLIP && (last.1 as i64) > i_min_clip && (first.1 as i64) > i_min_clip
}

/// Position after the reference-consuming part of the cigar (all ops but S, H, I), applied to
/// all but the last operation -- the "exact right clip position" computation used throughout.
#[inline]
pub(crate) fn right_clip_pos(l_cigar: &[(u32, u32)], map_pos: i64) -> i64 {
    let mut p = map_pos;
    for &(t, l) in &l_cigar[..l_cigar.len() - 1] {
        if t == BAM_CSOFT_CLIP || t == BAM_CHARD_CLIP || t == BAM_CINS {
            continue;
        }
        p += l as i64;
    }
    p
}

/// `XLog`
pub struct XLog {
    f: Option<std::io::BufWriter<std::fs::File>>,
}

impl XLog {
    /// `open_file`: truncate, write the header, and keep appending.
    pub fn open_file(sf_log: &str) -> Result<XLog> {
        utils::write_string(sf_log, "[Log starts:]\n")?;
        let f = utils::append(sf_log)?;
        Ok(XLog { f: Some(f) })
    }
    pub fn append(&mut self, s: &str) {
        if let Some(f) = self.f.as_mut() {
            let _ = f.write_all(s.as_bytes());
            let _ = f.write_all(b"\n");
            let _ = f.flush();
        }
    }
    pub fn close(&mut self) {
        if let Some(mut f) = self.f.take() {
            let _ = f.flush();
        }
    }
}

// ------------------------------------------------------------------------------------------------
// XClipDisc

/// Output of `collect_clipped_disc_reads_by_region` for one site.
struct SiteCollect {
    chrm: String,
    pos: i64,
    fq: String,
    disc_pos: String,
    n_lowq: i64,
}

pub struct XClipDisc {
    pub sf_bam: String,
    pub working_folder: String,
    pub n_jobs: usize,
    pub sf_reference: String,
}

impl XClipDisc {
    pub fn new(sf_bam: &str, working_folder: &str, n_jobs: usize, sf_ref: &str) -> Self {
        XClipDisc { sf_bam: sf_bam.to_string(), working_folder: working_folder.to_string(), n_jobs, sf_reference: sf_ref.to_string() }
    }

    fn pool(&self) -> rayon::ThreadPool {
        rayon::ThreadPoolBuilder::new().num_threads(self.n_jobs.max(1)).build().unwrap()
    }

    /// `collect_clipped_disc_reads_by_region` (returns the content of the per-site files).
    fn collect_clipped_disc_reads_by_region(&self, chrm: &str, insertion_pos: i64, extnd: i64) -> Result<SiteCollect> {
        let mut start_pos = insertion_pos - extnd;
        if start_pos <= 0 {
            start_pos = 1;
        }
        let end_pos = insertion_pos + extnd;
        let bam_info = BamInfo::new(&self.sf_bam, &self.sf_reference);
        let b_with_chr = bam_info.is_chrm_contain_chr()?;
        let chrm_in_bam = process_chrm_name(chrm, b_with_chr);

        let mut samfile = AlignmentFile::open(&self.sf_bam, &self.sf_reference)?;
        let names: Vec<String> = samfile.references().to_vec();
        let n_refs = names.len() as i32;
        let mut fq = String::new();
        let mut disc = String::new();
        let mut n_cnt_clip: i64 = 0;
        let mut n_cnt_low_mapq_clip: i64 = 0;
        let min_clip_mapq = gv::minimum_clip_mapq();
        let min_disc_mapq = gv::minimum_disc_mapq();
        let read_length = gv::read_length();
        let sep = gv::SEPERATOR;
        samfile.fetch(&chrm_in_bam, start_pos, end_pos, |a: &Record| {
            if a.is_dup() {
                return;
            }
            let b_first = !a.is_read2();
            if a.is_unmap() {
                return;
            }
            let l_cigar = a.cigar_tuples();
            if l_cigar.is_empty() {
                return;
            }
            let query_name = a.query_name();
            let query_seq = a.query_sequence().unwrap_or_default();
            let query_quality = a.query_qualities().unwrap_or_default();
            let mut map_pos = a.reference_start();
            let mut mate_chrm = "*".to_string();
            let mut mate_pos: i64 = 0;
            let is_rc = if a.is_rev() { 1 } else { 0 };
            let ntid = a.next_reference_id();
            if ntid >= 0 && ntid < n_refs && !a.mate_is_unmapped() {
                mate_chrm = names[ntid as usize].clone();
                mate_pos = a.next_reference_start();
            }
            let first = l_cigar[0];
            let last = l_cigar[l_cigar.len() - 1];
            if a.mapping_quality() < min_clip_mapq {
                if a.mapping_quality() < gv::CLIP_LOW_MAPQ {
                    if first.0 == BAM_CSOFT_CLIP {
                        if (map_pos - insertion_pos).abs() < gv::NEARBY_LOW_Q_CLIP {
                            n_cnt_low_mapq_clip += 1;
                        }
                    } else if last.0 == BAM_CSOFT_CLIP {
                        let tmp_map_lenth = read_length - last.1 as f64;
                        if ((map_pos as f64) + tmp_map_lenth - insertion_pos as f64).abs() < gv::NEARBY_LOW_Q_CLIP as f64 {
                            n_cnt_low_mapq_clip += 1;
                        }
                    }
                }
                return;
            }
            if is_two_side_clipped(&l_cigar, gv::TWO_SIDE_CLIP_MIN_LEN) {
                return;
            }
            let mut b_left_clip = 0;
            if first.0 == BAM_CSOFT_CLIP {
                if a.is_supp() || a.is_sec() {
                    return;
                }
                let clipped_seq = py_head(&query_seq, first.1 as i64);
                let len_clip_seq = clipped_seq.len() as i64;
                if len_clip_seq >= gv::MINIMUM_POLYA_CLIP && len_clip_seq <= gv::BWA_REALIGN_CUTOFF {
                    if !polya::contain_enough_a_t(clipped_seq, len_clip_seq - 1) {
                        return;
                    }
                } else if len_clip_seq < gv::BWA_REALIGN_CUTOFF {
                    return;
                }
                let l_quality_score = py_head_q(&query_quality, first.1 as i64);
                if !crate::clip_read::is_qualified_clip(l_quality_score) {
                    return;
                }
                let clipped_qulity = crate::clip_read::cvt_to_ascii_quality(l_quality_score);
                let rname = format!(
                    "{}{}{}{}{}{}{}{}{}{}{}",
                    chrm, sep, map_pos, sep, gv::FLAG_LEFT_CLIP, sep, is_rc, sep, insertion_pos, sep, n_cnt_clip
                );
                n_cnt_clip += 1;
                let _ = write!(fq, "@{}\n{}\n+\n{}\n", rname, clipped_seq, clipped_qulity);
                b_left_clip = 1;
            }
            let mut b_right_clip = 0;
            if last.0 == BAM_CSOFT_CLIP {
                map_pos = right_clip_pos(&l_cigar, map_pos);
                if a.is_supp() || a.is_sec() {
                    return;
                }
                let rname = format!(
                    "{}{}{}{}{}{}{}{}{}{}{}",
                    chrm, sep, map_pos, sep, gv::FLAG_RIGHT_CLIP, sep, is_rc, sep, insertion_pos, sep, n_cnt_clip
                );
                n_cnt_clip += 1;
                let clipped_seq = py_tail(&query_seq, last.1 as i64);
                let len_clip_seq = clipped_seq.len() as i64;
                if len_clip_seq >= gv::MINIMUM_POLYA_CLIP && len_clip_seq <= gv::BWA_REALIGN_CUTOFF {
                    if !polya::contain_enough_a_t(clipped_seq, len_clip_seq - 1) {
                        return;
                    }
                }
                if len_clip_seq < gv::BWA_REALIGN_CUTOFF {
                    return;
                }
                let l_quality_score = py_tail_q(&query_quality, last.1 as i64);
                if !crate::clip_read::is_qualified_clip(l_quality_score) {
                    return;
                }
                let clipped_qulity = crate::clip_read::cvt_to_ascii_quality(l_quality_score);
                let _ = write!(fq, "@{}\n{}\n+\n{}\n", rname, clipped_seq, clipped_qulity);
                b_right_clip = 1;
            }
            if mate_chrm == "*" {
                return;
            }
            if is_decoy_contig_chrms(&mate_chrm) {
                return;
            }
            if a.mapping_quality() < min_disc_mapq {
                return;
            }
            let is_mate_rc = if a.mate_is_rev() { 1 } else { 0 };
            if is_discordant(&chrm_in_bam, map_pos, &mate_chrm, mate_pos, gv::DISC_THRESHOLD) {
                let s_mate_first = if b_first { 0 } else { 1 };
                let _ = writeln!(
                    disc,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    chrm, map_pos, is_rc, is_mate_rc, mate_chrm, mate_pos, query_name, s_mate_first, insertion_pos, b_left_clip, b_right_clip
                );
            }
        })?;
        Ok(SiteCollect { chrm: chrm.to_string(), pos: insertion_pos, fq, disc_pos: disc, n_lowq: n_cnt_low_mapq_clip })
    }

    fn load_sites(sf_candidate_list: &str) -> Result<Vec<(String, i64)>> {
        let mut v = Vec::new();
        for line in utils::read_lines(sf_candidate_list)? {
            let fields = split_ws(&line);
            if fields.len() < 2 {
                if fields.is_empty() {
                    continue;
                }
                return Err(anyhow!("malformed candidate line: {}", line));
            }
            v.push((fields[0].to_string(), utils::pint(fields[1])?));
        }
        Ok(v)
    }

    fn collect_sites(&self, sites: &[(String, i64)], extnd: i64) -> Result<Vec<SiteCollect>> {
        let res: Vec<Result<SiteCollect>> = self.pool().install(|| {
            utils::par_map(sites.len(), |i| self.collect_clipped_disc_reads_by_region(&sites[i].0, sites[i].1, extnd))
        });
        res.into_iter().collect()
    }

    /// `collect_clipped_disc_reads_of_given_list` -> list of (chrm, pos, n_low_mapq_clip)
    pub fn collect_clipped_disc_reads_of_given_list(
        &self,
        sf_candidate_list: &str,
        extnd: i64,
        bin_size: i64,
        sf_all_clip_fq: &str,
        sf_disc_fa: &str,
    ) -> Result<Vec<(String, i64, i64)>> {
        let sites = Self::load_sites(sf_candidate_list)?;
        let collected = self.collect_sites(&sites, extnd)?;
        let l_low_mapq: Vec<(String, i64, i64)> = collected.iter().map(|c| (c.chrm.clone(), c.pos, c.n_lowq)).collect();
        let sf_all_disc_pos = format!("{}all_disc_pos{}", self.working_folder, gv::DISC_POS_SUFFIX);
        {
            let mut fout_clip = utils::create(sf_all_clip_fq)?;
            let mut fout_disc = utils::create(&sf_all_disc_pos)?;
            // the per-site temporary files are removed after merging, so a duplicated site is only
            // merged once
            let mut consumed: HashSet<(String, i64)> = HashSet::new();
            for c in collected.iter() {
                if !consumed.insert((c.chrm.clone(), c.pos)) {
                    continue;
                }
                fout_clip.write_all(c.fq.as_bytes())?;
                fout_disc.write_all(c.disc_pos.as_bytes())?;
            }
            fout_clip.flush()?;
            fout_disc.flush()?;
        }
        let bam_info = BamInfo::new(&self.sf_bam, &self.sf_reference);
        bam_info.extract_mate_reads_by_name(&sf_all_disc_pos, bin_size, &self.working_folder, self.n_jobs, sf_disc_fa)?;
        Ok(l_low_mapq)
    }

    /// `collect_clipped_reads_of_given_list` (clipped parts only; sites grouped by chromosome)
    pub fn collect_clipped_reads_of_given_list(&self, sf_candidate_list: &str, extnd: i64, sf_all_clip_fq: &str) -> Result<()> {
        let sites = Self::load_sites(sf_candidate_list)?;
        let mut m_chrm_pos: IndexMap<String, Vec<i64>> = IndexMap::new();
        for (c, p) in sites {
            m_chrm_pos.entry(c).or_default().push(p);
        }
        let sf_all_disc_pos = format!("{}all_disc_pos{}", self.working_folder, gv::DISC_POS_SUFFIX);
        let mut fout_clip = utils::create(sf_all_clip_fq)?;
        let mut fout_disc = utils::create(&sf_all_disc_pos)?;
        for (chrm, l_pos) in m_chrm_pos.iter() {
            let l_sites: Vec<(String, i64)> = l_pos.iter().map(|p| (chrm.clone(), *p)).collect();
            let collected = self.collect_sites(&l_sites, extnd)?;
            let mut consumed: HashSet<i64> = HashSet::new();
            for c in collected.iter() {
                if !consumed.insert(c.pos) {
                    continue;
                }
                fout_clip.write_all(c.fq.as_bytes())?;
                fout_disc.write_all(c.disc_pos.as_bytes())?;
            }
        }
        fout_clip.flush()?;
        fout_disc.flush()?;
        Ok(())
    }

    /// `collect_disc_reads_of_given_list`: re-uses the "all_disc_pos" file of the working folder.
    pub fn collect_disc_reads_of_given_list(&self, _sf_candidate_list: &str, _extnd: i64, bin_size: i64, sf_disc_fa: &str) -> Result<()> {
        let sf_all_disc_pos = format!("{}all_disc_pos{}", self.working_folder, gv::DISC_POS_SUFFIX);
        let bam_info = BamInfo::new(&self.sf_bam, &self.sf_reference);
        bam_info.extract_mate_reads_by_name(&sf_all_disc_pos, bin_size, &self.working_folder, self.n_jobs, sf_disc_fa)
    }

    /// `calc_AF_of_site` for one chromosome (returns the lines of the per-chromosome file).
    fn calc_af_of_chrm(&self, chrm: &str, l_bams: &[String], m_candidate_sites: &IndexMap<String, Vec<i64>>, extend: i64, clip_extnd: i64, af_cutoff: f64) -> Result<String> {
        // m_rslts: ins_pos -> [(cnt_clip, clip_range, cnt_all, irange, cnt_disc, cnt_all_2)]
        let mut m_rslts: IndexMap<i64, Vec<(i64, i64, i64, i64, i64, i64)>> = IndexMap::new();
        let empty: Vec<i64> = Vec::new();
        for sf_bam in l_bams {
            let bam_info = BamInfo::new(sf_bam, &self.sf_reference);
            let b_with_chr = bam_info.is_chrm_contain_chr()?;
            let chrm_in_bam = process_chrm_name(chrm, b_with_chr);
            let mut samfile = AlignmentFile::open(sf_bam, &self.sf_reference)?;
            let names: Vec<String> = samfile.references().to_vec();
            let n_refs = names.len() as i32;
            let l_sites = m_candidate_sites.get(chrm).ok_or_else(|| anyhow!("KeyError: {}", chrm)).unwrap_or(&empty);
            for &insertion_pos in l_sites {
                let mut start_pos = insertion_pos - extend;
                if start_pos <= 0 {
                    start_pos = 1;
                }
                let end_pos = insertion_pos + extend;
                let (mut cnt_all_reads, mut cnt_clip, mut cnt_disc, mut cnt_all_2) = (0i64, 0i64, 0i64, 0i64);
                let (mut lmost_region, mut rmost_region) = (-1i64, -1i64);
                samfile.fetch(&chrm_in_bam, start_pos, end_pos, |a: &Record| {
                    if a.is_dup() || a.is_unmap() {
                        return;
                    }
                    let l_cigar = a.cigar_tuples();
                    if l_cigar.is_empty() {
                        return;
                    }
                    if a.mapping_quality() < gv::minimum_disc_mapq() {
                        return;
                    }
                    if a.is_supp() || a.is_sec() {
                        return;
                    }
                    let map_pos = a.reference_start();
                    if lmost_region == -1 {
                        lmost_region = map_pos;
                    }
                    if rmost_region == -1 {
                        rmost_region = map_pos;
                    }
                    if map_pos < lmost_region {
                        lmost_region = map_pos;
                    }
                    if map_pos > rmost_region {
                        rmost_region = map_pos;
                    }
                    cnt_all_reads += 1;
                    if (insertion_pos - map_pos).abs() < extend {
                        cnt_all_2 += 1;
                    }
                    let ntid = a.next_reference_id();
                    if ntid >= 0 && ntid < n_refs && !a.mate_is_unmapped() {
                        let mate_chrm = &names[ntid as usize];
                        let mate_pos = a.next_reference_start();
                        if is_discordant(&chrm_in_bam, map_pos, mate_chrm, mate_pos, gv::DISC_THRESHOLD) {
                            cnt_disc += 1;
                        }
                    }
                    let first = l_cigar[0];
                    let last = l_cigar[l_cigar.len() - 1];
                    if first.0 == BAM_CSOFT_CLIP {
                        if (insertion_pos - map_pos).abs() <= clip_extnd {
                            cnt_clip += 1;
                        }
                    } else if last.0 == BAM_CSOFT_CLIP {
                        let clip_pos = right_clip_pos(&l_cigar, map_pos);
                        if (insertion_pos - clip_pos).abs() <= clip_extnd {
                            cnt_clip += 1;
                        }
                    }
                })?;
                let irange = (rmost_region - lmost_region).abs();
                m_rslts.entry(insertion_pos).or_default().push((cnt_clip, 2 * clip_extnd, cnt_all_reads, irange, cnt_disc, cnt_all_2));
            }
        }
        let mut out = String::new();
        for (ins_pos, rcds) in m_rslts.iter() {
            let mut cnt_clip = 0i64;
            let mut cnt_all = 0i64;
            let mut l_af: Vec<f64> = Vec::new();
            let mut l_af_disc: Vec<f64> = Vec::new();
            let mut l_af_s: Vec<String> = Vec::new();
            let mut l_af_disc_s: Vec<String> = Vec::new();
            for &(c_clip, clip_range, c_all, all_range, c_disc, c_disc_all) in rcds {
                cnt_clip = c_clip;
                cnt_all = c_all;
                if c_all <= 0 {
                    l_af.push(0.0);
                    l_af_disc.push(0.0);
                    l_af_s.push("0.0".into());
                    l_af_disc_s.push("0.0".into());
                } else {
                    let mut af = (c_clip * all_range) as f64 / (c_all * clip_range) as f64;
                    if af > 1.0 {
                        af = 1.0;
                    }
                    l_af.push(af);
                    l_af_s.push(crate::py::fmt::float_repr(af));
                    // ZeroDivisionError in the original when no read is within `extend`
                    let mut af_disc = if c_disc_all == 0 { f64::NAN } else { c_disc as f64 / c_disc_all as f64 };
                    if af_disc > 1.0 {
                        af_disc = 1.0;
                    }
                    l_af_disc.push(af_disc);
                    l_af_disc_s.push(crate::py::fmt::float_repr(af_disc));
                }
            }
            let b_clip_ok = l_af.iter().any(|&x| x <= af_cutoff);
            let b_disc_ok = l_af_disc.iter().any(|&x| x <= af_cutoff);
            if b_clip_ok && b_disc_ok {
                let _ = writeln!(out, "{}\t{}\t{}\t{}\t{}\t{}", chrm, ins_pos, cnt_clip, cnt_all, l_af_s.join("\t"), l_af_disc_s.join("\t"));
            }
        }
        Ok(out)
    }

    /// `calc_AF_by_clip_reads_of_given_list`
    pub fn calc_af_by_clip_reads_of_given_list(&self, sf_bam_list: &str, sf_candidate_list: &str, extnd: i64, clip_slack: i64, af: f64, sf_out: &str) -> Result<()> {
        let mut m_candidate_sites: IndexMap<String, Vec<i64>> = IndexMap::new();
        for (c, p) in Self::load_sites(sf_candidate_list)? {
            m_candidate_sites.entry(c).or_default().push(p);
        }
        let mut l_bams = Vec::new();
        for line in utils::read_lines(sf_bam_list)? {
            let f = split_ws(&line);
            if let Some(b) = f.first() {
                l_bams.push(b.to_string());
            }
        }
        let chrms: Vec<String> = m_candidate_sites.keys().cloned().collect();
        let res: Vec<Result<String>> =
            self.pool().install(|| utils::par_map(chrms.len(), |i| self.calc_af_of_chrm(&chrms[i], &l_bams, &m_candidate_sites, extnd, clip_slack, af)));
        let mut fout = utils::create(sf_out)?;
        for (i, r) in res.into_iter().enumerate() {
            let s = r?;
            // the per-chromosome temporary file is written in the working folder as well
            let sf_tmp_out = format!("{}{}{}", utils::with_slash(&self.working_folder), chrms[i], gv::ALLELE_FREQUENCY_SUFFIX);
            utils::write_string(&sf_tmp_out, &s)?;
            fout.write_all(s.as_bytes())?;
        }
        fout.flush()?;
        Ok(())
    }
}

// ------------------------------------------------------------------------------------------------
// XClipDiscFilter data types

/// m_ins_lpos / m_ins_rpos: chrm -> ins_pos -> ref clip position -> [position on consensus]
pub type ClipPosMap = SiteMap<IndexMap<i64, Vec<i64>>>;

/// Result of `parse_clip_realignment_consensus`.
pub struct ClipRealign {
    pub m_ins_lpos: ClipPosMap,
    pub m_ins_rpos: ClipPosMap,
    pub m_cns_lpos: SiteMap<IndexMap<i64, i64>>,
    pub m_cns_rpos: SiteMap<IndexMap<i64, i64>>,
    pub m_clip_sample: SiteMap<IndexMap<String, [i64; 2]>>,
    /// [n_left, n_right, n_left_polyA (not polyT), n_right_polyA]
    pub m_polya: SiteMap<[i64; 4]>,
}

/// A record of `check_clip_alone_consistency2`:
/// (lpeak, rpeak, lcns_start, lcns_end, rcns_start, rcns_end, b_no_lclip, b_no_rclip, b_consist, b_polyA)
#[derive(Clone, Debug)]
pub struct ClipChk {
    pub lpeak: Option<i64>,
    pub rpeak: Option<i64>,
    pub lc_start: i64,
    pub lc_end: i64,
    pub rc_start: i64,
    pub rc_end: i64,
    pub b_no_lclip: bool,
    pub b_no_rclip: bool,
    pub b_consist: bool,
    pub b_polya: bool,
}

/// A discordant read re-aligned to the consensus:
/// (anchor_pos, is_anchor_lclip, is_anchor_rclip, cns_pos, b_same_ori, b_anchor_rc)
#[derive(Clone, Debug)]
pub struct DiscHit {
    pub anchor_pos: i64,
    pub b_lclip: bool,
    pub b_rclip: bool,
    pub cns_pos: f64,
    pub b_same_ori: bool,
    pub b_anchor_rc: bool,
}

/// Result of `parse_disc_algnmt_consensus`.
pub struct DiscRealign {
    pub m_disc_pos: SiteMap<Vec<DiscHit>>,
    pub m_disc_sample: SiteMap<IndexMap<String, i64>>,
    pub m_polya: SiteMap<i64>,
}

/// A record of `check_seprt_disc_consistency`.
#[derive(Clone, Debug)]
pub struct DiscFlt {
    pub n_l_disc: i64,
    pub lcns_start: PyNum,
    pub lcns_end: PyNum,
    pub n_r_disc: i64,
    pub rcns_start: PyNum,
    pub rcns_end: PyNum,
    pub rc_rcd: [i64; 4],
    pub anchor_rc_rcd: [i64; 4],
}

/// Transduction information attached to a site:
/// `(s_trsdct_info, (b_ltrsdct, n_clip_trsdct, ref_clip_pos_transduct, n_disc_trsdct, n_polyA_trsdct))`
#[derive(Clone, Debug)]
pub struct TransductInfo {
    pub s_info: String,
    pub b_ltrsdct: bool,
    pub n_clip: i64,
    pub ref_clip_pos: i64,
    pub n_disc: i64,
    pub n_polya: i64,
}

/// The 34-field candidate record built by `_extract_detailed_info_TEI`.
#[derive(Clone, Debug)]
pub struct TeiRcd {
    pub ins_chrm: String,
    pub ins_pos: i64,
    pub refined_pos: i64,
    pub lpeak_pos: i64,
    pub rpeak_pos: i64,
    pub tsd: i64,
    pub n_total_lclip: i64,
    pub n_total_rclip: i64,
    pub n_total_ldisc: i64,
    pub n_total_rdisc: i64,
    pub n_total_lpolya: i64,
    pub n_total_rpolya: i64,
    pub flcov: PyNum,
    pub frcov: PyNum,
    pub nlclip: i64,
    pub nrclip: i64,
    pub nldisc: i64,
    pub nrdisc: i64,
    pub nlpolya: i64,
    pub nrpolya: i64,
    pub rc_rcd: [i64; 4],
    pub anchor_rc_rcd: [i64; 4],
    pub lc_cns_start: i64,
    pub lc_cns_end: i64,
    pub rc_cns_start: i64,
    pub rc_cns_end: i64,
    pub ldisc_cns_start: PyNum,
    pub ldisc_cns_end: PyNum,
    pub rdisc_cns_start: PyNum,
    pub rdisc_cns_end: PyNum,
    pub s_trsdct_info: String,
    pub n_clip_trsdct: i64,
    pub n_disc_trsdct: i64,
    pub n_polya_trsdct: i64,
}

impl TeiRcd {
    /// Common part of the dump functions (from "TSD" up to the transduction fields).
    fn fmt_body(&self, s_pos: &str) -> String {
        format!(
            "{}{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}:{}\t{}:{}\t{}:{}\t{}:{}\t{}\t{}\t{}\t{}\t",
            s_pos,
            self.n_total_lclip,
            self.n_total_rclip,
            self.n_total_ldisc,
            self.n_total_rdisc,
            self.n_total_lpolya,
            self.n_total_rpolya,
            self.flcov,
            self.frcov,
            self.nlclip,
            self.nrclip,
            self.nldisc,
            self.nrdisc,
            self.nlpolya,
            self.nrpolya,
            self.lc_cns_start,
            self.lc_cns_end,
            self.rc_cns_start,
            self.rc_cns_end,
            self.ldisc_cns_start,
            self.ldisc_cns_end,
            self.rdisc_cns_start,
            self.rdisc_cns_end,
            self.s_trsdct_info,
            self.n_clip_trsdct,
            self.n_disc_trsdct,
            self.n_polya_trsdct
        )
    }
    fn s_rc_disc(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t",
            self.rc_rcd[0],
            self.rc_rcd[1],
            self.rc_rcd[2],
            self.rc_rcd[3],
            is_3mer_inversion(&self.rc_rcd)
        )
    }
    fn s_anchor_rc(&self) -> String {
        format!("{}\t{}\t{}\t{}", self.anchor_rc_rcd[0], self.anchor_rc_rcd[1], self.anchor_rc_rcd[2], self.anchor_rc_rcd[3])
    }
}

/// `_is_3mer_inversion`
pub fn is_3mer_inversion(rc: &[i64; 4]) -> &'static str {
    let b_l_ff = !(rc[1] > rc[0]);
    let b_r_ff = !(rc[3] > rc[2]);
    if b_l_ff == b_r_ff {
        gv::NOT_FIVE_PRIME_INV
    } else {
        gv::FIVE_PRIME_INVERSION
    }
}

/// Header line used by the dump functions.
fn dump_title() -> String {
    let stitle1 = "#chrm\trefined-pos\tlclip-pos\trclip-pos\tTSD\tnalclip\tnarclip\tnaldisc\t";
    let stitle2 = "nardisc\tnalpolyA\tnarpolyA\tlcov\trcov\tnlclip\tnrclip\tnldisc\tnrdisc\tnlpolyA\tnrpolyA\t";
    let stitle3 = "lclip-cns-start:end\trclip_cns-start:end\tldisc-cns-start:end\trdisc-cns-start:end\t";
    let stitle4 = "Transduction-info\tntsd-clip\tntsd-disc\tntsd-polyA\t";
    let stitle5 = "ldisc-same\tldisc-diff\trdisc-same\trdisc-diff\n";
    format!("{}{}{}{}{}", stitle1, stitle2, stitle3, stitle4, stitle5)
}

#[inline]
fn site_has<T>(m: &SiteMap<T>, chrm: &str, pos: i64) -> bool {
    m.get(chrm).map_or(false, |x| x.contains_key(&pos))
}

#[inline]
fn site_get<'a, T>(m: &'a SiteMap<T>, chrm: &str, pos: i64) -> Option<&'a T> {
    m.get(chrm).and_then(|x| x.get(&pos))
}

/// `is_clipped_part_qualified_algnmt` -> (qualified, n_mapped_bases)
pub fn is_clipped_part_qualified_algnmt(l_cigar: &[(u32, u32)], ratio_cutoff: f64) -> (bool, i64) {
    if l_cigar.is_empty() {
        return (false, 0);
    }
    if l_cigar.len() > 2 {
        let first = l_cigar[0];
        let last = l_cigar[l_cigar.len() - 1];
        let (mut b_left_clip, mut i_left_clip_len) = (false, 0i64);
        if first.0 == BAM_CSOFT_CLIP || first.0 == BAM_CHARD_CLIP {
            b_left_clip = true;
            i_left_clip_len = first.1 as i64;
        }
        let (mut b_right_clip, mut i_right_clip_len) = (false, 0i64);
        if last.0 == BAM_CSOFT_CLIP || last.0 == BAM_CHARD_CLIP {
            b_right_clip = true;
            i_right_clip_len = last.1 as i64;
        }
        if b_left_clip && b_right_clip && i_left_clip_len > gv::MAX_CLIP_CLIP_LEN && i_right_clip_len > gv::MAX_CLIP_CLIP_LEN {
            return (false, 0);
        }
    }
    let mut n_total = 0i64;
    let mut n_map = 0i64;
    for &(t, l) in l_cigar {
        if t == BAM_CMATCH {
            n_map += l as i64;
        }
        if t != BAM_CDEL {
            n_total += l as i64;
        }
    }
    if (n_map as f64) < (n_total as f64) * ratio_cutoff {
        return (false, 0);
    }
    (true, n_map)
}

// ------------------------------------------------------------------------------------------------
// XClipDiscFilter

pub struct XClipDiscFilter {
    pub sf_bam_list: String,
    pub working_folder: String,
    pub n_jobs: usize,
    pub sf_reference: String,
}

/// Category sets of `_first_stage_classify`.
pub type InsCategories = HashMap<&'static str, HashSet<String>>;

impl XClipDiscFilter {
    pub fn new(sf_bam_list: &str, working_folder: &str, n_jobs: usize, sf_reference: &str) -> Self {
        XClipDiscFilter {
            sf_bam_list: sf_bam_list.to_string(),
            working_folder: working_folder.to_string(),
            n_jobs,
            sf_reference: sf_reference.to_string(),
        }
    }

    fn bam_list(&self) -> Result<Vec<String>> {
        let mut v = Vec::new();
        for line in utils::read_lines(&self.sf_bam_list)? {
            let f = split_ws(&line);
            if let Some(b) = f.first() {
                v.push(b.to_string());
            }
        }
        Ok(v)
    }

    /// `_get_n_bams`
    fn get_n_bams(&self) -> Result<i64> {
        let mut n = 0;
        for line in utils::read_lines(&self.sf_bam_list)? {
            if !utils::rstrip(&line).is_empty() {
                n += 1;
            }
        }
        Ok(n)
    }

    /// `cnt_n_bams`
    pub fn cnt_n_bams(sf_bam: &str) -> Result<i64> {
        let mut n = 0;
        for line in utils::read_lines(sf_bam)? {
            if utils::rstrip(&line).len() > 1 {
                n += 1;
            }
        }
        Ok(n)
    }

    pub fn clean_file_by_path(sf_path: &str) {
        if utils::is_file(sf_path) {
            utils::remove_file(sf_path);
        }
    }

    /// `_cnt_clip_reads`
    fn cnt_clip_reads(ins_chrm: &str, ins_pos: i64, peak_clip_pos: i64, m_clip_pos: &ClipPosMap) -> i64 {
        let mut nclip = 0i64;
        match site_get(m_clip_pos, ins_chrm, ins_pos) {
            None => {
                println!("Error in cnt clip step, {}:{} not in the dict!!!", ins_chrm, ins_pos);
                nclip
            }
            Some(m) => {
                for (tmp_pos, l) in m.iter() {
                    if (tmp_pos - peak_clip_pos).abs() <= gv::CLIP_SEARCH_WINDOW {
                        nclip += l.len() as i64;
                    }
                }
                nclip
            }
        }
    }

    /// `_refine_represent_clip_pos`
    pub fn refine_represent_clip_pos(m_list: &SiteMap<ClipChk>) -> SiteMap<i64> {
        let mut m_refined: SiteMap<i64> = IndexMap::new();
        for (ins_chrm, m) in m_list.iter() {
            for (&ins_pos, record) in m.iter() {
                let lpeak_pos = record.lpeak.unwrap_or(-1);
                let rpeak_pos = record.rpeak.unwrap_or(-1);
                let mut refined_pos = lpeak_pos;
                if lpeak_pos == -1 {
                    refined_pos = rpeak_pos;
                }
                if refined_pos == -1 {
                    refined_pos = ins_pos;
                }
                m_refined.entry(ins_chrm.clone()).or_default().insert(ins_pos, refined_pos);
            }
        }
        m_refined
    }

    /// `dump_TEI_info`
    pub fn dump_tei_info(l_candidates: &[TeiRcd], sf_final_list: &str, b_with_head: bool) -> Result<()> {
        let mut fout = utils::create(sf_final_list)?;
        if b_with_head {
            fout.write_all(dump_title().as_bytes())?;
        }
        for r in l_candidates {
            let s_pos = format!("{}\t{}\t{}\t{}\t{}\t", r.ins_chrm, r.refined_pos, r.lpeak_pos, r.rpeak_pos, r.tsd);
            let sinfo = format!("{}{}{}\n", r.fmt_body(&s_pos), r.s_rc_disc(), r.s_anchor_rc());
            fout.write_all(sinfo.as_bytes())?;
        }
        fout.flush()?;
        Ok(())
    }

    /// `dump_TEI_info_with_ori_pos`
    pub fn dump_tei_info_with_ori_pos(l_candidates: &[TeiRcd], sf_final_list: &str, b_with_head: bool) -> Result<()> {
        let mut fout = utils::create(sf_final_list)?;
        if b_with_head {
            fout.write_all(dump_title().as_bytes())?;
        }
        for r in l_candidates {
            let s_pos = format!("{}\t{}\t{}\t{}\t{}\t{}\t", r.ins_chrm, r.ins_pos, r.refined_pos, r.lpeak_pos, r.rpeak_pos, r.tsd);
            let sinfo = format!("{}{}{}\n", r.fmt_body(&s_pos), r.s_rc_disc(), r.s_anchor_rc());
            fout.write_all(sinfo.as_bytes())?;
        }
        fout.flush()?;
        Ok(())
    }

    /// `load_TEI_info_from_file`
    pub fn load_tei_info_from_file(sf_tei: &str) -> Result<SiteMap<(i64, i64, i64, i64, i64, i64, String, String, String, String)>> {
        let mut m: SiteMap<(i64, i64, i64, i64, i64, i64, String, String, String, String)> = IndexMap::new();
        for line in utils::read_lines(sf_tei)? {
            let f = split_ws(&line);
            if f.len() < 24 {
                return Err(anyhow!("IndexError: malformed line in {}", sf_tei));
            }
            let rec = (
                utils::pint(f[6])?,
                utils::pint(f[7])?,
                utils::pint(f[8])?,
                utils::pint(f[9])?,
                utils::pint(f[10])?,
                utils::pint(f[11])?,
                f[20].to_string(),
                f[21].to_string(),
                f[22].to_string(),
                f[23].to_string(),
            );
            m.entry(f[0].to_string()).or_default().insert(utils::pint(f[1])?, rec);
        }
        Ok(m)
    }

    /// `sprt_TEI_to_td_non_td`
    pub fn sprt_tei_to_td_non_td(sf_cns: &str, sf_td: &str, sf_non_td: &str) -> Result<()> {
        let mut fout_td = utils::create(sf_td)?;
        let mut fout_non_td = utils::create(sf_non_td)?;
        for line in utils::read_raw_lines(sf_cns)? {
            let fields = split_ws(&line);
            // NOTE: `len(fields) >= 23` then `fields[23]` as in the original
            if fields.len() >= 23 && fields.get(23).copied() == Some(gv::NOT_TRANSDUCTION) {
                if fields[23].split(':').count() > 1 {
                    fout_td.write_all(line.as_bytes())?;
                } else {
                    fout_non_td.write_all(line.as_bytes())?;
                }
            } else {
                fout_td.write_all(line.as_bytes())?;
            }
        }
        Ok(())
    }

    /// `sprt_TEI_to_td_orphan_non_td`
    pub fn sprt_tei_to_td_orphan_non_td(sf_cns: &str, sf_non_td: &str, sf_td: &str, sf_td_sibling: &str, sf_orphan: &str) -> Result<()> {
        let mut fout_td = utils::create(sf_td)?;
        let mut fout_orphan = utils::create(sf_orphan)?;
        let mut fout_sibling = utils::create(sf_td_sibling)?;
        let mut fout_non_td = utils::create(sf_non_td)?;
        for line in utils::read_raw_lines(sf_cns)? {
            let fields = split_ws(&line);
            if line.contains("orphan") || (fields.len() >= 32 && fields.get(32).copied() == Some(gv::ORPHAN_TRANSDUCTION)) {
                fout_orphan.write_all(line.as_bytes())?;
                fout_td.write_all(line.as_bytes())?;
            } else if fields.len() >= 23 && fields.get(23).copied() == Some(gv::NOT_TRANSDUCTION) {
                if fields[23].split(':').count() > 1 {
                    fout_td.write_all(line.as_bytes())?;
                } else {
                    fout_non_td.write_all(line.as_bytes())?;
                }
            } else {
                if line.contains(gv::SIBLING_LABEL) {
                    fout_sibling.write_all(line.as_bytes())?;
                }
                fout_td.write_all(line.as_bytes())?;
            }
        }
        Ok(())
    }

    /// `_get_insertion_length`
    fn get_insertion_length(ldisc_cns_start: PyNum, ldisc_cns_end: PyNum, rdisc_cns_start: PyNum, rdisc_cns_end: PyNum, icns_lth: i64, stsdc: &str) -> PyNum {
        let sub = |a: i64, b: PyNum| -> PyNum {
            match b {
                PyNum::Int(x) => PyNum::Int(a - x),
                PyNum::Float(x) => PyNum::Float(a as f64 - x),
            }
        };
        let sub2 = |a: PyNum, b: PyNum| -> PyNum {
            match (a, b) {
                (PyNum::Int(x), PyNum::Int(y)) => PyNum::Int(x - y),
                _ => PyNum::Float(a.f() - b.f()),
            }
        };
        let mut i_tei_len = PyNum::Int(-1);
        if ldisc_cns_end.f() == -1.0 && rdisc_cns_start.f() > 0.0 {
            i_tei_len = sub(icns_lth, rdisc_cns_start);
        } else if rdisc_cns_end.f() == -1.0 && ldisc_cns_start.f() > 0.0 {
            i_tei_len = sub(icns_lth, ldisc_cns_start);
        } else if rdisc_cns_end.f() < 0.0 && ldisc_cns_start.f() < 0.0 {
            i_tei_len = PyNum::Int(-1);
        } else if stsdc == gv::NOT_TRANSDUCTION {
            if ldisc_cns_start.f() > 0.0 && rdisc_cns_start.f() > 0.0 {
                i_tei_len = sub2(rdisc_cns_end, ldisc_cns_start);
                if i_tei_len.f() < 0.0 {
                    i_tei_len = sub2(ldisc_cns_end, rdisc_cns_start);
                }
            }
        } else if ldisc_cns_start.f() < rdisc_cns_start.f() {
            i_tei_len = sub(icns_lth, ldisc_cns_start);
        } else {
            i_tei_len = sub(icns_lth, rdisc_cns_start);
        }
        i_tei_len
    }

    /// `_second_stage_classify_dump`
    #[allow(clippy::too_many_arguments)]
    pub fn second_stage_classify_dump(
        m_ins_categories: &InsCategories,
        l_candidates: &[TeiRcd],
        m_high_confident: &SiteMap<()>,
        m_one_side_low_confident: &SiteMap<()>,
        _m_hcov: &SiteMap<[i64; 2]>,
        m_not_hit_end: &SiteMap<()>,
        m_gntp_features: &SiteMap<GntpFeat>,
        i_cns_lth: i64,
        b_with_head: bool,
        sf_final_list: &str,
    ) -> Result<()> {
        let sf_final_high_confident = format!("{}{}", sf_final_list, gv::HIGH_CONFIDENT_SUFFIX);
        let mut fout_final = utils::create(sf_final_list)?;
        let mut fout_hc = utils::create(&sf_final_high_confident)?;
        if b_with_head {
            let stitle = format!(
                "{}{}{}{}{}{}",
                "#chrm\trefined-pos\tlclip-pos\trclip-pos\tTSD\tnalclip\tnarclip\tnaldisc\t",
                "nardisc\tnalpolyA\tnarpolyA\tlcov\trcov\tnlclip\tnrclip\tnldisc\tnrdisc\tnlpolyA\tnrpolyA\t",
                "lclip-cns-start:end\trclip_cns-start:end\tldisc-cns-start:end\trdisc-cns-start:end\t",
                "Transduction-info\tntsd-clip\tntsd-disc\tntsd-polyA\t",
                "ldisc-same\tldisc-diff\trdisc-same\trdisc-diff\t3mer-inversion\t",
                "confidential\tinsertion-length\n"
            );
            fout_final.write_all(stitle.as_bytes())?;
            fout_hc.write_all(stitle.as_bytes())?;
        }
        let empty: HashSet<String> = HashSet::new();
        let cat = |k: &str| -> &HashSet<String> { m_ins_categories.get(k).unwrap_or(&empty) };
        for r in l_candidates {
            let ins_chrm = r.ins_chrm.as_str();
            let ins_pos = r.ins_pos;
            let s_pos = format!("{}\t{}\t{}\t{}\t{}\t", ins_chrm, r.refined_pos, r.lpeak_pos, r.rpeak_pos, r.tsd);
            let mut s_trsdct_info = r.s_trsdct_info.clone();
            let s_anchor_rc = format!("\t{}", r.s_anchor_rc());
            let mut b_polya_dominant = false;
            let s_hit_end = if site_has(m_not_hit_end, ins_chrm, ins_pos) {
                format!("\t{}", gv::NOT_HIT_END_OF_CNS)
            } else {
                format!("\t{}", gv::HIT_END_OF_CNS)
            };
            let s_both_consistent = if site_has(m_high_confident, ins_chrm, ins_pos) {
                format!("\t{}", gv::BOTH_END_CONSISTNT)
            } else {
                format!("\t{}", gv::ONE_END_CONSISTNT)
            };
            let g = match site_get(m_gntp_features, ins_chrm, ins_pos) {
                Some(g) => g,
                None => continue,
            };
            let n_raw_lclip = g.n_l_raw_clip;
            let n_raw_rclip = g.n_r_raw_clip;
            let sf_gntp_info = format!(
                "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                g.n_af_clip, g.n_full_map, g.n_l_raw_clip, g.n_r_raw_clip, g.n_disc_pairs, g.n_concd_pairs, g.n_disc_large_indel, g.s_clip_lens
            );
            if n_raw_lclip + n_raw_rclip > 0
                && ((r.n_total_lpolya + r.n_total_rpolya) as f64) / ((n_raw_lclip + n_raw_rclip) as f64) > gv::MAX_POLYA_RATIO
            {
                b_polya_dominant = true;
            }
            let s_tmp_id = format!("{}{}{}", ins_chrm, gv::SEPERATOR, ins_pos);
            let mut b_high_confident = false;
            let mut s_extra_info: &str;
            if cat(gv::TWO_SIDE_TPRT_BOTH).contains(&s_tmp_id) {
                s_extra_info = if b_polya_dominant { gv::TWO_SIDE_POLYA_DOMINANT } else { gv::TWO_SIDE_TPRT_BOTH };
                b_high_confident = true;
            } else if cat(gv::TWO_SIDE_TPRT).contains(&s_tmp_id) {
                s_extra_info = if b_polya_dominant { gv::TWO_SIDE_POLYA_DOMINANT } else { gv::TWO_SIDE_TPRT };
                b_high_confident = true;
            } else if cat(gv::TWO_SIDE).contains(&s_tmp_id) {
                b_high_confident = true;
                if b_polya_dominant {
                    s_extra_info = gv::TWO_SIDE_POLYA_DOMINANT;
                    b_high_confident = false;
                } else {
                    s_extra_info = gv::TWO_SIDE;
                }
            } else if s_trsdct_info != gv::NOT_TRANSDUCTION {
                s_extra_info = gv::ONE_SIDE_TRSDCT;
            } else if cat(gv::ONE_HALF_SIDE).contains(&s_tmp_id) {
                s_extra_info = gv::ONE_HALF_SIDE;
                b_high_confident = true;
                if b_polya_dominant {
                    s_extra_info = gv::ONE_HALF_SIDE_POLYA_DOMINANT;
                    b_high_confident = false;
                } else if cat(gv::ONE_HALF_SIDE_TRPT_BOTH).contains(&s_tmp_id) {
                    s_extra_info = gv::ONE_HALF_SIDE_TRPT_BOTH;
                } else if cat(gv::ONE_HALF_SIDE_TRPT).contains(&s_tmp_id) {
                    s_extra_info = gv::ONE_HALF_SIDE_TRPT;
                }
            } else if cat(gv::ONE_SIDE).contains(&s_tmp_id) {
                s_extra_info = gv::ONE_SIDE;
                if site_has(m_one_side_low_confident, ins_chrm, ins_pos) {
                    s_extra_info = gv::ONE_SIDE_COVERAGE_CONFLICT;
                }
            } else if cat(gv::HIGH_COV_ISD).contains(&s_tmp_id) {
                s_extra_info = gv::HIGH_COV_ISD;
            } else {
                s_extra_info = gv::OTHER_TYPE;
            }
            if s_trsdct_info == gv::ONE_SIDE_FLANKING {
                if gv::is_call_sva() {
                    b_high_confident = true;
                }
                s_extra_info = gv::ONE_SIDE_FLANKING;
                s_trsdct_info = gv::NOT_TRANSDUCTION.to_string();
            }
            let i_tei_lth =
                Self::get_insertion_length(r.ldisc_cns_start, r.ldisc_cns_end, r.rdisc_cns_start, r.rdisc_cns_end, i_cns_lth, &s_trsdct_info);
            let s_lenth_info = format!("\t{}", i_tei_lth);
            // s_transduct_out is built from the record's original transduction info
            let sinfo = format!(
                "{}{}{}{}{}{}{}{}\n",
                r.fmt_body(&s_pos),
                r.s_rc_disc(),
                s_extra_info,
                s_both_consistent,
                s_hit_end,
                sf_gntp_info,
                s_anchor_rc,
                s_lenth_info
            );
            fout_final.write_all(sinfo.as_bytes())?;
            if b_high_confident {
                fout_hc.write_all(sinfo.as_bytes())?;
            }
        }
        fout_final.flush()?;
        fout_hc.flush()?;
        Ok(())
    }

    /// `_extract_detailed_info_TEI`
    #[allow(clippy::too_many_arguments)]
    pub fn extract_detailed_info_tei(
        m_disc_filtered: &SiteMap<DiscFlt>,
        m_disc_bf_filter: &SiteMap<ClipChk>,
        m_clip_lpos: &ClipPosMap,
        m_clip_rpos: &ClipPosMap,
        m_polya: &SiteMap<[i64; 4]>,
        m_clip_checked_list: &SiteMap<ClipChk>,
        m_transduct_info: &SiteMap<TransductInfo>,
        m_rescued_n_src: &SiteMap<(i64, i64, i64)>,
        m_depth: Option<&SiteMap<[PyNum; 4]>>,
    ) -> Result<Vec<TeiRcd>> {
        let mut l_candidates = Vec::new();
        for (ins_chrm, m) in m_disc_filtered.iter() {
            for (&ins_pos, disc_record) in m.iter() {
                let mut nlclip = 0i64;
                let mut nrclip = 0i64;
                let mut nlpolya = 0i64;
                let mut nrpolya = 0i64;
                let mut n_total_lclip = 0i64;
                let mut n_total_rclip = 0i64;
                let mut n_total_ldisc = 0i64;
                let mut n_total_rdisc = 0i64;
                let mut n_total_lpolya = 0i64;
                let mut n_total_rpolya = 0i64;
                let record = site_get(m_disc_bf_filter, ins_chrm, ins_pos).ok_or_else(|| anyhow!("KeyError: {}:{}", ins_chrm, ins_pos))?;
                let mut lpeak_pos = record.lpeak.unwrap_or(-1);
                let mut rpeak_pos = record.rpeak.unwrap_or(-1);
                if lpeak_pos > 0 {
                    nlclip = Self::cnt_clip_reads(ins_chrm, ins_pos, lpeak_pos, m_clip_lpos);
                    n_total_lclip += nlclip;
                }
                if rpeak_pos > 0 {
                    nrclip = Self::cnt_clip_reads(ins_chrm, ins_pos, rpeak_pos, m_clip_rpos);
                    n_total_rclip += nrclip;
                }
                if let Some(p) = site_get(m_polya, ins_chrm, ins_pos) {
                    nlpolya = p[0];
                    nrpolya = p[1];
                    n_total_lpolya += nlpolya;
                    n_total_rpolya += nrpolya;
                    let nlpolya_seq_a = p[2];
                    let nrpolya_seq_a = p[3];
                    if nlpolya_seq_a < nlpolya - nlpolya_seq_a {
                        nlpolya = -nlpolya;
                    }
                    if nrpolya_seq_a < nrpolya - nrpolya_seq_a {
                        nrpolya = -nrpolya;
                    }
                }
                let (mut lc_cns_start, mut lc_cns_end, mut rc_cns_start, mut rc_cns_end) = (-1i64, -1i64, -1i64, -1i64);
                if let Some(c) = site_get(m_clip_checked_list, ins_chrm, ins_pos) {
                    lc_cns_start = c.lc_start;
                    lc_cns_end = c.lc_end;
                    rc_cns_start = c.rc_start;
                    rc_cns_end = c.rc_end;
                }
                let mut n_clip_trsdct = 0i64;
                let mut n_disc_trsdct = 0i64;
                let mut n_polya_trsdct = 0i64;
                let mut s_trsdct_info = gv::NOT_TRANSDUCTION.to_string();
                if let Some(t) = site_get(m_transduct_info, ins_chrm, ins_pos) {
                    s_trsdct_info = t.s_info.clone();
                    n_clip_trsdct = t.n_clip;
                    n_disc_trsdct = t.n_disc;
                    n_polya_trsdct = t.n_polya;
                    if t.b_ltrsdct {
                        if t.n_clip > nlclip {
                            lpeak_pos = t.ref_clip_pos;
                        }
                        n_total_lclip += t.n_clip;
                        n_total_rdisc += t.n_disc;
                        n_total_lpolya += t.n_polya;
                    } else {
                        if t.n_clip > nrclip {
                            rpeak_pos = t.ref_clip_pos;
                        }
                        n_total_rclip += t.n_clip;
                        n_total_ldisc += t.n_disc;
                        n_total_rpolya += t.n_polya;
                    }
                }
                if let Some(&(n_rsc_lclip, n_rsc_rclip, n_rsc_disc)) = site_get(m_rescued_n_src, ins_chrm, ins_pos) {
                    s_trsdct_info = gv::ONE_SIDE_FLANKING.to_string();
                    n_disc_trsdct = n_rsc_disc;
                    if n_rsc_lclip > n_rsc_rclip {
                        n_clip_trsdct = n_rsc_lclip;
                        n_total_lclip += n_rsc_lclip;
                        n_total_rdisc += n_rsc_disc;
                    } else {
                        n_clip_trsdct = n_rsc_rclip;
                        n_total_rclip += n_rsc_lclip;
                        n_total_ldisc += n_rsc_disc;
                    }
                }
                let mut tsd = -1i64;
                if lpeak_pos > 0 && rpeak_pos > 0 {
                    tsd = (rpeak_pos - lpeak_pos).abs();
                }
                let mut refined_pos = lpeak_pos;
                if lpeak_pos == -1 || n_total_rclip > n_total_lclip {
                    refined_pos = rpeak_pos;
                }
                let nldisc = disc_record.n_l_disc;
                n_total_ldisc += nldisc;
                let nrdisc = disc_record.n_r_disc;
                n_total_rdisc += nrdisc;
                let mut flcov = PyNum::Float(0.0);
                let mut frcov = PyNum::Float(0.0);
                if let Some(md) = m_depth {
                    if let Some(d) = site_get(md, ins_chrm, ins_pos) {
                        flcov = d[0];
                        frcov = d[1];
                    }
                }
                l_candidates.push(TeiRcd {
                    ins_chrm: ins_chrm.clone(),
                    ins_pos,
                    refined_pos,
                    lpeak_pos,
                    rpeak_pos,
                    tsd,
                    n_total_lclip,
                    n_total_rclip,
                    n_total_ldisc,
                    n_total_rdisc,
                    n_total_lpolya,
                    n_total_rpolya,
                    flcov,
                    frcov,
                    nlclip,
                    nrclip,
                    nldisc,
                    nrdisc,
                    nlpolya,
                    nrpolya,
                    rc_rcd: disc_record.rc_rcd,
                    anchor_rc_rcd: disc_record.anchor_rc_rcd,
                    lc_cns_start,
                    lc_cns_end,
                    rc_cns_start,
                    rc_cns_end,
                    ldisc_cns_start: disc_record.lcns_start,
                    ldisc_cns_end: disc_record.lcns_end,
                    rdisc_cns_start: disc_record.rcns_start,
                    rdisc_cns_end: disc_record.rcns_end,
                    s_trsdct_info,
                    n_clip_trsdct,
                    n_disc_trsdct,
                    n_polya_trsdct,
                });
            }
        }
        Ok(l_candidates)
    }

    /// `_is_distance_consistency`
    fn is_distance_consistency(lstart: f64, lend: f64, rstart: f64, rend: f64, i_is: i64) -> bool {
        let idist = i_is as f64 / 2.0 + 1.0;
        let lstart = lstart - idist;
        let lend = lend + idist;
        let rstart = rstart - idist;
        let rend = rend + idist;
        (lstart >= rstart && lstart <= rend) || (lend >= rstart && lend <= rend)
    }

    /// `_is_depth_differ`
    fn is_depth_differ(flcov: f64, frcov: f64, b_left: bool, ratio: f64) -> bool {
        if b_left {
            if frcov <= 0.0 {
                false
            } else {
                (frcov - flcov) / frcov >= ratio
            }
        } else if flcov <= 0.0 {
            false
        } else {
            (flcov - frcov) / flcov >= ratio
        }
    }

    /// `filter_by_clip_disc_polyA_consistency` -> (selected, m_high_confident(by refined pos),
    /// m_one_side_low_confident, m_not_hit_end).  The coverage fields of the records are updated.
    #[allow(clippy::too_many_arguments)]
    pub fn filter_by_clip_disc_polya_consistency(
        l_candidates: &mut [TeiRcd],
        m_coverage: &SiteMap<[PyNum; 4]>,
        i_is: i64,
        nclip_cutoff: f64,
        ndisc_cutoff: f64,
        min_cns_end: i64,
        depth_ratio: f64,
        xlog: &mut XLog,
    ) -> Result<(Vec<TeiRcd>, SiteMap<()>, SiteMap<()>, SiteMap<()>)> {
        let mut l_selected = Vec::new();
        let mut m_high_confident: SiteMap<()> = IndexMap::new();
        let mut m_one_side_low_confident: SiteMap<()> = IndexMap::new();
        let mut m_not_hit_end: SiteMap<()> = IndexMap::new();
        let min_cns_end = min_cns_end as f64;
        for record in l_candidates.iter_mut() {
            let ins_chrm = record.ins_chrm.clone();
            let ins_pos = record.ins_pos;
            let refined_pos = record.refined_pos;
            let nlclip = record.n_total_lclip as f64;
            let nrclip = record.n_total_rclip as f64;
            let nldisc = record.n_total_ldisc as f64;
            let nrdisc = record.n_total_rdisc as f64;
            let nlpolya = record.n_total_lpolya;
            let nrpolya = record.n_total_rpolya;
            let ilclip_start = record.lc_cns_start as f64;
            let ilclip_end = record.lc_cns_end as f64;
            let irclip_start = record.rc_cns_start as f64;
            let irclip_end = record.rc_cns_end as f64;
            let ildisc_start = record.ldisc_cns_start.f();
            let ildisc_end = record.ldisc_cns_end.f();
            let irdisc_start = record.rdisc_cns_start.f();
            let irdisc_end = record.rdisc_cns_end.f();
            let cov = site_get(m_coverage, &ins_chrm, ins_pos).ok_or_else(|| anyhow!("KeyError: {}:{} not in coverage", ins_chrm, ins_pos))?;
            let flcov = cov[0];
            let frcov = cov[1];
            record.flcov = flcov;
            record.frcov = frcov;
            let (flcov, frcov) = (flcov.f(), frcov.f());
            let mut bl_depth_differ = Self::is_depth_differ(flcov, frcov, true, depth_ratio);
            let mut br_depth_differ = Self::is_depth_differ(flcov, frcov, false, depth_ratio);
            let b_transduction = record.s_trsdct_info != gv::NOT_TRANSDUCTION;
            if b_transduction {
                if bl_depth_differ || br_depth_differ {
                    xlog.append(&format!(
                        "{}:{} is filtered out: it's called as candidate transduction, but left and right coverage are different!\n",
                        ins_chrm, ins_pos
                    ));
                    continue;
                } else {
                    l_selected.push(record.clone());
                    xlog.append(&format!("{}:{} is selected as transduction (or one side in flanking region) event!\n", ins_chrm, ins_pos));
                    continue;
                }
            }
            let bl_partial = nldisc < ndisc_cutoff && nrclip <= 0.0;
            let br_partial = nrdisc < ndisc_cutoff && nlclip <= 0.0;
            if bl_partial && !bl_depth_differ {
                xlog.append(&format!(
                    "{}:{} is filtered out: no left side feature, but coverage for left and right are similar!\n",
                    ins_chrm, ins_pos
                ));
                m_one_side_low_confident.entry(ins_chrm.clone()).or_default().insert(ins_pos, ());
            }
            if br_partial && !br_depth_differ {
                xlog.append(&format!(
                    "{}:{} is filtered out: no right side feature, but coverage for left and right are similar!\n",
                    ins_chrm, ins_pos
                ));
                m_one_side_low_confident.entry(ins_chrm.clone()).or_default().insert(ins_pos, ());
            }
            let mut bl_consist = false;
            let mut br_consist = false;
            let mut b_both_end_consist = false;
            if nlclip >= nclip_cutoff && nrdisc >= ndisc_cutoff {
                bl_consist = Self::is_distance_consistency(ilclip_start, ilclip_end, irdisc_start, irdisc_end, i_is);
            }
            if nrclip >= nclip_cutoff && nldisc >= ndisc_cutoff {
                br_consist = Self::is_distance_consistency(irclip_start, irclip_end, ildisc_start, ildisc_end, i_is);
            }
            if !bl_consist && !br_consist && !gv::is_call_sva() {
                xlog.append(&format!("{}:{} is filtered out: Neither of the sides are consistent!\n", ins_chrm, ins_pos));
                continue;
            } else if bl_consist && br_consist {
                b_both_end_consist = true;
            }
            if (nlclip >= nclip_cutoff && nrdisc < ndisc_cutoff) && (nrclip >= nclip_cutoff && nldisc < ndisc_cutoff) {
                xlog.append(&format!(
                    "{}:{} is filtered out: Have both side clipped reads, but don't have enough discordant support!\n",
                    ins_chrm, ins_pos
                ));
                continue;
            }
            let b_tailed_disc = ildisc_end >= min_cns_end || irdisc_end >= min_cns_end;
            let b_tailed_clip = ilclip_end >= min_cns_end || irclip_end >= min_cns_end;
            let b_partial_disc = nldisc < ndisc_cutoff || nrdisc < ndisc_cutoff;
            if nlpolya == 0 && nrpolya == 0 && gv::rna_mediated() {
                b_both_end_consist = false;
                if nlclip >= nclip_cutoff && nrclip >= nclip_cutoff && (!b_tailed_disc && !b_tailed_clip) {
                    xlog.append(&format!("{}:{} is filtered out: Have both side clipped reads, but no polyA!\n", ins_chrm, ins_pos));
                    continue;
                }
            }
            if !bl_depth_differ && nrclip < nclip_cutoff {
                bl_depth_differ = Self::is_depth_differ(flcov, frcov, true, depth_ratio);
            }
            if !br_depth_differ && nlclip < nclip_cutoff {
                br_depth_differ = Self::is_depth_differ(flcov, frcov, false, depth_ratio);
            }
            let b_partial_clip = (nlclip < nclip_cutoff && br_depth_differ) || (nrclip < nclip_cutoff && bl_depth_differ);
            if !b_tailed_disc && !b_partial_disc && !b_tailed_clip && !b_partial_clip {
                xlog.append(&format!(
                    "{}:{} is filtered out: {} {} {} {}!\n",
                    ins_chrm, ins_pos, min_cns_end as i64, record.ldisc_cns_end, record.rdisc_cns_end, record.lc_cns_end
                ));
                xlog.append(&format!("{}:{} is filtered out: Disc or clip cluster doesn't reach the end of cns!\n", ins_chrm, ins_pos));
                continue;
            }
            l_selected.push(record.clone());
            if !b_tailed_clip && !b_tailed_disc {
                m_not_hit_end.entry(ins_chrm.clone()).or_default().insert(ins_pos, ());
            }
            if b_both_end_consist {
                m_high_confident.entry(ins_chrm.clone()).or_default().insert(refined_pos, ());
            }
        }
        Ok((l_selected, m_high_confident, m_one_side_low_confident, m_not_hit_end))
    }

    /// `_first_stage_classify`
    pub fn first_stage_classify(l_candidates: &[TeiRcd], m_hcov: &SiteMap<[i64; 2]>) -> InsCategories {
        let mut m: InsCategories = HashMap::new();
        for k in [
            gv::TWO_SIDE,
            gv::TWO_SIDE_TPRT,
            gv::TWO_SIDE_TPRT_BOTH,
            gv::ONE_HALF_SIDE,
            gv::ONE_HALF_SIDE_TRPT_BOTH,
            gv::ONE_HALF_SIDE_TRPT,
            gv::ONE_SIDE,
            gv::ONE_SIDE_TRSDCT,
            gv::ONE_SIDE_WEAK,
            gv::ONE_SIDE_OTHER,
            gv::ONE_SIDE_SV,
            gv::HIGH_COV_ISD,
            gv::OTHER_TYPE,
            gv::TWO_SIDE_POLYA_DOMINANT,
        ] {
            m.insert(k, HashSet::new());
        }
        let hit = |a: f64, b: f64| -> (bool, bool) {
            if a != -1.0 && b != -1.0 {
                (true, false)
            } else if a != -1.0 || b != -1.0 {
                (false, true)
            } else {
                (false, false)
            }
        };
        for r in l_candidates {
            let s_id = format!("{}{}{}", r.ins_chrm, gv::SEPERATOR, r.ins_pos);
            if site_has(m_hcov, &r.ins_chrm, r.ins_pos) {
                m.get_mut(gv::HIGH_COV_ISD).unwrap().insert(s_id);
                continue;
            }
            let (b_lc_hit, _) = hit(r.lc_cns_start as f64, r.lc_cns_end as f64);
            let (b_rc_hit, _) = hit(r.rc_cns_start as f64, r.rc_cns_end as f64);
            let (b_ld_hit, _) = hit(r.ldisc_cns_start.f(), r.ldisc_cns_end.f());
            let (b_rd_hit, _) = hit(r.rdisc_cns_start.f(), r.rdisc_cns_end.f());
            let has_polya = (r.n_total_rpolya + r.n_total_lpolya) > 0;
            let has_tsd = r.tsd > 0 && r.tsd < gv::TSD_CUTOFF;
            if b_lc_hit && b_rc_hit && b_ld_hit && b_rd_hit {
                m.get_mut(gv::TWO_SIDE).unwrap().insert(s_id.clone());
                if has_polya && has_tsd {
                    m.get_mut(gv::TWO_SIDE_TPRT_BOTH).unwrap().insert(s_id);
                } else if has_polya || has_tsd {
                    m.get_mut(gv::TWO_SIDE_TPRT).unwrap().insert(s_id);
                }
            } else if (b_lc_hit && b_rd_hit) && (!b_rc_hit && !b_ld_hit) {
                m.get_mut(gv::ONE_SIDE).unwrap().insert(s_id);
            } else if (b_rc_hit && b_ld_hit) && (!b_lc_hit && !b_rd_hit) {
                m.get_mut(gv::ONE_SIDE).unwrap().insert(s_id);
            } else if (b_lc_hit && b_rd_hit) && (b_rc_hit || b_ld_hit) || (b_rc_hit && b_ld_hit) && (b_lc_hit || b_rd_hit) {
                m.get_mut(gv::ONE_HALF_SIDE).unwrap().insert(s_id.clone());
                if has_polya && has_tsd {
                    m.get_mut(gv::ONE_HALF_SIDE_TRPT_BOTH).unwrap().insert(s_id);
                } else if has_polya || has_tsd {
                    m.get_mut(gv::ONE_HALF_SIDE_TRPT).unwrap().insert(s_id);
                }
            } else {
                m.get_mut(gv::OTHER_TYPE).unwrap().insert(s_id);
            }
        }
        m
    }

    /// `_get_min_cns_length`
    pub fn get_min_cns_length(sf_cns: &str) -> Result<i64> {
        let mut min_len: i64 = 100000000000000000;
        let mut i_acm_len: i64 = 0;
        for line in utils::read_raw_lines(sf_cns)? {
            if line.starts_with('>') {
                if i_acm_len != 0 && min_len > i_acm_len {
                    min_len = i_acm_len;
                }
                i_acm_len = 0;
                continue;
            }
            i_acm_len += utils::rstrip(&line).len() as i64;
        }
        if min_len > i_acm_len {
            min_len = i_acm_len;
        }
        Ok(min_len)
    }

    /// `check_clip_rep_consistency`
    fn check_clip_rep_consistency(mpos: &IndexMap<i64, Vec<i64>>, l_ngbr_pos: &[i64], imax_dist: i64, ratio_cutoff: f64) -> (bool, i64, i64) {
        let mut m_cns_pos: HashMap<i64, i64> = HashMap::new();
        for p in l_ngbr_pos {
            if let Some(l) = mpos.get(p) {
                for &c in l {
                    *m_cns_pos.entry(c).or_insert(0) += 1;
                }
            }
        }
        let mut peak_win_start = 0i64;
        let mut peak_win_end = 0i64;
        if m_cns_pos.is_empty() {
            return (false, peak_win_start, peak_win_end);
        }
        let mut l_pos: Vec<i64> = m_cns_pos.keys().cloned().collect();
        l_pos.sort();
        let mut peak_win_cnt = 0i64;
        let mut cur_win_cnt = 0i64;
        let mut cur_win_start = l_pos[0];
        let mut cur_win_end = 0i64;
        let mut pre_pos = 0i64;
        for &cur_pos in &l_pos {
            if cur_pos - pre_pos > imax_dist && pre_pos != 0 {
                if cur_win_cnt > peak_win_cnt {
                    peak_win_cnt = cur_win_cnt;
                    peak_win_start = cur_win_start;
                    peak_win_end = cur_win_end;
                }
                cur_win_start = cur_pos;
                cur_win_cnt = 0;
            }
            cur_win_end = cur_pos;
            cur_win_cnt += m_cns_pos[&cur_pos];
            pre_pos = cur_pos;
        }
        if cur_win_cnt > peak_win_cnt {
            peak_win_start = cur_win_start;
            peak_win_end = cur_win_end;
        }
        let n_hit = l_pos.iter().filter(|&&p| p >= peak_win_start && p <= peak_win_end).count();
        let ok = (n_hit as f64) / (l_pos.len() as f64) >= ratio_cutoff;
        (ok, peak_win_start, peak_win_end)
    }

    /// `check_one_side_clip_consistency` -> (peak_pos, cns_start, cns_end, b_consist, n_clip)
    fn check_one_side_clip_consistency(m_clip_pos: &ClipPosMap, idist: i64, ratio: f64, bsize: i64) -> SiteMap<(i64, i64, i64, bool, i64)> {
        let mut m_checked: SiteMap<(i64, i64, i64, bool, i64)> = IndexMap::new();
        for (ins_chrm, m) in m_clip_pos.iter() {
            for (&ins_pos, mpos) in m.iter() {
                let cnt: IndexMap<i64, usize> = mpos.iter().map(|(k, v)| (*k, v.len())).collect();
                let (l1_pos, l1_nbr, l1_acm, _, _, _) = find_first_second_peak(&cnt, bsize);
                let (b_lconsist, s, e) = Self::check_clip_rep_consistency(mpos, &l1_nbr, idist, ratio);
                m_checked.entry(ins_chrm.clone()).or_default().insert(ins_pos, (l1_pos, s, e, b_lconsist, l1_acm));
            }
        }
        m_checked
    }

    /// `_is_high_lowq_clip`
    fn is_high_lowq_clip(n_lowq_clip: i64, n_all_clip: i64, f_ratio: f64) -> bool {
        if n_all_clip == 0 {
            return true;
        }
        (n_lowq_clip as f64) / (n_all_clip as f64) > f_ratio
    }

    /// `check_clip_alone_consistency2`
    pub fn check_clip_alone_consistency2(
        m_left_clip_pos: &ClipPosMap,
        m_right_clip_pos: &ClipPosMap,
        idist: i64,
        ratio: f64,
        n_cutoff: f64,
        bsize: i64,
        m_lowq_clip: Option<&SiteMap<i64>>,
    ) -> SiteMap<ClipChk> {
        let mut m_checked: SiteMap<ClipChk> = IndexMap::new();
        let m_left = Self::check_one_side_clip_consistency(m_left_clip_pos, idist, ratio, bsize);
        let m_right = Self::check_one_side_clip_consistency(m_right_clip_pos, idist, ratio, bsize);
        let lowq_skip = |chrm: &str, pos: i64, cnt_clip: i64| -> bool {
            if let Some(mq) = m_lowq_clip {
                if let Some(&n_lowq) = site_get(mq, chrm, pos) {
                    return Self::is_high_lowq_clip(n_lowq, n_lowq + cnt_clip, gv::MAX_LOWQ_CLIP_RATIO);
                }
            }
            false
        };
        for (ins_chrm, m) in m_left.iter() {
            for (&ins_pos, lrecord) in m.iter() {
                let mut cnt_clip = lrecord.4;
                if let Some(rrecord) = site_get(&m_right, ins_chrm, ins_pos) {
                    cnt_clip += rrecord.4;
                    if (cnt_clip as f64) < n_cutoff {
                        continue;
                    }
                    if lowq_skip(ins_chrm, ins_pos, cnt_clip) {
                        continue;
                    }
                    m_checked.entry(ins_chrm.clone()).or_default().insert(
                        ins_pos,
                        ClipChk {
                            lpeak: Some(lrecord.0),
                            rpeak: Some(rrecord.0),
                            lc_start: lrecord.1,
                            lc_end: lrecord.2,
                            rc_start: rrecord.1,
                            rc_end: rrecord.2,
                            b_no_lclip: false,
                            b_no_rclip: false,
                            b_consist: true,
                            b_polya: false,
                        },
                    );
                } else {
                    if (cnt_clip as f64) < n_cutoff / 2.0 + 1.0 {
                        continue;
                    }
                    if lowq_skip(ins_chrm, ins_pos, cnt_clip) {
                        continue;
                    }
                    m_checked.entry(ins_chrm.clone()).or_default().insert(
                        ins_pos,
                        ClipChk {
                            lpeak: Some(lrecord.0),
                            rpeak: None,
                            lc_start: lrecord.1,
                            lc_end: lrecord.2,
                            rc_start: -1,
                            rc_end: -1,
                            b_no_lclip: false,
                            b_no_rclip: true,
                            b_consist: true,
                            b_polya: false,
                        },
                    );
                }
            }
        }
        for (ins_chrm, m) in m_right.iter() {
            for (&ins_pos, rrecord) in m.iter() {
                let cnt_clip = rrecord.4;
                if (cnt_clip as f64) < n_cutoff / 2.0 + 1.0 {
                    continue;
                }
                if lowq_skip(ins_chrm, ins_pos, cnt_clip) {
                    continue;
                }
                if site_has(&m_left, ins_chrm, ins_pos) {
                    continue;
                }
                m_checked.entry(ins_chrm.clone()).or_default().insert(
                    ins_pos,
                    ClipChk {
                        lpeak: None,
                        rpeak: Some(rrecord.0),
                        lc_start: -1,
                        lc_end: -1,
                        rc_start: rrecord.1,
                        rc_end: rrecord.2,
                        b_no_lclip: true,
                        b_no_rclip: false,
                        b_consist: true,
                        b_polya: false,
                    },
                );
            }
        }
        m_checked
    }

    /// `check_clip_alone_consistency_for_somatic` -> (checked list, sites forming a clip cluster)
    #[allow(clippy::too_many_arguments)]
    pub fn check_clip_alone_consistency_for_somatic(
        m_left_clip_pos: &ClipPosMap,
        m_right_clip_pos: &ClipPosMap,
        sf_ori_list: &str,
        idist: i64,
        ratio: f64,
        n_cutoff: f64,
        bsize: i64,
    ) -> Result<(SiteMap<ClipChk>, SiteMap<ClipChk>)> {
        let empty_rcd = ClipChk {
            lpeak: None,
            rpeak: None,
            lc_start: -1,
            lc_end: -1,
            rc_start: -1,
            rc_end: -1,
            b_no_lclip: false,
            b_no_rclip: false,
            b_consist: false,
            b_polya: false,
        };
        let mut m_checked: SiteMap<ClipChk> = IndexMap::new();
        let mut m_form: SiteMap<ClipChk> = IndexMap::new();
        let m_left = Self::check_one_side_clip_consistency(m_left_clip_pos, idist, ratio, bsize);
        let m_right = Self::check_one_side_clip_consistency(m_right_clip_pos, idist, ratio, bsize);
        for (ins_chrm, m) in m_left.iter() {
            for (&ins_pos, lrecord) in m.iter() {
                let mut cnt_clip = lrecord.4;
                if let Some(rrecord) = site_get(&m_right, ins_chrm, ins_pos) {
                    cnt_clip += rrecord.4;
                    if (cnt_clip as f64) < n_cutoff {
                        m_checked.entry(ins_chrm.clone()).or_default().insert(ins_pos, empty_rcd.clone());
                        continue;
                    }
                    let r = ClipChk {
                        lpeak: Some(lrecord.0),
                        rpeak: Some(rrecord.0),
                        lc_start: lrecord.1,
                        lc_end: lrecord.2,
                        rc_start: rrecord.1,
                        rc_end: rrecord.2,
                        b_no_lclip: false,
                        b_no_rclip: false,
                        b_consist: true,
                        b_polya: false,
                    };
                    m_checked.entry(ins_chrm.clone()).or_default().insert(ins_pos, r.clone());
                    m_form.entry(ins_chrm.clone()).or_default().insert(ins_pos, r);
                } else {
                    if (cnt_clip as f64) < n_cutoff / 2.0 + 1.0 {
                        m_checked.entry(ins_chrm.clone()).or_default().insert(ins_pos, empty_rcd.clone());
                        continue;
                    }
                    let r = ClipChk {
                        lpeak: Some(lrecord.0),
                        rpeak: None,
                        lc_start: lrecord.1,
                        lc_end: lrecord.2,
                        rc_start: -1,
                        rc_end: -1,
                        b_no_lclip: false,
                        b_no_rclip: true,
                        b_consist: true,
                        b_polya: false,
                    };
                    m_checked.entry(ins_chrm.clone()).or_default().insert(ins_pos, r.clone());
                    m_form.entry(ins_chrm.clone()).or_default().insert(ins_pos, r);
                }
            }
        }
        for (ins_chrm, m) in m_right.iter() {
            for (&ins_pos, rrecord) in m.iter() {
                let cnt_clip = rrecord.4;
                if site_has(&m_left, ins_chrm, ins_pos) {
                    continue;
                }
                if (cnt_clip as f64) < n_cutoff / 2.0 + 1.0 {
                    m_checked.entry(ins_chrm.clone()).or_default().insert(ins_pos, empty_rcd.clone());
                    continue;
                }
                let r = ClipChk {
                    lpeak: None,
                    rpeak: Some(rrecord.0),
                    lc_start: -1,
                    lc_end: -1,
                    rc_start: rrecord.1,
                    rc_end: rrecord.2,
                    b_no_lclip: true,
                    b_no_rclip: false,
                    b_consist: true,
                    b_polya: false,
                };
                m_checked.entry(ins_chrm.clone()).or_default().insert(ins_pos, r.clone());
                m_form.entry(ins_chrm.clone()).or_default().insert(ins_pos, r);
            }
        }
        for line in utils::read_lines(sf_ori_list)? {
            let f = split_ws(&line);
            if f.len() < 2 {
                continue;
            }
            let ins_pos = utils::pint(f[1])?;
            let e = m_checked.entry(f[0].to_string()).or_default();
            if !e.contains_key(&ins_pos) {
                e.insert(ins_pos, empty_rcd.clone());
            }
        }
        Ok((m_checked, m_form))
    }

    /// `check_disc_sample_consistency`
    pub fn check_disc_sample_consistency(&self, m_clip_checked: &SiteMap<ClipChk>, m_sample_disc: &SiteMap<IndexMap<String, i64>>, n_disc_cutoff: f64) -> Result<SiteMap<ClipChk>> {
        let mut m_selected: SiteMap<ClipChk> = IndexMap::new();
        let n_bams_all = self.get_n_bams()?;
        for (ins_chrm, m) in m_clip_checked.iter() {
            for (&ins_pos, rcd) in m.iter() {
                let mut b_sample_consist = true;
                let mut tmp_n_bams = 0i64;
                if let Some(ms) = site_get(m_sample_disc, ins_chrm, ins_pos) {
                    for (_, &n) in ms.iter() {
                        tmp_n_bams += 1;
                        if (n as f64) < n_disc_cutoff {
                            b_sample_consist = false;
                            break;
                        }
                    }
                }
                if !b_sample_consist {
                    continue;
                }
                if tmp_n_bams < n_bams_all {
                    continue;
                }
                let e = m_selected.entry(ins_chrm.clone()).or_default();
                if !e.contains_key(&ins_pos) {
                    e.insert(ins_pos, rcd.clone());
                }
            }
        }
        Ok(m_selected)
    }

    /// Split the disc reads of one site into left/right groups (shared by the two versions of
    /// `check_seprt_disc_consistency`).
    fn split_disc(l_disc: &[DiscHit], clip_pos: i64, somatic: bool) -> (Vec<f64>, Vec<f64>, [i64; 4], [i64; 4]) {
        let mut l_l = Vec::new();
        let mut l_r = Vec::new();
        let (mut n_l_rc, mut n_l_not_rc, mut n_r_rc, mut n_r_not_rc) = (0i64, 0i64, 0i64, 0i64);
        let (mut n_l_rc2, mut n_l_non_rc, mut n_r_rc2, mut n_r_non_rc) = (0i64, 0i64, 0i64, 0i64);
        for d in l_disc {
            let mut b_lclip_disc = false;
            let mut b_rclip_disc = false;
            if (d.anchor_pos - clip_pos).abs() < gv::NEARBY_CLIP {
                if d.b_lclip && !d.b_rclip {
                    b_lclip_disc = true;
                }
                if d.b_rclip && !d.b_lclip {
                    b_rclip_disc = true;
                }
            }
            let to_right = if b_lclip_disc {
                true
            } else if b_rclip_disc {
                false
            } else {
                d.anchor_pos >= clip_pos
            };
            if to_right {
                l_r.push(d.cns_pos);
                if d.b_same_ori {
                    n_r_rc += 1;
                } else {
                    n_r_not_rc += 1;
                }
                if d.b_anchor_rc {
                    if somatic {
                        n_r_rc += 1;
                    } else {
                        n_r_rc2 += 1;
                    }
                } else {
                    n_r_non_rc += 1;
                }
            } else {
                l_l.push(d.cns_pos);
                if d.b_same_ori {
                    n_l_rc += 1;
                } else {
                    n_l_not_rc += 1;
                }
                if d.b_anchor_rc {
                    if somatic {
                        n_l_rc += 1;
                    } else {
                        n_l_rc2 += 1;
                    }
                } else {
                    n_l_non_rc += 1;
                }
            }
        }
        if somatic {
            // the somatic version re-uses n_l_rc/n_r_rc for both records
            (l_l, l_r, [n_l_rc, n_l_not_rc, n_r_rc, n_r_not_rc], [n_l_rc, n_l_non_rc, n_r_rc, n_r_non_rc])
        } else {
            (l_l, l_r, [n_l_rc, n_l_not_rc, n_r_rc, n_r_not_rc], [n_l_rc2, n_l_non_rc, n_r_rc2, n_r_non_rc])
        }
    }

    /// `check_seprt_disc_consistency`
    pub fn check_seprt_disc_consistency(
        m_repsnt_pos: &SiteMap<i64>,
        m_disc: &SiteMap<Vec<DiscHit>>,
        icorcord: i64,
        ratio: f64,
        ndisc_cutoff: f64,
        b_transduct: bool,
    ) -> SiteMap<DiscFlt> {
        let mut m_selected: SiteMap<DiscFlt> = IndexMap::new();
        for (ins_chrm, m) in m_repsnt_pos.iter() {
            for (&ins_pos, &clip_pos) in m.iter() {
                let l_disc = match site_get(m_disc, ins_chrm, ins_pos) {
                    Some(l) => l,
                    None => continue,
                };
                let (mut l_l, mut l_r, rc_rcd, anchor_rc_rcd) = Self::split_disc(l_disc, clip_pos, false);
                let n_l_disc = l_l.len() as i64;
                let n_r_disc = l_r.len() as i64;
                let (b_lcluster, lcns_start, lcns_end) = is_disc_cluster_f(&mut l_l, icorcord as f64, ratio);
                let (b_rcluster, rcns_start, rcns_end) = is_disc_cluster_f(&mut l_r, icorcord as f64, ratio);
                if !b_lcluster && !b_rcluster {
                    continue;
                }
                if b_lcluster && b_rcluster && ((n_l_disc + n_r_disc) as f64) < ndisc_cutoff {
                    continue;
                }
                if b_transduct {
                    if b_lcluster && !b_rcluster && (n_l_disc as f64) < ndisc_cutoff / 2.0 {
                        continue;
                    } else if !b_lcluster && b_rcluster && (n_r_disc as f64) < ndisc_cutoff / 2.0 {
                        continue;
                    }
                }
                m_selected
                    .entry(ins_chrm.clone())
                    .or_default()
                    .insert(ins_pos, DiscFlt { n_l_disc, lcns_start, lcns_end, n_r_disc, rcns_start, rcns_end, rc_rcd, anchor_rc_rcd });
            }
        }
        m_selected
    }

    /// `check_seprt_disc_consistency_for_somatic`
    pub fn check_seprt_disc_consistency_for_somatic(m_repsnt_pos: &SiteMap<i64>, m_disc: &SiteMap<Vec<DiscHit>>, icorcord: i64, ratio: f64) -> SiteMap<DiscFlt> {
        let mut m_selected: SiteMap<DiscFlt> = IndexMap::new();
        for (ins_chrm, m) in m_repsnt_pos.iter() {
            for (&ins_pos, &clip_pos) in m.iter() {
                let l_disc = match site_get(m_disc, ins_chrm, ins_pos) {
                    Some(l) => l,
                    None => {
                        m_selected.entry(ins_chrm.clone()).or_default().insert(
                            ins_pos,
                            DiscFlt {
                                n_l_disc: 0,
                                lcns_start: PyNum::Int(-1),
                                lcns_end: PyNum::Int(-1),
                                n_r_disc: 0,
                                rcns_start: PyNum::Int(-1),
                                rcns_end: PyNum::Int(-1),
                                rc_rcd: [0; 4],
                                anchor_rc_rcd: [0; 4],
                            },
                        );
                        continue;
                    }
                };
                let (mut l_l, mut l_r, rc_rcd, anchor_rc_rcd) = Self::split_disc(l_disc, clip_pos, true);
                let n_l_disc = l_l.len() as i64;
                let n_r_disc = l_r.len() as i64;
                let (_, lcns_start, lcns_end) = is_disc_cluster_f(&mut l_l, icorcord as f64, ratio);
                let (_, rcns_start, rcns_end) = is_disc_cluster_f(&mut l_r, icorcord as f64, ratio);
                m_selected
                    .entry(ins_chrm.clone())
                    .or_default()
                    .insert(ins_pos, DiscFlt { n_l_disc, lcns_start, lcns_end, n_r_disc, rcns_start, rcns_end, rc_rcd, anchor_rc_rcd });
            }
        }
        m_selected
    }

    /// `parse_clip_realignment_consensus`
    pub fn parse_clip_realignment_consensus(&self, sf_clip_alignmt: &str, bmapped_cutoff: f64) -> Result<ClipRealign> {
        let mut r = ClipRealign {
            m_ins_lpos: IndexMap::new(),
            m_ins_rpos: IndexMap::new(),
            m_cns_lpos: IndexMap::new(),
            m_cns_rpos: IndexMap::new(),
            m_clip_sample: IndexMap::new(),
            m_polya: IndexMap::new(),
        };
        let mut samfile = AlignmentFile::open(sf_clip_alignmt, &self.sf_reference)?;
        let mut err: Option<anyhow::Error> = None;
        samfile.fetch_all(|a: &Record| {
            if err.is_some() {
                return;
            }
            let read_info = a.query_name();
            let f: Vec<&str> = read_info.split(gv::SEPERATOR).collect();
            if f.len() < 7 {
                err = Some(anyhow!("IndexError: unexpected read name {}", read_info));
                return;
            }
            let ori_chrm = f[0].to_string();
            let ref_mpos = match utils::pint(f[1]) {
                Ok(v) => v,
                Err(e) => {
                    err = Some(e);
                    return;
                }
            };
            let ori_clip_flag = f[2];
            let ori_insertion_pos = match utils::pint(f[4]) {
                Ok(v) => v,
                Err(e) => {
                    err = Some(e);
                    return;
                }
            };
            let sample_id = f[6].to_string();
            if (ref_mpos - ori_insertion_pos).abs() > gv::NEARBY_CLIP {
                return;
            }
            let b_left = ori_clip_flag != gv::FLAG_RIGHT_CLIP;
            if a.is_unmap() {
                return;
            }
            let l_cigar = a.cigar_tuples();
            let (b_ok, n_map_bases) = is_clipped_part_qualified_algnmt(&l_cigar, bmapped_cutoff);
            if !b_ok {
                return;
            }
            let clipped_seq = a.query_sequence().unwrap_or_default();
            let b_clip_part_rc = a.is_rev();
            let n = gv::CK_POLYA_SEQ_MAX;
            let s_clip_seq_ck = if !b_clip_part_rc {
                if !b_left {
                    py_head(&clipped_seq, n)
                } else {
                    py_tail(&clipped_seq, n)
                }
            } else if b_left {
                py_head(&clipped_seq, n)
            } else {
                py_tail(&clipped_seq, n)
            };
            let b_polya = polya::is_consecutive_polya(s_clip_seq_ck);
            if b_polya {
                let p = r.m_polya.entry(ori_chrm.clone()).or_default().entry(ori_insertion_pos).or_insert([0; 4]);
                if b_left {
                    p[0] += 1;
                    if !b_clip_part_rc {
                        p[2] += 1;
                    }
                } else {
                    p[1] += 1;
                    if !b_clip_part_rc {
                        p[3] += 1;
                    }
                }
            }
            let mut pos_in_consensus = a.reference_start();
            if (b_clip_part_rc && b_left) || (!b_clip_part_rc && !b_left) {
                pos_in_consensus += n_map_bases;
            }
            let (m_ins, m_cns) = if b_left { (&mut r.m_ins_lpos, &mut r.m_cns_lpos) } else { (&mut r.m_ins_rpos, &mut r.m_cns_rpos) };
            m_ins
                .entry(ori_chrm.clone())
                .or_default()
                .entry(ori_insertion_pos)
                .or_default()
                .entry(ref_mpos)
                .or_default()
                .push(pos_in_consensus);
            *m_cns.entry(ori_chrm.clone()).or_default().entry(ori_insertion_pos).or_default().entry(pos_in_consensus).or_insert(0) += 1;
            let s = r.m_clip_sample.entry(ori_chrm).or_default().entry(ori_insertion_pos).or_default().entry(sample_id).or_insert([0; 2]);
            if b_left {
                s[0] += 1;
            } else {
                s[1] += 1;
            }
        })?;
        if let Some(e) = err {
            return Err(e);
        }
        Ok(r)
    }

    /// `parse_disc_algnmt_consensus`
    pub fn parse_disc_algnmt_consensus(&self, sf_disc_alignmt: &str, bmapped_cutoff: f64) -> Result<DiscRealign> {
        let mut r = DiscRealign { m_disc_pos: IndexMap::new(), m_disc_sample: IndexMap::new(), m_polya: IndexMap::new() };
        let mut samfile = AlignmentFile::open(sf_disc_alignmt, &self.sf_reference)?;
        let mut err: Option<anyhow::Error> = None;
        samfile.fetch_all(|a: &Record| {
            if err.is_some() {
                return;
            }
            if a.is_unmap() {
                return;
            }
            let l_cigar = a.cigar_tuples();
            let (b_ok, _) = is_clipped_part_qualified_algnmt(&l_cigar, bmapped_cutoff);
            if !b_ok {
                return;
            }
            let read_info = a.query_name();
            let f: Vec<&str> = read_info.split(gv::SEPERATOR).collect();
            let n = f.len();
            if n < 8 {
                err = Some(anyhow!("IndexError: unexpected read name {}", read_info));
                return;
            }
            let s_anchor_lclip = f[n - 8];
            let s_anchor_rclip = f[n - 7];
            let s_anchor_rc = f[n - 6];
            let s_anchor_mate_rc = f[n - 5];
            let anchor_map_pos = match utils::pint(f[n - 4]) {
                Ok(v) => v,
                Err(e) => {
                    err = Some(e);
                    return;
                }
            };
            let ins_chrm = f[n - 3].to_string();
            let ins_pos = match utils::pint(f[n - 2]) {
                Ok(v) => v,
                Err(e) => {
                    err = Some(e);
                    return;
                }
            };
            let sample_id = f[n - 1].to_string();
            let read_seq = a.query_sequence().unwrap_or_default();
            if polya::is_consecutive_polya_t(&read_seq) {
                *r.m_polya.entry(ins_chrm.clone()).or_default().entry(ins_pos).or_insert(0) += 1;
            }
            let pos_in_consensus = a.reference_start() as f64 + (read_seq.len() as f64) / 2.0;
            let b_anchor_mate_rc = s_anchor_mate_rc == "1";
            let mut b_masked_rc = a.is_rev();
            if b_anchor_mate_rc {
                b_masked_rc = !a.is_rev();
            }
            let b_anchor_rc = s_anchor_rc == "1";
            let b_anchor_lclip = s_anchor_lclip == "1";
            let b_anchor_rclip = s_anchor_rclip == "1";
            let b_same_dir_xor = b_masked_rc != b_anchor_rc;
            r.m_disc_pos.entry(ins_chrm.clone()).or_default().entry(ins_pos).or_default().push(DiscHit {
                anchor_pos: anchor_map_pos,
                b_lclip: b_anchor_lclip,
                b_rclip: b_anchor_rclip,
                cns_pos: pos_in_consensus,
                b_same_ori: b_same_dir_xor,
                b_anchor_rc,
            });
            *r.m_disc_sample.entry(ins_chrm).or_default().entry(ins_pos).or_default().entry(sample_id).or_insert(0) += 1;
        })?;
        if let Some(e) = err {
            return Err(e);
        }
        Ok(r)
    }

    /// Append `src` to `out`, adding "~sample_id" to every `every`-th line (starting with the first).
    fn append_with_sample_id(out: &mut dyn Write, src: &str, every: usize, sample_cnt: usize) -> Result<()> {
        for (n_cnt, line) in utils::read_raw_lines(src)?.into_iter().enumerate() {
            if n_cnt % every == 0 {
                let l = format!("{}{}{}\n", utils::rstrip(&line), gv::SEPERATOR, sample_cnt);
                out.write_all(l.as_bytes())?;
            } else {
                out.write_all(line.as_bytes())?;
            }
        }
        Ok(())
    }

    /// `collect_clipped_disc_reads` -> m_lowq_clip
    pub fn collect_clipped_disc_reads(&self, sf_candidate_list: &str, extnd: i64, bin_size: i64, sf_clip_fq: &str, sf_disc_fa: &str) -> Result<SiteMap<i64>> {
        let mut m_lowq_clip: SiteMap<i64> = IndexMap::new();
        let mut fout_clip_fq = utils::create(sf_clip_fq)?;
        let mut fout_disc_fa = utils::create(sf_disc_fa)?;
        for (sample_cnt, sf_bam) in self.bam_list()?.iter().enumerate() {
            let xclip_disc = XClipDisc::new(sf_bam, &self.working_folder, self.n_jobs, &self.sf_reference);
            let sf_disc_fa_tmp = format!("{}temp_disc.fa{}", self.working_folder, sample_cnt);
            let sf_clip_fa_tmp = format!("{}temp_clip.fq{}", self.working_folder, sample_cnt);
            let l_low = xclip_disc.collect_clipped_disc_reads_of_given_list(sf_candidate_list, extnd, bin_size, &sf_clip_fa_tmp, &sf_disc_fa_tmp)?;
            for (c, p, n) in l_low {
                *m_lowq_clip.entry(c).or_default().entry(p).or_insert(0) += n;
            }
            Self::append_with_sample_id(&mut fout_disc_fa, &sf_disc_fa_tmp, 2, sample_cnt)?;
            Self::append_with_sample_id(&mut fout_clip_fq, &sf_clip_fa_tmp, 4, sample_cnt)?;
            Self::clean_file_by_path(&sf_disc_fa_tmp);
            Self::clean_file_by_path(&sf_clip_fa_tmp);
        }
        fout_clip_fq.flush()?;
        fout_disc_fa.flush()?;
        Ok(m_lowq_clip)
    }

    /// `collect_clipped_reads`
    pub fn collect_clipped_reads(&self, sf_candidate_list: &str, extnd: i64, sf_clip_fq: &str) -> Result<()> {
        let mut fout_clip_fq = utils::create(sf_clip_fq)?;
        for (sample_cnt, sf_bam) in self.bam_list()?.iter().enumerate() {
            let xclip_disc = XClipDisc::new(sf_bam, &self.working_folder, self.n_jobs, &self.sf_reference);
            let sf_clip_fa_tmp = format!("{}temp_clip.fq{}", self.working_folder, sample_cnt);
            xclip_disc.collect_clipped_reads_of_given_list(sf_candidate_list, extnd, &sf_clip_fa_tmp)?;
            Self::append_with_sample_id(&mut fout_clip_fq, &sf_clip_fa_tmp, 4, sample_cnt)?;
            Self::clean_file_by_path(&sf_clip_fa_tmp);
        }
        fout_clip_fq.flush()?;
        Ok(())
    }

    /// `collect_disc_reads`
    pub fn collect_disc_reads(&self, sf_candidate_list: &str, extnd: i64, bin_size: i64, sf_disc_fa: &str) -> Result<()> {
        let mut fout_disc_fa = utils::create(sf_disc_fa)?;
        for (sample_cnt, sf_bam) in self.bam_list()?.iter().enumerate() {
            let xclip_disc = XClipDisc::new(sf_bam, &self.working_folder, self.n_jobs, &self.sf_reference);
            let sf_disc_fa_tmp = format!("{}temp_disc.fa{}", self.working_folder, sample_cnt);
            xclip_disc.collect_disc_reads_of_given_list(sf_candidate_list, extnd, bin_size, &sf_disc_fa_tmp)?;
            Self::append_with_sample_id(&mut fout_disc_fa, &sf_disc_fa_tmp, 2, sample_cnt)?;
            Self::clean_file_by_path(&sf_disc_fa_tmp);
        }
        fout_disc_fa.flush()?;
        Ok(())
    }

    /// `call_MEIs_consensus`: the "-N" step.
    #[allow(clippy::too_many_arguments)]
    pub fn call_meis_consensus(
        &self,
        sf_candidate_list: &str,
        _sf_raw_disc: &str,
        extnd: i64,
        bin_size: i64,
        sf_rep_cns: &str,
        _sf_flank: &str,
        _i_flank_lenth: i64,
        bmapped_cutoff: f64,
        i_concord_dist: i64,
        f_concord_ratio: f64,
        nclip_cutoff: i64,
        ndisc_cutoff: i64,
        sf_final_list: &str,
    ) -> Result<()> {
        let sf_cns_log = format!("{}filtering_log.txt", self.working_folder);
        let mut xlog = XLog::open_file(&sf_cns_log)?;
        let sf_clip_fq = format!("{}candidate_sites_all_clip.fq", self.working_folder);
        let sf_disc_fa = format!("{}candidate_sites_all_disc.fa", self.working_folder);
        xlog.append("[Filtering:collect_clipped_disc_reads:Starts...]\n");
        let m_lowq_clip_cnt = self.collect_clipped_disc_reads(sf_candidate_list, extnd, bin_size, &sf_clip_fq, &sf_disc_fa)?;
        xlog.append("[Filtering:collect_clipped_disc_reads:Finished]\n");

        let bwa_align = BwAlign::new(gv::BWA_PATH, gv::BWA_REALIGN_CUTOFF, self.n_jobs as i64);
        let sf_clip_algnmt = format!("{}temp_clip.sam", self.working_folder);
        xlog.append("[Filtering:realign_clipped_read_with_polyA:Starts...]\n");
        bwa_align.realign_clipped_read_with_polya(sf_rep_cns, &sf_clip_fq, &sf_clip_algnmt)?;
        Self::clean_file_by_path(&sf_clip_fq);
        xlog.append("[Filtering:realign_clipped_read_with_polyA:Finished]\n");
        let sf_disc_algnmt = format!("{}temp_disc.sam", self.working_folder);
        xlog.append("[Filtering:realign_disc_reads:Starts...]\n");
        bwa_align.realign_disc_reads(sf_rep_cns, &sf_disc_fa, &sf_disc_algnmt);
        Self::clean_file_by_path(&sf_disc_fa);
        xlog.append("[Filtering:realign_disc_reads:Finished...]\n");

        xlog.append("[Filtering:Clip-Consistency-Checking:Starts...]\n");
        let clip = self.parse_clip_realignment_consensus(&sf_clip_algnmt, bmapped_cutoff)?;
        let idist = gv::CLIP_SEARCH_WINDOW;
        let ratio = gv::CLIP_CONSISTENT_RATIO;
        let exd_bin_size = idist;
        let _n_bam = Self::cnt_n_bams(&self.sf_bam_list)?;
        let m_clip_checked_list = Self::check_clip_alone_consistency2(
            &clip.m_ins_lpos,
            &clip.m_ins_rpos,
            idist,
            ratio,
            nclip_cutoff as f64,
            exd_bin_size,
            Some(&m_lowq_clip_cnt),
        );
        xlog.append("[Filtering:Clip-Consistency-Checking:Finished...]\n");

        xlog.append("[Filtering:Disc-Consistency-Checking:Starts...]\n");
        let disc = self.parse_disc_algnmt_consensus(&sf_disc_algnmt, bmapped_cutoff)?;
        let n_disc_half_cutoff = ndisc_cutoff as f64 / 2.0;
        let nclip_half_cutoff = nclip_cutoff as f64 / 2.0;
        let m_disc_consist_list = if gv::check_by_sample() {
            self.check_disc_sample_consistency(&m_clip_checked_list, &disc.m_disc_sample, n_disc_half_cutoff)?
        } else {
            m_clip_checked_list.clone()
        };
        xlog.append(&format!("[Filtering:Disc-Consistency-Checking: m_disc_consist_list size is {}]\n", m_disc_consist_list.len()));
        let m_represent_clip_pos = Self::refine_represent_clip_pos(&m_disc_consist_list);
        let m_disc_filtered = Self::check_seprt_disc_consistency(
            &m_represent_clip_pos,
            &disc.m_disc_pos,
            i_concord_dist,
            f_concord_ratio,
            ndisc_cutoff as f64,
            false,
        );
        xlog.append("[Filtering:Disc-Consistency-Checking:Finished...]\n");
        let m_transdct_info: SiteMap<TransductInfo> = IndexMap::new();
        let m_rescued: SiteMap<(i64, i64, i64)> = IndexMap::new();
        let mut l_candidates = Self::extract_detailed_info_tei(
            &m_disc_filtered,
            &m_disc_consist_list,
            &clip.m_ins_lpos,
            &clip.m_ins_rpos,
            &clip.m_polya,
            &m_clip_checked_list,
            &m_transdct_info,
            &m_rescued,
            None,
        )?;

        xlog.append("[Filtering:calc_coverage_samples:Starts...]\n");
        let rd = ReadDepth::new(&self.working_folder, self.n_jobs, &self.sf_reference);
        let i_total_cov = gv::average_cov();
        xlog.append(&format!(
            "[Filtering:calc_coverage_samples:Finished with estimated coverage:{}]\n",
            gv::average_cov_str()
        ));
        let i_max_cov = gv::MAX_COV_TIMES as f64 * i_total_cov;
        let sf_before_filtering = format!("{}.before_calling_transduction", sf_final_list);
        Self::dump_tei_info(&l_candidates, &sf_before_filtering, false)?;

        let search_win = gv::COV_SEARCH_WINDOW;
        let focal_win = gv::LOCAL_COV_WIN;
        let focal_win2 = gv::COV_ISD_CHK_WIN;
        let l_sites: Vec<(String, i64)> =
            m_disc_filtered.iter().flat_map(|(c, m)| m.keys().map(move |p| (c.clone(), *p))).collect();
        let m_read_depth = rd.calc_coverage_of_two_regions(&l_sites, &self.sf_bam_list, search_win, focal_win, focal_win2)?;
        let sf_cov = format!("{}.sites_cov", sf_before_filtering);
        ReadDepth::dump_coverage_info(&m_read_depth, &sf_cov)?;
        let mut m_hcov: SiteMap<[i64; 2]> = IndexMap::new();
        for (ins_chrm, m) in m_read_depth.iter() {
            for (&ins_pos, cov) in m.iter() {
                let i_cur_lcov = cov[2].int();
                let i_cur_rcov = cov[3].int();
                if (i_cur_lcov as f64) > i_max_cov || (i_cur_rcov as f64) > i_max_cov {
                    let e = m_hcov.entry(ins_chrm.clone()).or_default();
                    if !e.contains_key(&ins_pos) {
                        e.insert(ins_pos, [i_cur_lcov, i_cur_rcov]);
                    }
                }
            }
        }
        let i_cns_lth = Self::get_min_cns_length(sf_rep_cns)?;
        let m_ins_categories = Self::first_stage_classify(&l_candidates, &m_hcov);

        xlog.append("[Filtering:Filter-Dump-Insertion-Info:Starts...]\n");
        let sf_before_filtering = format!("{}.before_filtering", sf_final_list);
        Self::dump_tei_info_with_ori_pos(&l_candidates, &sf_before_filtering, false)?;

        let x_gntper = XGenotyper::new(&self.sf_reference, &self.working_folder, self.n_jobs);
        let is_extnd = gv::dft_is() as i64;
        xlog.append(&format!("[extension length is: {}]\n", is_extnd));
        let sf_gntp_feature = format!("{}.gntp.features", sf_final_list);
        x_gntper.call_genotype(&self.sf_bam_list, &sf_before_filtering, is_extnd, &sf_gntp_feature)?;
        let m_gntp_info = XGenotyper::load_in_features_from_file_with_cov_cutoff(&sf_gntp_feature, i_max_cov)?;

        let min_l1_end = i_cns_lth - 150;
        let depth_ratio = 0.35;
        let (l_final_candidates, m_high_confidt, m_one_side_low_confident, m_not_hit_end) = Self::filter_by_clip_disc_polya_consistency(
            &mut l_candidates,
            &m_read_depth,
            i_concord_dist,
            nclip_half_cutoff,
            n_disc_half_cutoff,
            min_l1_end,
            depth_ratio,
            &mut xlog,
        )?;
        Self::second_stage_classify_dump(
            &m_ins_categories,
            &l_final_candidates,
            &m_high_confidt,
            &m_one_side_low_confident,
            &m_hcov,
            &m_not_hit_end,
            &m_gntp_info,
            i_cns_lth,
            false,
            sf_final_list,
        )?;
        xlog.append("[Filtering:Filter-Dump-Insertion-Info:Finished...]\n");
        xlog.close();
        Ok(())
    }

    /// `call_clip_disc_cluster` (used by the case-control somatic calling) -> sites forming clip clusters
    #[allow(clippy::too_many_arguments)]
    pub fn call_clip_disc_cluster(
        &self,
        sf_candidate_list: &str,
        extnd: i64,
        bin_size: i64,
        sf_rep_cns: &str,
        bmapped_cutoff: f64,
        i_concord_dist: i64,
        f_concord_ratio: f64,
        nclip_cutoff: f64,
        ndisc_cutoff: f64,
        working_folder: &str,
        sf_out: &str,
    ) -> Result<SiteMap<ClipChk>> {
        let sf_clip_fq = format!("{}candidate_sites_all_clip_from_control.fq", working_folder);
        let sf_disc_fa = format!("{}candidate_sites_all_disc_from_control.fa", working_folder);
        self.collect_clipped_disc_reads(sf_candidate_list, extnd, bin_size, &sf_clip_fq, &sf_disc_fa)?;
        let bwa_align = BwAlign::new(gv::BWA_PATH, gv::BWA_REALIGN_CUTOFF, self.n_jobs as i64);
        let sf_clip_algnmt = format!("{}temp_clip.sam", working_folder);
        bwa_align.realign_clipped_read_with_polya(sf_rep_cns, &sf_clip_fq, &sf_clip_algnmt)?;
        Self::clean_file_by_path(&sf_clip_fq);
        let sf_disc_algnmt = format!("{}temp_disc.sam", working_folder);
        bwa_align.realign_disc_reads(sf_rep_cns, &sf_disc_fa, &sf_disc_algnmt);
        Self::clean_file_by_path(&sf_disc_fa);

        let clip = self.parse_clip_realignment_consensus(&sf_clip_algnmt, bmapped_cutoff)?;
        let idist = gv::CLIP_SEARCH_WINDOW;
        let ratio = gv::CLIP_CONSISTENT_RATIO;
        let (m_clip_checked_list, m_clip_only_cluster) = Self::check_clip_alone_consistency_for_somatic(
            &clip.m_ins_lpos,
            &clip.m_ins_rpos,
            sf_candidate_list,
            idist,
            ratio,
            nclip_cutoff,
            idist,
        )?;
        {
            let mut fout = utils::create(&format!("{}.tmp_clip_chk.txt", sf_out))?;
            for (ins_chrm, m) in m_clip_checked_list.iter() {
                for (ins_pos, r) in m.iter() {
                    writeln!(fout, "{}\t{}\t{}:{}\t{}:{}", ins_chrm, ins_pos, r.lc_start, r.lc_end, r.rc_start, r.rc_end)?;
                }
            }
        }
        let disc = self.parse_disc_algnmt_consensus(&sf_disc_algnmt, bmapped_cutoff)?;
        let n_disc_half_cutoff = ndisc_cutoff / 2.0;
        let m_disc_consist_list = if gv::check_by_sample() {
            self.check_disc_sample_consistency(&m_clip_checked_list, &disc.m_disc_sample, n_disc_half_cutoff)?
        } else {
            m_clip_checked_list.clone()
        };
        let m_represent_clip_pos = Self::refine_represent_clip_pos(&m_disc_consist_list);
        let m_disc_filtered =
            Self::check_seprt_disc_consistency(&m_represent_clip_pos, &disc.m_disc_pos, i_concord_dist, f_concord_ratio, ndisc_cutoff, false);
        let m_transdct_info: SiteMap<TransductInfo> = IndexMap::new();
        let m_rescued: SiteMap<(i64, i64, i64)> = IndexMap::new();
        let l_candidates = Self::extract_detailed_info_tei(
            &m_disc_filtered,
            &m_disc_consist_list,
            &clip.m_ins_lpos,
            &clip.m_ins_rpos,
            &clip.m_polya,
            &m_clip_checked_list,
            &m_transdct_info,
            &m_rescued,
            None,
        )?;
        Self::dump_tei_info_with_ori_pos(&l_candidates, sf_out, false)?;
        Ok(m_clip_only_cluster)
    }
}
