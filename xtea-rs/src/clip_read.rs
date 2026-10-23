//! Port of `clip_read.py` (ClipReadInfo): collection of clipped reads/positions and counting
//! of clipped parts re-aligned to the repeat library.

use crate::annotation::XAnnotation;
use crate::bam::{AlignmentFile, PyAln, Record};
use crate::chrom::is_decoy_contig_chrms;
use crate::globals as gv;
use crate::polya;
use crate::utils::{self, split_ws};
use anyhow::Result;
use indexmap::IndexMap;
use std::io::Write;
use std::sync::Arc;

pub struct ClipReadInfo {
    pub sf_bam: String,
    pub working_folder: String,
    pub n_jobs: usize,
    pub sf_reference: String,
}

/// `_cvt_to_Ascii_quality`
pub fn cvt_to_ascii_quality(q: &[u8]) -> String {
    crate::alignments::cvt_to_ascii_quality(q)
}

/// `_is_qualified_clip`: false when more than half of the bases have quality < 15.
pub fn is_qualified_clip(l_score: &[u8]) -> bool {
    let n_half = l_score.len() as f64 / 2.0;
    let mut n_cnt = 0i64;
    for &s in l_score {
        if (s as i64) < gv::CLIP_PHRED_SCORE_CUTOFF {
            n_cnt += 1;
        }
        if n_cnt as f64 > n_half {
            return false;
        }
    }
    true
}

fn clipped_rname(query_name: &str, mate_chrm: &str, mate_pos: i64, chrm: &str, map_pos: i64, flag: &str, b_first: bool) -> String {
    let s = gv::SEPERATOR;
    format!(
        "{q}{s}{mc}{s}{mp}{s}{c}{s}{p}{s}{f}{s}{n}",
        q = query_name,
        mc = mate_chrm,
        mp = mate_pos,
        c = chrm,
        p = map_pos,
        f = flag,
        n = if b_first { 1 } else { 2 },
        s = s
    )
}

impl ClipReadInfo {
    pub fn new(sf_bam: &str, n_jobs: usize, sf_ref: &str) -> Self {
        ClipReadInfo { sf_bam: sf_bam.to_string(), working_folder: "./".to_string(), n_jobs, sf_reference: sf_ref.to_string() }
    }

    pub fn set_working_folder(&mut self, wf: &str) {
        self.working_folder = utils::with_slash(wf);
    }

    fn pool(&self) -> rayon::ThreadPool {
        rayon::ThreadPoolBuilder::new().num_threads(self.n_jobs.max(1)).build().unwrap()
    }

    fn mate_info(a: &Record, names: &[String]) -> Option<(String, i64)> {
        let ntid = a.next_reference_id();
        if ntid < -1 || ntid >= names.len() as i32 {
            return None; // not in m_chrm_id
        }
        if !a.mate_is_unmapped() && ntid >= 0 {
            Some((names[ntid as usize].clone(), a.next_reference_start()))
        } else {
            Some(("*".to_string(), 0))
        }
    }

    /// `collect_clip_positions_by_chrm`
    pub fn collect_clip_positions_by_chrm(&self, chrm: &str, working_folder: &str, xannotation: &XAnnotation, i_clip_cutoff: i64, b_se: bool) -> Result<()> {
        let mut samfile = AlignmentFile::open(&self.sf_bam, &self.sf_reference)?;
        let names: Vec<String> = samfile.references().to_vec();
        let mut m_clip_pos: IndexMap<i64, [i64; 3]> = IndexMap::new();
        let min_mapq = gv::minimum_clip_mapq();
        samfile.fetch_contig(chrm, |a| {
            if a.is_dup() || a.is_unmap() {
                return;
            }
            let cig = a.cigar_tuples();
            if cig.is_empty() {
                return;
            }
            if cig.len() == 1 && cig[0].0 == 0 {
                return;
            }
            if a.mapping_quality() < min_mapq {
                return;
            }
            let (mate_chrm, mate_pos) = match Self::mate_info(a, &names) {
                Some(x) => x,
                None => return,
            };
            let mut map_pos = a.reference_start();
            let mut b_mate_in_rep = false;
            if !b_se {
                b_mate_in_rep = xannotation.in_repeat(&mate_chrm, mate_pos);
            }
            if cig[0].0 == 4 || cig[0].0 == 5 {
                let e = m_clip_pos.entry(map_pos).or_insert([0, 0, 0]);
                e[0] += 1;
                if b_mate_in_rep {
                    e[2] += 1;
                }
                if a.is_supp() || a.is_sec() {
                    return;
                }
            }
            let last = cig[cig.len() - 1];
            if last.0 == 4 || last.0 == 5 {
                for &(t, l) in &cig[..cig.len() - 1] {
                    if t == 4 || t == 5 || t == 1 {
                        continue;
                    }
                    map_pos += l as i64;
                }
                let e = m_clip_pos.entry(map_pos).or_insert([0, 0, 0]);
                e[1] += 1;
                if b_mate_in_rep {
                    e[2] += 1;
                }
            }
        })?;
        let sf_clip_pos = format!("{}{}{}", working_folder, chrm, gv::CLIP_POS_SUFFIX);
        let mut f = utils::create(&sf_clip_pos)?;
        for (pos, v) in m_clip_pos.iter() {
            if v[0] + v[1] >= i_clip_cutoff {
                writeln!(f, "{}\t{}\t{}\t{}", pos, v[0], v[1], v[2])?;
            }
        }
        f.flush()?;
        Ok(())
    }

    /// `collect_clip_positions`
    pub fn collect_clip_positions(&self, sf_annotation: &str, i_clip_cutoff: i64, b_se: bool, sf_pub_folder: &str) -> Result<()> {
        let bam_info = crate::alignments::BamInfo::new(&self.sf_bam, &self.sf_reference);
        let b_with_chr = bam_info.is_chrm_contain_chr()?;
        let samfile = AlignmentFile::open(&self.sf_bam, &self.sf_reference)?;
        let chrms: Vec<String> = samfile.references().iter().filter(|c| !is_decoy_contig_chrms(c)).cloned().collect();
        drop(samfile);
        // Every worker of the original loads the same annotation; load it once and share.
        let mut xannotation = XAnnotation::new(sf_annotation);
        xannotation.set_with_chr(b_with_chr);
        let i_boundary_extnd = if gv::is_call_sva() { gv::SVA_ANNOTATION_EXTND } else { 0 };
        xannotation.load_rmsk_annotation_with_extnd_with_lenth_cutoff(i_boundary_extnd, 0)?;
        xannotation.index_rmsk_annotation_interval_tree();
        let xann = Arc::new(xannotation);
        let wf = self.working_folder.clone();
        let results: Vec<Result<()>> = self.pool().install(|| {
            utils::par_map(chrms.len(), |i| self.collect_clip_positions_by_chrm(&chrms[i], &wf, &xann, i_clip_cutoff, b_se))
        });
        for r in results {
            r?;
        }
        for chrm in chrms.iter() {
            let sf_clip_pos = format!("{}{}{}", self.working_folder, chrm, gv::CLIP_POS_SUFFIX);
            let sf_pub_pos = format!("{}{}{}", sf_pub_folder, chrm, gv::CLIP_POS_SUFFIX);
            if utils::is_link(&sf_pub_pos) || utils::is_file(&sf_pub_pos) {
                utils::remove_file(&sf_pub_pos);
            }
            utils::run_cmd_quiet(&format!("ln -s {} {}", sf_clip_pos, sf_pub_folder));
        }
        Ok(())
    }

    /// `collect_clipped_parts_by_chrm`
    pub fn collect_clipped_parts_by_chrm(&self, chrm: &str, working_folder: &str) -> Result<()> {
        let sf_clip_pos = format!("{}{}{}", working_folder, chrm, gv::CLIP_POS_SUFFIX);
        if !utils::exists(&sf_clip_pos) {
            return Ok(());
        }
        let mut m_pos = std::collections::HashSet::new();
        for line in utils::read_lines(&sf_clip_pos)? {
            let f = split_ws(&line);
            if f.is_empty() {
                continue;
            }
            m_pos.insert(utils::pint(f[0])?);
        }
        let sf_clip_fq = format!("{}{}{}", working_folder, chrm, gv::CLIP_FQ_SUFFIX);
        let mut fq = utils::create(&sf_clip_fq)?;
        let mut samfile = AlignmentFile::open(&self.sf_bam, &self.sf_reference)?;
        let names: Vec<String> = samfile.references().to_vec();
        let mut werr: Option<std::io::Error> = None;
        samfile.fetch_contig(chrm, |a| {
            if a.is_dup() || a.is_unmap() {
                return;
            }
            let b_first = !a.is_read2();
            let cig = a.cigar_tuples();
            if cig.is_empty() {
                return;
            }
            if cig.len() == 1 && cig[0].0 == 0 {
                return;
            }
            let query_name = a.query_name();
            let query_seq = a.query_sequence().unwrap_or_default();
            let query_quality = a.query_qualities().unwrap_or_default();
            let mut map_pos = a.reference_start();
            let (mate_chrm, mate_pos) = match Self::mate_info(a, &names) {
                Some(x) => x,
                None => return,
            };
            if cig[0].0 == 4 {
                if a.is_supp() || a.is_sec() {
                    return;
                }
                if m_pos.contains(&map_pos) {
                    let n = cig[0].1 as usize;
                    let clipped_seq = &query_seq[..n.min(query_seq.len())];
                    let len_clip = clipped_seq.len() as i64;
                    // NOTE: a `continue` in the original skips the whole read (incl. its right clip)
                    if len_clip >= gv::MINIMUM_POLYA_CLIP && len_clip <= gv::BWA_REALIGN_CUTOFF {
                        if !polya::contain_enough_a_t(clipped_seq, len_clip) {
                            return;
                        }
                    } else if len_clip < gv::BWA_REALIGN_CUTOFF {
                        return;
                    }
                    let lq = &query_quality[..n.min(query_quality.len())];
                    if !is_qualified_clip(lq) {
                        return;
                    }
                    let qs = cvt_to_ascii_quality(lq);
                    let rname = clipped_rname(&query_name, &mate_chrm, mate_pos, chrm, map_pos, gv::FLAG_LEFT_CLIP, b_first);
                    if let Err(e) = write!(fq, "@{}\n{}\n+\n{}\n", rname, clipped_seq, qs) {
                        werr = Some(e);
                    }
                }
            }
            let last = cig[cig.len() - 1];
            if last.0 == 4 {
                for &(t, l) in &cig[..cig.len() - 1] {
                    if t == 4 || t == 5 || t == 1 {
                        continue;
                    }
                    map_pos += l as i64;
                }
                if a.is_supp() || a.is_sec() {
                    return;
                }
                if m_pos.contains(&map_pos) {
                    let n = last.1 as usize;
                    // Python: query_seq[-n:] (n == 0 gives the whole sequence)
                    let start = if n == 0 { 0 } else { query_seq.len().saturating_sub(n) };
                    let clipped_seq = &query_seq[start..];
                    let len_clip = clipped_seq.len() as i64;
                    if len_clip >= gv::MINIMUM_POLYA_CLIP && len_clip <= gv::BWA_REALIGN_CUTOFF {
                        if !polya::contain_enough_a_t(clipped_seq, len_clip) {
                            return;
                        }
                    } else if len_clip < gv::BWA_REALIGN_CUTOFF {
                        return;
                    }
                    let qstart = if n == 0 { 0 } else { query_quality.len().saturating_sub(n) };
                    let lq = &query_quality[qstart..];
                    if !is_qualified_clip(lq) {
                        return;
                    }
                    let qs = cvt_to_ascii_quality(lq);
                    let rname = clipped_rname(&query_name, &mate_chrm, mate_pos, chrm, map_pos, gv::FLAG_RIGHT_CLIP, b_first);
                    if let Err(e) = write!(fq, "@{}\n{}\n+\n{}\n", rname, clipped_seq, qs) {
                        werr = Some(e);
                    }
                }
            }
        })?;
        if let Some(e) = werr {
            return Err(e.into());
        }
        fq.flush()?;
        Ok(())
    }

    /// `collect_clipped_parts`
    pub fn collect_clipped_parts(&self, sf_all_clip_fq: &str) -> Result<()> {
        let samfile = AlignmentFile::open(&self.sf_bam, &self.sf_reference)?;
        let references: Vec<String> = samfile.references().to_vec();
        drop(samfile);
        let mut chrms = Vec::new();
        for c in references.iter() {
            if is_decoy_contig_chrms(c) {
                println!("Skip chromosome {}", c);
                continue;
            }
            chrms.push(c.clone());
        }
        let wf = self.working_folder.clone();
        let results: Vec<Result<()>> = self.pool().install(|| utils::par_map(chrms.len(), |i| self.collect_clipped_parts_by_chrm(&chrms[i], &wf)));
        for r in results {
            r?;
        }
        let mut fout = utils::create(sf_all_clip_fq)?;
        for c in references.iter() {
            let sf = format!("{}{}{}", self.working_folder, c, gv::CLIP_FQ_SUFFIX);
            if !utils::is_file(&sf) {
                continue;
            }
            let data = std::fs::read(&sf)?;
            fout.write_all(&data)?;
            utils::remove_file(&sf);
        }
        fout.flush()?;
        Ok(())
    }

    /// Shared filter of `run_cnt_clip_part_aligned_to_rep_by_chrm*`: returns the parsed
    /// (ori_chrm, ori_mpos, flag) when the re-aligned clipped part qualifies.
    fn qualified_realigned_clip(a: &Record) -> Option<(String, i64, String)> {
        if a.is_sec() || a.is_supp() || a.is_unmap() {
            return None;
        }
        let qname = a.query_name();
        let f: Vec<&str> = qname.split(gv::SEPERATOR).collect();
        if f.len() < 4 {
            return None;
        }
        let ori_chrm = f[f.len() - 4].to_string();
        let ori_mpos: i64 = f[f.len() - 3].parse().ok()?;
        let flag = f[f.len() - 2].to_string();
        let cig = a.cigar_tuples();
        if cig.is_empty() {
            return None;
        }
        if cig.len() > 2 {
            let (mut bl, mut ll, mut br, mut lr) = (false, 0i64, false, 0i64);
            if cig[0].0 == 4 || cig[0].0 == 5 {
                bl = true;
                ll = cig[0].1 as i64;
            }
            let last = cig[cig.len() - 1];
            if last.0 == 4 || last.0 == 5 {
                br = true;
                lr = last.1 as i64;
            }
            if bl && br && ll > gv::MAX_CLIP_CLIP_LEN && lr > gv::MAX_CLIP_CLIP_LEN {
                return None;
            }
        }
        let mut n_total: i64 = 0;
        let mut n_map: i64 = 0;
        for &(t, l) in &cig {
            if t == 0 {
                n_map += l as i64;
            }
            if t != 2 {
                n_total += l as i64;
            }
        }
        if (n_map as f64) < (n_total as f64) * 3.0 / 4.0 {
            return None;
        }
        Some((ori_chrm, ori_mpos, flag))
    }

    /// `cnt_clip_part_aligned_to_rep` (+ per chromosome `run_cnt_clip_part_aligned_to_rep_by_chrm`).
    /// The SAM is read once and the records are dispatched per chromosome (equivalent output).
    pub fn cnt_clip_part_aligned_to_rep(&self, sf_clip_sam: &str) -> Result<()> {
        self.cnt_clip_part_impl(sf_clip_sam, false)
    }

    /// `cnt_clip_part_aligned_to_rep_mosaic`
    pub fn cnt_clip_part_aligned_to_rep_mosaic(&self, sf_clip_sam: &str) -> Result<()> {
        self.cnt_clip_part_impl(sf_clip_sam, true)
    }

    fn cnt_clip_part_impl(&self, sf_clip_sam: &str, mosaic: bool) -> Result<()> {
        let samfile = AlignmentFile::open(&self.sf_bam, &self.sf_reference)?;
        let chrms: Vec<String> = samfile.references().iter().filter(|c| !is_decoy_contig_chrms(c)).cloned().collect();
        drop(samfile);
        let wanted: std::collections::HashSet<&String> = chrms.iter().collect();
        // per chromosome: pos -> [left, right], polyA -> [l, r]
        let mut sites: IndexMap<String, IndexMap<i64, [i64; 2]>> = IndexMap::new();
        let mut polya_cnt: IndexMap<String, IndexMap<i64, [i64; 2]>> = IndexMap::new();
        let mut sam = AlignmentFile::open(sf_clip_sam, "")?;
        sam.fetch_all(|a| {
            let (ori_chrm, ori_mpos, flag) = match Self::qualified_realigned_clip(a) {
                Some(x) => x,
                None => return,
            };
            if !wanted.contains(&ori_chrm) {
                return;
            }
            let b_left = flag != gv::FLAG_RIGHT_CLIP;
            if mosaic {
                let b_clip_rc = a.is_rev();
                let clipped_seq = a.query_sequence().unwrap_or_default();
                let mx = gv::CK_POLYA_SEQ_MAX as usize;
                let take_head = |s: &str| -> String { if s.len() > mx { s[..mx].to_string() } else { s.to_string() } };
                let take_tail = |s: &str| -> String { if s.len() > mx { s[s.len() - mx..].to_string() } else { s.to_string() } };
                let ck = if !b_clip_rc {
                    if !b_left { take_head(&clipped_seq) } else { take_tail(&clipped_seq) }
                } else if b_left {
                    take_head(&clipped_seq)
                } else {
                    take_tail(&clipped_seq)
                };
                let b_polya = polya::is_consecutive_polya_t(&ck);
                let e = polya_cnt.entry(ori_chrm.clone()).or_default().entry(ori_mpos).or_insert([0, 0]);
                if b_polya {
                    if b_left {
                        e[0] += 1;
                    } else {
                        e[1] += 1;
                    }
                }
            }
            let e = sites.entry(ori_chrm).or_default().entry(ori_mpos).or_insert([0, 0]);
            if b_left {
                e[0] += 1;
            } else {
                e[1] += 1;
            }
        })?;
        for chrm in chrms.iter() {
            let sf_out = format!("{}{}{}", self.working_folder, chrm, gv::CLIP_RE_ALIGN_POS_SUFFIX);
            let mut f = utils::create(&sf_out)?;
            if let Some(m) = sites.get(chrm) {
                for (pos, v) in m.iter() {
                    if mosaic {
                        let p = polya_cnt.get(chrm).and_then(|x| x.get(pos)).cloned().unwrap_or([0, 0]);
                        write!(f, "{}\t{}\t{}\t{}\t{}\n", pos, v[0], v[1], p[0], p[1])?;
                    } else {
                        write!(f, "{}\t{}\t{}\t\n", pos, v[0], v[1])?;
                    }
                }
            }
            f.flush()?;
        }
        Ok(())
    }

    fn load_realign_pos(&self, chrm: &str) -> Result<std::collections::HashMap<i64, String>> {
        let mut m = std::collections::HashMap::new();
        let sf = format!("{}{}{}", self.working_folder, chrm, gv::CLIP_RE_ALIGN_POS_SUFFIX);
        if utils::is_file(&sf) {
            for line in utils::read_lines(&sf)? {
                let f = split_ws(&line);
                if f.is_empty() {
                    continue;
                }
                m.insert(utils::pint(f[0])?, f[1..].join("\t"));
            }
        }
        Ok(m)
    }

    /// `merge_clip_positions`
    pub fn merge_clip_positions(&self, sf_pclip: &str, sf_out: &str) -> Result<()> {
        self.merge_impl(sf_pclip, sf_out, None, "0\t0")
    }

    /// `merge_clip_positions_with_cutoff`
    pub fn merge_clip_positions_with_cutoff(&self, cutoff_left_clip: i64, cutoff_right_clip: i64, max_cov_cutoff: i64, sf_pclip: &str, sf_out: &str) -> Result<()> {
        self.merge_impl(sf_pclip, sf_out, Some((cutoff_left_clip, cutoff_right_clip, max_cov_cutoff)), "0\t0")
    }

    /// `merge_clip_positions_with_cutoff_polyA`
    pub fn merge_clip_positions_with_cutoff_polya(&self, cutoff_left_clip: i64, cutoff_right_clip: i64, max_cov_cutoff: i64, sf_pclip: &str, sf_out: &str) -> Result<()> {
        self.merge_impl(sf_pclip, sf_out, Some((cutoff_left_clip, cutoff_right_clip, max_cov_cutoff)), "0\t0\t0\t0")
    }

    fn merge_impl(&self, sf_pclip: &str, sf_out: &str, cut: Option<(i64, i64, i64)>, empty: &str) -> Result<()> {
        let samfile = AlignmentFile::open(&self.sf_bam, &self.sf_reference)?;
        let references: Vec<String> = samfile.references().to_vec();
        drop(samfile);
        let mut fout = utils::create(sf_out)?;
        for chrm in references.iter() {
            let m_realign = self.load_realign_pos(chrm)?;
            let sf_clip_pos = format!("{}{}{}", sf_pclip, chrm, gv::CLIP_POS_SUFFIX);
            if !utils::is_file(&sf_clip_pos) {
                continue;
            }
            for line in utils::read_raw_lines(&sf_clip_pos)? {
                let f = split_ws(&line);
                if f.is_empty() {
                    continue;
                }
                let pos = utils::pint(f[0])?;
                if let Some((cl, cr, mx)) = cut {
                    let il = utils::pint(f[1])?;
                    let ir = utils::pint(f[2])?;
                    if il < cl && ir < cr {
                        continue;
                    }
                    if il + ir > mx {
                        continue;
                    }
                }
                let tail = match m_realign.get(&pos) {
                    Some(s) => s.as_str(),
                    None => empty,
                };
                write!(fout, "{}\t{}\t{}\n", chrm, utils::rstrip(&line), tail)?;
            }
        }
        fout.flush()?;
        Ok(())
    }
}
