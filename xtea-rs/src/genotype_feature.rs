//! Port of `x_genotype_feature.py` (`XGenotyper`): clip/discordant features used for genotyping.

use crate::alignments::BamInfo;
use crate::annotation::process_chrm_name;
use crate::bam::{AlignmentFile, PyAln, Record, BAM_CDEL, BAM_CINS, BAM_CMATCH, BAM_CSOFT_CLIP};
use crate::chrom::is_decoy_contig_chrms;
use crate::clip_disc_filter::{is_discordant, py_head, py_tail, right_clip_pos};
use crate::globals as gv;
use crate::intermediate_sites::SiteMap;
use crate::polya;
use crate::utils::{self, split_ws};
use anyhow::{anyhow, Result};
use indexmap::IndexMap;
use std::collections::HashSet;
use std::io::Write;

/// (n_af_clip, n_full_map, n_l_raw_clip, n_r_raw_clip, n_disc_pairs, n_concd_pairs,
///  n_disc_large_indel, s_clip_lens, n_polyA, n_disc_chrms)
#[derive(Clone, Debug)]
pub struct GntpFeat {
    pub n_af_clip: i64,
    pub n_full_map: i64,
    pub n_l_raw_clip: i64,
    pub n_r_raw_clip: i64,
    pub n_disc_pairs: i64,
    pub n_concd_pairs: i64,
    pub n_disc_large_indel: i64,
    pub s_clip_lens: String,
    pub n_polya: i64,
    pub n_disc_chrms: i64,
}

impl GntpFeat {
    fn fmt_fields(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.n_af_clip,
            self.n_full_map,
            self.n_l_raw_clip,
            self.n_r_raw_clip,
            self.n_disc_pairs,
            self.n_concd_pairs,
            self.n_disc_large_indel,
            self.s_clip_lens,
            self.n_polya,
            self.n_disc_chrms
        )
    }
}

pub struct XGenotyper {
    pub sf_reference: String,
    pub working_folder: String,
    pub n_jobs: usize,
}

/// `_has_large_indel_in_read` (the last operation is not checked)
fn has_large_indel_in_read(l_cigar: &[(u32, u32)]) -> bool {
    for &(t, l) in &l_cigar[..l_cigar.len().saturating_sub(1)] {
        if (t == BAM_CINS || t == BAM_CDEL) && (l as i64) >= gv::LARGE_INDEL_IN_READ {
            return true;
        }
    }
    false
}

/// `is_concrdant`
fn is_concrdant(chrm: &str, map_pos: i64, mate_chrm: &str, mate_pos: i64, ins_pos: i64, i_is: f64) -> bool {
    if chrm == mate_chrm {
        let rl = gv::read_length();
        let (i_start, i_end) = if map_pos > mate_pos {
            (mate_pos as f64, map_pos as f64 + rl)
        } else {
            (map_pos as f64, mate_pos as f64 + rl)
        };
        let ip = ins_pos as f64;
        if i_start <= ip && i_end >= ip && (ip - i_start) <= i_is && (i_end - ip) <= i_is {
            return true;
        }
    }
    false
}

impl XGenotyper {
    pub fn new(sf_ref: &str, working_folder: &str, n_jobs: usize) -> Self {
        XGenotyper { sf_reference: sf_ref.to_string(), working_folder: working_folder.to_string(), n_jobs }
    }

    /// `collect_features_one_site` -> the feature line of the site
    pub fn collect_features_one_site(&self, chrm: &str, ins_pos: i64, extnd: f64, sf_bam: &str) -> Result<String> {
        // the extension may be a float (DFT_IS set to the mean insert size); pysam truncates the
        // region coordinates
        let mut start_pos = ins_pos as f64 - extnd;
        if start_pos <= 0.0 {
            start_pos = 1.0;
        }
        let end_pos = ins_pos as f64 + extnd;
        let (start_pos, end_pos) = (start_pos.trunc() as i64, end_pos.trunc() as i64);
        let bam_info = BamInfo::new(sf_bam, &self.sf_reference);
        let b_with_chr = bam_info.is_chrm_contain_chr()?;
        let chrm_in_bam = process_chrm_name(chrm, b_with_chr);
        let mut samfile = AlignmentFile::open(sf_bam, &self.sf_reference)?;
        let names: Vec<String> = samfile.references().to_vec();
        let n_refs = names.len() as i32;
        let mut n_disc_pairs = 0i64;
        let mut n_concd_pairs = 0i64;
        let mut n_l_raw_clip = 0i64;
        let mut n_r_raw_clip = 0i64;
        let mut n_l_af_clip = 0i64;
        let mut n_r_af_clip = 0i64;
        let mut n_full_map = 0i64;
        let mut n_l_full_map = 0i64;
        let mut n_r_full_map = 0i64;
        let mut l_lclip_lens: Vec<String> = Vec::new();
        let mut l_rclip_lens: Vec<String> = Vec::new();
        let mut n_disc_large_indel = 0i64;
        let mut n_polya = 0i64;
        let mut m_mate_chrms: HashSet<String> = HashSet::new();
        let mut l_check_concord: Vec<(String, i64, String, i64)> = Vec::new();
        let mut m_clip_qname: HashSet<String> = HashSet::new();
        let dft_is = gv::dft_is();
        samfile.fetch(&chrm_in_bam, start_pos, end_pos, |a: &Record| {
            if a.is_dup() {
                return;
            }
            if a.is_unmap() {
                return;
            }
            let l_cigar = a.cigar_tuples();
            if l_cigar.is_empty() {
                return;
            }
            let query_name = a.query_name();
            let query_seq = a.query_sequence().unwrap_or_default();
            let mut map_pos = a.reference_start();
            let mut mate_chrm = "*".to_string();
            let mut mate_pos = 0i64;
            let b_clip_part_rc = a.is_rev();
            let ntid = a.next_reference_id();
            if ntid >= 0 && ntid < n_refs && !a.mate_is_unmapped() {
                mate_chrm = names[ntid as usize].clone();
                mate_pos = a.next_reference_start();
            }
            let first = l_cigar[0];
            let last = l_cigar[l_cigar.len() - 1];
            if l_cigar.len() == 1 && first.0 == BAM_CMATCH {
                let i_map_start = map_pos;
                let i_map_end = i_map_start + first.1 as i64;
                if ins_pos >= i_map_start && ins_pos <= i_map_end {
                    n_full_map += 1;
                    if (i_map_start - ins_pos).abs() < gv::BWA_HALF_READ_MIN_SCORE {
                        n_r_full_map += 1;
                    } else if (i_map_end - ins_pos).abs() < gv::BWA_HALF_READ_MIN_SCORE {
                        n_l_full_map += 1;
                    }
                }
            }
            let mut s_clip_seq_ck: &str = "";
            let n_ck = gv::CK_POLYA_SEQ_MAX;
            if first.0 == BAM_CSOFT_CLIP {
                if a.is_supp() || a.is_sec() {
                    return;
                }
                let clipped_seq = py_head(&query_seq, first.1 as i64);
                let len_clip_seq = clipped_seq.len();
                if (map_pos - ins_pos).abs() < gv::TSD_CUTOFF {
                    n_l_raw_clip += 1;
                    m_clip_qname.insert(query_name.clone());
                }
                if (map_pos - ins_pos).abs() < gv::CK_POLYA_CLIP_WIN {
                    s_clip_seq_ck = if !b_clip_part_rc { py_tail(clipped_seq, n_ck) } else { py_head(clipped_seq, n_ck) };
                }
                if (map_pos - ins_pos).abs() < gv::CLIP_EXACT_CLIP_SLACK {
                    n_l_af_clip += 1;
                    l_lclip_lens.push(len_clip_seq.to_string());
                }
            }
            if last.0 == BAM_CSOFT_CLIP {
                map_pos = right_clip_pos(&l_cigar, map_pos);
                if a.is_supp() || a.is_sec() {
                    return;
                }
                let clipped_seq = py_tail(&query_seq, last.1 as i64);
                let len_clip_seq = clipped_seq.len();
                if (map_pos - ins_pos).abs() < gv::TSD_CUTOFF {
                    n_r_raw_clip += 1;
                    m_clip_qname.insert(query_name.clone());
                }
                if (map_pos - ins_pos).abs() < gv::CK_POLYA_CLIP_WIN {
                    s_clip_seq_ck = if !b_clip_part_rc { py_head(clipped_seq, n_ck) } else { py_tail(clipped_seq, n_ck) };
                }
                if (map_pos - ins_pos).abs() < gv::CLIP_EXACT_CLIP_SLACK {
                    n_r_af_clip += 1;
                    l_rclip_lens.push(len_clip_seq.to_string());
                }
            }
            if polya::is_consecutive_polya_t(s_clip_seq_ck) {
                n_polya += 1;
            }
            if mate_chrm == "*" {
                return;
            }
            if is_decoy_contig_chrms(&mate_chrm) {
                return;
            }
            m_mate_chrms.insert(mate_chrm.clone());
            if is_discordant(&chrm_in_bam, map_pos, &mate_chrm, mate_pos, gv::DISC_THRESHOLD) {
                if has_large_indel_in_read(&l_cigar) {
                    n_disc_large_indel += 1;
                }
                if ((map_pos - ins_pos).abs() as f64) <= dft_is {
                    n_disc_pairs += 1;
                }
            } else {
                l_check_concord.push((query_name, map_pos, mate_chrm, mate_pos));
            }
        })?;
        for (qname, map_pos, mate_chrm, mate_pos) in l_check_concord.iter() {
            if m_clip_qname.contains(qname) {
                continue;
            }
            if is_concrdant(&chrm_in_bam, *map_pos, mate_chrm, *mate_pos, ins_pos, dft_is) {
                n_concd_pairs += 1;
            }
        }
        let mut n_af_clip = n_l_af_clip;
        let mut n_full = n_full_map - n_l_full_map;
        let mut s_clip_lens = l_lclip_lens.join(":");
        if n_l_af_clip < n_r_af_clip {
            n_af_clip = n_r_af_clip;
            n_full -= n_r_full_map;
            s_clip_lens = l_rclip_lens.join(":");
        }
        if s_clip_lens.is_empty() {
            s_clip_lens = "none".to_string();
        }
        let n_disc_chrms = m_mate_chrms.len();
        Ok(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            chrm, ins_pos, n_af_clip, n_full, n_l_raw_clip, n_r_raw_clip, n_disc_pairs, n_concd_pairs, n_disc_large_indel, s_clip_lens, n_polya, n_disc_chrms
        ))
    }

    /// `call_genotype`: collect the features of the sites of `sf_candidate_list` for each BAM and
    /// merge them (per-BAM outputs are kept as `sf_out{i}.out`).
    pub fn call_genotype(&self, sf_bam_list: &str, sf_candidate_list: &str, extnd: i64, sf_out: &str) -> Result<()> {
        self.call_genotype_f(sf_bam_list, sf_candidate_list, extnd as f64, sf_out)
    }

    /// `call_genotype` with a (possibly float) extension
    pub fn call_genotype_f(&self, sf_bam_list: &str, sf_candidate_list: &str, extnd: f64, sf_out: &str) -> Result<()> {
        let mut l_bams = Vec::new();
        for line in utils::read_lines(sf_bam_list)? {
            let f = split_ws(&line);
            if let Some(b) = f.first() {
                l_bams.push(b.to_string());
            }
        }
        let mut l_sites: Vec<(String, i64)> = Vec::new();
        for line in utils::read_lines(sf_candidate_list)? {
            let f = split_ws(&line);
            if f.len() < 2 {
                continue;
            }
            l_sites.push((f[0].to_string(), utils::pint(f[1])?));
        }
        let pool = rayon::ThreadPoolBuilder::new().num_threads(self.n_jobs.max(1)).build().unwrap();
        for (i_cnt, sf_bam) in l_bams.iter().enumerate() {
            let sf_out_features = format!("{}{}.out", sf_out, i_cnt);
            let res: Vec<Result<String>> =
                pool.install(|| utils::par_map(l_sites.len(), |i| self.collect_features_one_site(&l_sites[i].0, l_sites[i].1, extnd, sf_bam)));
            let mut fout = utils::create(&sf_out_features)?;
            let mut consumed: HashSet<(String, i64)> = HashSet::new();
            for (i, r) in res.into_iter().enumerate() {
                let s = r?;
                if !consumed.insert(l_sites[i].clone()) {
                    continue;
                }
                fout.write_all(s.as_bytes())?;
            }
            fout.flush()?;
        }
        // merge the features of the different BAMs.  NOTE: as in the original, all but the first
        // merged field are taken from the *last candidate site* of the list (leaked loop
        // variables `chrm`/`pos`); this only matters with more than one BAM.
        let last_site = l_sites.last().cloned();
        let mut m_merged: SiteMap<GntpFeat> = IndexMap::new();
        for i_cnt in 0..l_bams.len() {
            let sf_out_features = format!("{}{}.out", sf_out, i_cnt);
            let m_info = Self::load_in_features_from_file(&sf_out_features)?;
            for (ins_chrm, m) in m_info.iter() {
                for (&ins_pos, rcd) in m.iter() {
                    let existing = m_merged.get(ins_chrm).and_then(|x| x.get(&ins_pos)).cloned();
                    if let Some(old) = existing {
                        let (lc, lp) = last_site.clone().ok_or_else(|| anyhow!("KeyError"))?;
                        let last = m_info.get(&lc).and_then(|x| x.get(&lp)).ok_or_else(|| anyhow!("KeyError: {}:{}", lc, lp))?;
                        let merged = GntpFeat {
                            n_af_clip: rcd.n_af_clip + old.n_af_clip,
                            n_full_map: last.n_full_map + old.n_full_map,
                            n_l_raw_clip: last.n_l_raw_clip + old.n_l_raw_clip,
                            n_r_raw_clip: last.n_r_raw_clip + old.n_r_raw_clip,
                            n_disc_pairs: last.n_disc_pairs + old.n_disc_pairs,
                            n_concd_pairs: last.n_concd_pairs + old.n_concd_pairs,
                            n_disc_large_indel: last.n_disc_large_indel + old.n_disc_large_indel,
                            s_clip_lens: format!("{}:{}", last.s_clip_lens, old.s_clip_lens),
                            n_polya: last.n_polya + old.n_polya,
                            n_disc_chrms: last.n_disc_chrms + old.n_disc_chrms,
                        };
                        m_merged.get_mut(ins_chrm).unwrap().insert(ins_pos, merged);
                    } else {
                        m_merged.entry(ins_chrm.clone()).or_default().insert(ins_pos, rcd.clone());
                    }
                }
            }
        }
        let mut fout = utils::create(sf_out)?;
        for (ins_chrm, m) in m_merged.iter() {
            for (ins_pos, rcd) in m.iter() {
                writeln!(fout, "{}\t{}\t{}", ins_chrm, ins_pos, rcd.fmt_fields())?;
            }
        }
        fout.flush()?;
        Ok(())
    }

    fn parse_feature_line(fields: &[&str]) -> Result<GntpFeat> {
        if fields.len() < 12 {
            return Err(anyhow!("IndexError: malformed genotype feature line"));
        }
        Ok(GntpFeat {
            n_af_clip: utils::pint(fields[2])?,
            n_full_map: utils::pint(fields[3])?,
            n_l_raw_clip: utils::pint(fields[4])?,
            n_r_raw_clip: utils::pint(fields[5])?,
            n_disc_pairs: utils::pint(fields[6])?,
            n_concd_pairs: utils::pint(fields[7])?,
            n_disc_large_indel: utils::pint(fields[8])?,
            s_clip_lens: fields[9].to_string(),
            n_polya: utils::pint(fields[10])?,
            n_disc_chrms: utils::pint(fields[11])?,
        })
    }

    /// `load_in_features_from_file`
    pub fn load_in_features_from_file(sf_features: &str) -> Result<SiteMap<GntpFeat>> {
        let mut m: SiteMap<GntpFeat> = IndexMap::new();
        for line in utils::read_lines(sf_features)? {
            let fields = split_ws(&line);
            if fields.is_empty() {
                continue;
            }
            let chrm = fields[0].to_string();
            let pos = utils::pint(fields[1])?;
            let g = Self::parse_feature_line(&fields)?;
            let e = m.entry(chrm).or_default();
            if !e.contains_key(&pos) {
                e.insert(pos, g);
            }
        }
        Ok(m)
    }

    /// `load_in_features_from_file_with_cov_cutoff`: skip sites with more clipped reads (at the
    /// exact breakpoint) than `i_cov_cutoff`.
    pub fn load_in_features_from_file_with_cov_cutoff(sf_features: &str, i_cov_cutoff: f64) -> Result<SiteMap<GntpFeat>> {
        let mut m: SiteMap<GntpFeat> = IndexMap::new();
        for line in utils::read_lines(sf_features)? {
            let fields = split_ws(&line);
            if fields.is_empty() {
                continue;
            }
            let chrm = fields[0].to_string();
            let pos = utils::pint(fields[1])?;
            if fields.len() < 10 {
                return Err(anyhow!("IndexError: malformed genotype feature line"));
            }
            let n_clip_len = fields[9].split(':').count();
            if (n_clip_len as f64) > i_cov_cutoff {
                continue;
            }
            let g = Self::parse_feature_line(&fields)?;
            let e = m.entry(chrm).or_default();
            if !e.contains_key(&pos) {
                e.insert(pos, g);
            }
        }
        Ok(m)
    }

    /// `filter_by_af`
    pub fn filter_by_af(m_features: &SiteMap<GntpFeat>, f_af_min: f64, f_af_max: f64, sf_out: &str) -> Result<()> {
        let mut fout = utils::create(sf_out)?;
        for (chrm, m) in m_features.iter() {
            for (pos, g) in m.iter() {
                let mut f_clip_af = 0.0;
                if g.n_full_map + g.n_af_clip != 0 {
                    f_clip_af = g.n_af_clip as f64 / (g.n_full_map + g.n_af_clip) as f64;
                }
                let mut f_disc_af = 0.0;
                let n_all_pairs = g.n_disc_pairs + g.n_concd_pairs;
                if n_all_pairs != 0 {
                    f_disc_af = g.n_disc_pairs as f64 / n_all_pairs as f64;
                }
                if (f_clip_af >= f_af_min && f_clip_af <= f_af_max) && (f_disc_af >= f_af_min && f_disc_af <= f_af_max) {
                    writeln!(
                        fout,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        chrm, pos, g.n_af_clip, g.n_full_map, g.n_l_raw_clip, g.n_r_raw_clip, g.n_disc_pairs, g.n_concd_pairs, g.n_disc_large_indel, g.s_clip_lens
                    )?;
                }
            }
        }
        Ok(())
    }
}

