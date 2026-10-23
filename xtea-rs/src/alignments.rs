//! Port of `x_alignments.py` (BamInfo, XAlignment, XBamInfo helpers).

use crate::annotation::{process_chrm_name, XAnnotation};
use crate::bam::{AlignmentFile, PyAln, Record};
use crate::chrom::is_decoy_contig_chrms;
use crate::cluster::form_one_side_cluster;
use crate::globals as gv;
use crate::utils;
use anyhow::Result;
use indexmap::IndexMap;
use std::collections::HashMap;

/// `BamInfo`
pub struct BamInfo {
    pub sf_bam: String,
    pub sf_reference: String,
}

impl BamInfo {
    pub fn new(sf_bam: &str, sf_ref: &str) -> Self {
        BamInfo { sf_bam: sf_bam.to_string(), sf_reference: sf_ref.to_string() }
    }

    fn open(&self) -> Result<AlignmentFile> {
        AlignmentFile::open(&self.sf_bam, &self.sf_reference)
    }

    /// `index_reference_name_id`
    pub fn index_reference_name_id(&self) -> Result<HashMap<i32, String>> {
        let f = self.open()?;
        let mut m = HashMap::new();
        for (i, n) in f.references().iter().enumerate() {
            m.insert(i as i32, n.clone());
        }
        Ok(m)
    }

    /// `get_all_reference_names` (ordered like the header)
    pub fn get_all_reference_names(&self) -> Result<Vec<String>> {
        let f = self.open()?;
        // dict semantics: later duplicates keep the first position
        let mut seen = std::collections::HashSet::new();
        let mut v = Vec::new();
        for n in f.references() {
            if seen.insert(n.clone()) {
                v.push(n.clone());
            }
        }
        Ok(v)
    }

    /// `get_all_chrom_name_length`
    pub fn get_all_chrom_name_length(&self) -> Result<IndexMap<String, i64>> {
        let f = self.open()?;
        let mut m = IndexMap::new();
        for (n, l) in f.references().iter().zip(f.lengths().iter()) {
            m.insert(n.clone(), *l);
        }
        Ok(m)
    }

    /// `is_chrm_contain_chr`
    pub fn is_chrm_contain_chr(&self) -> Result<bool> {
        let f = self.open()?;
        Ok(f.references().iter().any(|n| n == "chr1"))
    }

    pub fn process_chrm_name(&self, chrm: &str, b_with_chr: bool) -> String {
        process_chrm_name(chrm, b_with_chr)
    }
}

/// Result of `BamInfo.cnt_discordant_pairs`: (n_in_rep, n_raw, (is_cluster, chrm, pos)).
pub type DiscCnt = (i64, i64, (bool, Option<String>, Option<i64>));

/// `BamInfo.cnt_discordant_pairs`.  NOTE: the distance test uses the window `start`, exactly as
/// in the original.
pub fn cnt_discordant_pairs(
    bamfile: &mut AlignmentFile,
    m_chrm_ids: &HashMap<i32, String>,
    chrm: &str,
    start: i64,
    end: i64,
    i_is: i64,
    f_dev: i64,
    xannotation: &XAnnotation,
) -> Result<DiscCnt> {
    let mut n_cnt: i64 = 0;
    let mut n_raw_cnt: i64 = 0;
    let mut i_max_is: i64 = 2500;
    let tmp = crate::py::fmt::py_int((i_is as f64) + 3.0 * (f_dev as f64));
    if tmp > i_max_is {
        i_max_is = tmp;
    }
    let mut m_mate_pos: IndexMap<String, Vec<i64>> = IndexMap::new();
    let min_mapq = gv::minimum_disc_mapq();
    bamfile.fetch(chrm, start, end, |a: &Record| {
        if a.is_dup() || a.is_supp() {
            return;
        }
        if a.is_unmap() || a.mate_is_unmapped() {
            return;
        }
        if a.is_sec() {
            return;
        }
        if a.mapping_quality() < min_mapq {
            return;
        }
        let ntid = a.next_reference_id();
        if ntid < 0 {
            return;
        }
        let mate_chrm = match m_chrm_ids.get(&ntid) {
            Some(n) => n.clone(),
            None => return,
        };
        if is_decoy_contig_chrms(&mate_chrm) {
            return;
        }
        let mate_pos = a.next_reference_start();
        if chrm != mate_chrm || (mate_pos - start).abs() > i_is {
            if chrm != mate_chrm || (mate_pos - start).abs() > i_max_is {
                m_mate_pos.entry(mate_chrm.clone()).or_default().push(mate_pos);
                n_raw_cnt += 1;
            }
            if xannotation.in_repeat(&mate_chrm, mate_pos) {
                n_cnt += 1;
            }
        }
    })?;
    let cl = form_one_side_cluster(&mut m_mate_pos, i_is, gv::MIN_RAW_DISC_CLUSTER_RATIO);
    Ok((n_cnt, n_raw_cnt, cl))
}

/// `XAlignment._cvt_to_Ascii_quality`
pub fn cvt_to_ascii_quality(q: &[u8]) -> String {
    q.iter().map(|&x| char::from_u32(x as u32 + 33).unwrap_or('?')).collect()
}

/// `XAlignment._change_barcode_format`: "xxx-1" -> "xxx_1"
pub fn change_barcode_format(s: &str) -> String {
    s.split('-').collect::<Vec<_>>().join("_")
}

/// `XAlignment._cvt_algnmt_to_fq`
pub fn cvt_algnmt_to_fq(a: &Record) -> String {
    let mut name = a.query_name();
    name.push_str(if a.is_read1() { "_1" } else { "_2" });
    if let Some(bx) = a.get_tag_string("BX") {
        if !bx.is_empty() {
            name.push_str(&format!(" BX:{}", bx));
        }
    }
    let seq = a.query_sequence().unwrap_or_default();
    let qual = cvt_to_ascii_quality(&a.query_qualities().unwrap_or_default());
    format!("@{}\n{}\n+\n{}\n", name, seq, qual)
}

/// `XAlignment._cvt_algnmt_to_fa`
pub fn cvt_algnmt_to_fa(a: &Record) -> String {
    let mut name = a.query_name();
    name.push_str(if a.is_read1() { "_1" } else { "_2" });
    if let Some(bx) = a.get_tag_string("BX") {
        if !bx.is_empty() {
            name.push_str(&format!(" BX:{}", bx));
        }
    }
    let seq = a.query_sequence().unwrap_or_default();
    format!(">{}\n{}\n", name, seq)
}

/// `XAlignment.is_discordant_pair_no_unmap`
pub fn is_discordant_pair_no_unmap(chrm: &str, mate_chrm: &str, template_lth: i64, i_is: i64, f_dev: i64) -> bool {
    if chrm == "*" || mate_chrm == "*" {
        return false;
    }
    if mate_chrm != chrm {
        return true;
    }
    let i_dist = template_lth.abs();
    let i_is1 = i_is + 3 * f_dev;
    i_dist > i_is1
}

/// `XAlignment.is_TE_caused_discordant_pair`
pub fn is_te_caused_discordant_pair(chrm: &str, map_pos: i64, mate_chrm: &str, mate_map_pos: i64, xann: &XAnnotation) -> (bool, bool) {
    let (a, _) = xann.is_within_repeat_region(chrm, map_pos);
    let (b, _) = xann.is_within_repeat_region(mate_chrm, mate_map_pos);
    (a || b, a && b)
}

/// `OneAlignment.is_full_map`
pub fn one_alignment_is_full_map(cigar: &[(u32, u32)], min_clip_len: i64) -> bool {
    if cigar.is_empty() {
        return false;
    }
    if cigar.len() == 1 && cigar[0].0 == 0 {
        return true;
    }
    if (cigar[0].0 == 4 || cigar[0].0 == 5) && (cigar[0].1 as i64) <= min_clip_len {
        return true;
    }
    let l = cigar[cigar.len() - 1];
    if (l.0 == 4 || l.0 == 5) && (l.1 as i64) <= min_clip_len {
        return true;
    }
    false
}

/// `BarcodeAlignment` accessors.
pub fn barcode_get_chrm(a: &Record) -> String {
    a.get_tag_string("ch").unwrap_or_else(|| "*".to_string())
}
pub fn barcode_get_mate_chrm(a: &Record) -> String {
    a.get_tag_string("nh").unwrap_or_else(|| "*".to_string())
}
pub fn barcode_get_map_pos(a: &Record) -> i64 {
    a.get_tag_string("cp").and_then(|s| s.parse::<i64>().ok()).unwrap_or(0)
}

/// `BamInfo.load_read_names`
pub fn load_read_names(sf_names: &str) -> Result<Vec<String>> {
    let mut v: Vec<String> = utils::read_lines(sf_names)?;
    if let Some(i) = v.iter().position(|s| s.is_empty()) {
        v.remove(i);
    }
    Ok(v)
}
