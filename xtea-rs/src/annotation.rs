//! Port of `x_annotation.py` (RepeatMasker / bed annotation loading and queries).

use crate::fasta::{reverse_complement_xtea, FastaFile};
use crate::globals as gv;
use crate::intervaltree::IntervalTree;
use crate::utils::{self, split_ws};
use anyhow::{anyhow, Result};
use indexmap::IndexMap;
use std::collections::HashMap;
use std::io::Write;
use std::sync::Arc;

/// One annotation record.  Different Python loaders store tuples of different shapes; the
/// `Layout` of the owning `XAnnotation` tells which fields are meaningful.
#[derive(Clone, Debug)]
pub struct RmskRec {
    pub end: i64,
    pub div: f64,
    pub rc: bool,
    pub sub_type: Arc<str>,
    pub family: Arc<str>,
    pub csn_start: Option<i64>,
    pub csn_end: Option<i64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Layout {
    /// (end_pos, b_rc, sub_type, csn_start, csn_end)
    Basic,
    /// (end_pos, mismatch_rate, b_rc, sub_type, csn_start, csn_end)
    Div6,
    /// (end_pos, mismatch_rate, b_rc, sub_type, family, csn_start, csn_end)
    Div7,
}

#[derive(Default)]
struct Interner {
    map: HashMap<String, Arc<str>>,
}

impl Interner {
    fn get(&mut self, s: &str) -> Arc<str> {
        if let Some(a) = self.map.get(s) {
            return a.clone();
        }
        let a: Arc<str> = Arc::from(s);
        self.map.insert(s.to_string(), a.clone());
        a
    }
}

/// `XAnnotation`
pub struct XAnnotation {
    pub sf_annotation: String,
    /// m_rmsk_annotation[chrm][start_pos] = [records]
    pub m_rmsk_annotation: IndexMap<String, IndexMap<i64, Vec<RmskRec>>>,
    pub m_rmsk_chrm_regions: HashMap<String, Vec<i64>>,
    pub b_with_chr: bool,
    pub m_interval_tree: HashMap<String, IntervalTree>,
    pub layout: Layout,
    interner: Interner,
}

/// `_process_chrm_name` shared logic (also used by XBlackList and others).
pub fn process_chrm_name(chrm: &str, b_with_chr: bool) -> String {
    let b_chrm_with_chr = chrm.len() > 3 && chrm.starts_with("chr");
    if b_with_chr && b_chrm_with_chr {
        chrm.to_string()
    } else if b_with_chr && !b_chrm_with_chr {
        format!("chr{}", chrm)
    } else if !b_with_chr && b_chrm_with_chr {
        chrm[3..].to_string()
    } else {
        chrm.to_string()
    }
}

fn parse_i(s: &str) -> Result<i64> {
    utils::pint(s)
}

impl XAnnotation {
    pub fn new(sf_rmsk: &str) -> Self {
        XAnnotation {
            sf_annotation: sf_rmsk.to_string(),
            m_rmsk_annotation: IndexMap::new(),
            m_rmsk_chrm_regions: HashMap::new(),
            b_with_chr: true,
            m_interval_tree: HashMap::new(),
            layout: Layout::Basic,
            interner: Interner::default(),
        }
    }

    pub fn set_with_chr(&mut self, b: bool) {
        self.b_with_chr = b;
    }

    pub fn process_chrm_name(&self, chrm: &str) -> String {
        process_chrm_name(chrm, self.b_with_chr)
    }

    pub fn is_ref_chrm_with_chr(&self, sf_ref: &str) -> Result<bool> {
        let fa = FastaFile::open(sf_ref)?;
        Ok(fa.has_reference("chr1"))
    }

    fn push(&mut self, chrm: String, start: i64, rec: RmskRec) {
        let m = self.m_rmsk_annotation.entry(chrm.clone()).or_default();
        if m.contains_key(&start) {
            // matches the "has more than 1 annotation" message in the original (printed for
            // the un-extended start position check)
        }
        m.entry(start).or_default().push(rec);
    }

    fn warn_dup(&self, chrm: &str, start_pos: i64) {
        if let Some(m) = self.m_rmsk_annotation.get(chrm) {
            if m.contains_key(&start_pos) {
                println!("Position {}:{} has more than  1 annotation!!!!", chrm, start_pos);
            }
        }
    }

    /// Parse the common RepeatMasker columns.
    fn parse_rmsk_fields(&mut self, fields: &[&str]) -> Result<(String, i64, i64, bool, Arc<str>, Option<i64>, Option<i64>)> {
        if fields.len() < 14 {
            return Err(anyhow!("malformed RepeatMasker line: {}", fields.join(" ")));
        }
        let chrm = self.process_chrm_name(fields[4]);
        let start_pos = parse_i(fields[5])?;
        let end_pos = parse_i(fields[6])?;
        let b_rc = fields[8] == "C";
        let sub_type = self.interner.get(fields[9]);
        let (cs, ce) = if b_rc {
            (parse_i(fields[13])?, parse_i(fields[12])?)
        } else {
            (parse_i(fields[11])?, parse_i(fields[12])?)
        };
        Ok((chrm, start_pos, end_pos, b_rc, sub_type, Some(cs), Some(ce)))
    }

    fn empty_arc(&mut self) -> Arc<str> {
        self.interner.get("")
    }

    /// `load_rmsk_annotation`
    pub fn load_rmsk_annotation(&mut self) -> Result<()> {
        self.layout = Layout::Basic;
        let lines = utils::read_lines(&self.sf_annotation)?;
        let fam = self.empty_arc();
        for line in lines {
            let fields = split_ws(&line);
            let (chrm, start_pos, end_pos, b_rc, sub_type, cs, ce) = self.parse_rmsk_fields(&fields)?;
            self.warn_dup(&chrm, start_pos);
            let extd = start_pos - gv::load_rmsk_left_extnd();
            self.push(chrm, extd, RmskRec { end: end_pos, div: 0.0, rc: b_rc, sub_type, family: fam.clone(), csn_start: cs, csn_end: ce });
        }
        Ok(())
    }

    /// `load_rmsk_annotation_no_extnd`
    pub fn load_rmsk_annotation_no_extnd(&mut self) -> Result<()> {
        self.layout = Layout::Basic;
        let lines = utils::read_lines(&self.sf_annotation)?;
        let fam = self.empty_arc();
        for line in lines {
            let fields = split_ws(&line);
            let (chrm, start_pos, end_pos, b_rc, sub_type, cs, ce) = self.parse_rmsk_fields(&fields)?;
            self.warn_dup(&chrm, start_pos);
            self.push(chrm, start_pos, RmskRec { end: end_pos, div: 0.0, rc: b_rc, sub_type, family: fam.clone(), csn_start: cs, csn_end: ce });
        }
        Ok(())
    }

    /// `load_annotation_no_extnd_from_bed`
    pub fn load_annotation_no_extnd_from_bed(&mut self) -> Result<()> {
        self.layout = Layout::Basic;
        let lines = utils::read_lines(&self.sf_annotation)?;
        let fam = self.empty_arc();
        for line in lines {
            let fields = split_ws(&line);
            if fields.len() < 3 {
                return Err(anyhow!("malformed bed line: {}", line));
            }
            let chrm = self.process_chrm_name(fields[0]);
            let start_pos = parse_i(fields[1])?;
            let end_pos = parse_i(fields[2])?;
            let s_rc = fields[fields.len() - 2];
            let b_rc = s_rc == "C";
            let mut sub = fields[fields.len() - 1].to_string();
            if !sub.contains("SVA") {
                sub = String::new();
            }
            let sub_type = self.interner.get(&sub);
            self.warn_dup(&chrm, start_pos);
            self.push(chrm, start_pos, RmskRec { end: end_pos, div: 0.0, rc: b_rc, sub_type, family: fam.clone(), csn_start: None, csn_end: None });
        }
        Ok(())
    }

    /// `load_rmsk_annotation_with_extnd_with_lenth_cutoff`
    pub fn load_rmsk_annotation_with_extnd_with_lenth_cutoff(&mut self, i_extnd: i64, i_min_len: i64) -> Result<()> {
        self.layout = Layout::Basic;
        let lines = utils::read_lines(&self.sf_annotation)?;
        let fam = self.empty_arc();
        for line in lines {
            let fields = split_ws(&line);
            let (chrm, start_pos, end_pos, b_rc, sub_type, cs, ce) = self.parse_rmsk_fields(&fields)?;
            if (end_pos - start_pos).abs() < i_min_len {
                continue;
            }
            self.warn_dup(&chrm, start_pos);
            self.push(
                chrm,
                start_pos - i_extnd,
                RmskRec { end: end_pos + i_extnd, div: 0.0, rc: b_rc, sub_type, family: fam.clone(), csn_start: cs, csn_end: ce },
            );
        }
        Ok(())
    }

    /// `load_rmsk_annotation_with_divgnt`
    pub fn load_rmsk_annotation_with_divgnt(&mut self) -> Result<()> {
        self.layout = Layout::Div6;
        let lines = utils::read_lines(&self.sf_annotation)?;
        let fam = self.empty_arc();
        for line in lines {
            let fields = split_ws(&line);
            if fields.len() < 2 {
                return Err(anyhow!("malformed RepeatMasker line: {}", line));
            }
            let mismatch = utils::pfloat(fields[1])?;
            let (chrm, start_pos, end_pos, b_rc, sub_type, cs, ce) = self.parse_rmsk_fields(&fields)?;
            self.warn_dup(&chrm, start_pos);
            self.push(chrm, start_pos - 100, RmskRec { end: end_pos, div: mismatch, rc: b_rc, sub_type, family: fam.clone(), csn_start: cs, csn_end: ce });
        }
        Ok(())
    }

    fn skip_rmsk_line(line: &str, fields: &[&str]) -> bool {
        if utils::rstrip(line).is_empty() {
            return true;
        }
        if fields.len() < 15 {
            return true;
        }
        let f0 = utils::rstrip(fields[0]);
        f0 == "SW" || f0 == "score"
    }

    /// `load_rmsk_annotation_with_divgnt_extnd`
    pub fn load_rmsk_annotation_with_divgnt_extnd(&mut self, iextnd: i64) -> Result<()> {
        self.layout = Layout::Div7;
        let lines = utils::read_lines(&self.sf_annotation)?;
        for line in lines {
            let fields = split_ws(&line);
            if Self::skip_rmsk_line(&line, &fields) {
                continue;
            }
            let mismatch = utils::pfloat(fields[1])?;
            let (chrm, s0, e0, b_rc, sub_type, cs, ce) = self.parse_rmsk_fields(&fields)?;
            let start_pos = s0 - iextnd;
            let end_pos = e0 + iextnd;
            let family = self.interner.get(fields[10]);
            self.warn_dup(&chrm, start_pos);
            self.push(chrm, start_pos - 100, RmskRec { end: end_pos, div: mismatch, rc: b_rc, sub_type, family, csn_start: cs, csn_end: ce });
        }
        Ok(())
    }

    /// `load_rmsk_annotation_with_divgnt_no_slack_with_extnd`
    pub fn load_rmsk_annotation_with_divgnt_no_slack_with_extnd(&mut self, iextnd: i64) -> Result<()> {
        self.layout = Layout::Div7;
        let lines = utils::read_lines(&self.sf_annotation)?;
        for line in lines {
            let fields = split_ws(&line);
            if Self::skip_rmsk_line(&line, &fields) {
                continue;
            }
            let mismatch = utils::pfloat(fields[1])?;
            let (chrm, s0, e0, b_rc, sub_type, cs, ce) = self.parse_rmsk_fields(&fields)?;
            let start_pos = s0 - iextnd;
            let end_pos = e0 + iextnd;
            let family = self.interner.get(fields[10]);
            self.warn_dup(&chrm, start_pos);
            self.push(chrm, start_pos, RmskRec { end: end_pos, div: mismatch, rc: b_rc, sub_type, family, csn_start: cs, csn_end: ce });
        }
        Ok(())
    }

    /// `load_rmsk_annotation_with_extnd_div_with_lenth_cutoff`
    pub fn load_rmsk_annotation_with_extnd_div_with_lenth_cutoff(&mut self, i_extnd: i64, i_min_len: i64) -> Result<()> {
        self.layout = Layout::Div7;
        let lines = utils::read_lines(&self.sf_annotation)?;
        for line in lines {
            let fields = split_ws(&line);
            if Self::skip_rmsk_line(&line, &fields) {
                continue;
            }
            let mismatch = utils::pfloat(fields[1])?;
            let (chrm, start_pos, end_pos, b_rc, sub_type, cs, ce) = self.parse_rmsk_fields(&fields)?;
            let family = self.interner.get(fields[10]);
            if (end_pos - start_pos).abs() < i_min_len {
                continue;
            }
            self.warn_dup(&chrm, start_pos);
            self.push(
                chrm,
                start_pos - i_extnd,
                RmskRec { end: end_pos + i_extnd, div: mismatch, rc: b_rc, sub_type, family, csn_start: cs, csn_end: ce },
            );
        }
        Ok(())
    }

    /// `load_rmsk_annotation_L1`
    pub fn load_rmsk_annotation_l1(&mut self) -> Result<()> {
        self.layout = Layout::Basic;
        let lines = utils::read_lines(&self.sf_annotation)?;
        let fam = self.empty_arc();
        for line in lines {
            let fields = split_ws(&line);
            if fields.len() < 11 {
                return Err(anyhow!("malformed RepeatMasker line: {}", line));
            }
            if fields[10] != "LINE/L1" {
                continue;
            }
            let (chrm, start_pos, end_pos, b_rc, sub_type, cs, ce) = self.parse_rmsk_fields(&fields)?;
            self.warn_dup(&chrm, start_pos);
            self.push(chrm, start_pos - 100, RmskRec { end: end_pos, div: 0.0, rc: b_rc, sub_type, family: fam.clone(), csn_start: cs, csn_end: ce });
        }
        Ok(())
    }

    /// `index_rmsk_annotation`: sorted start positions per chromosome (binary search index).
    pub fn index_rmsk_annotation(&mut self) {
        for (chrm, m) in self.m_rmsk_annotation.iter() {
            let mut v: Vec<i64> = m.keys().cloned().collect();
            v.sort();
            self.m_rmsk_chrm_regions.insert(chrm.clone(), v);
        }
    }

    /// `is_within_repeat_region` (binary search version).  NOTE: the chromosome name is not
    /// processed here, exactly like the original.
    pub fn is_within_repeat_region(&self, chrm: &str, pos: i64) -> (bool, i64) {
        let regions = match self.m_rmsk_chrm_regions.get(chrm) {
            Some(r) => r,
            None => return (false, -1),
        };
        let ann = &self.m_rmsk_annotation[chrm];
        let mut lo: i64 = 0;
        let mut hi: i64 = regions.len() as i64 - 1;
        while lo <= hi {
            // Python 3: (lo + hi) / 2 is a float -> list index would fail; the original
            // Python 2 semantics (integer division) are the only working ones.
            let mid = (lo + hi) / 2;
            let mid_start = regions[mid as usize];
            let mid_end = ann[&mid_start][0].end;
            if pos >= mid_start && pos <= mid_end {
                return (true, mid_start);
            } else if pos > mid_start {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        (false, -1)
    }

    /// `index_rmsk_annotation_interval_tree`
    pub fn index_rmsk_annotation_interval_tree(&mut self) {
        for (chrm, m) in self.m_rmsk_annotation.iter() {
            let mut t = IntervalTree::new();
            for (pos, recs) in m.iter() {
                t.addi(*pos, recs[0].end);
            }
            self.m_interval_tree.insert(chrm.clone(), t);
        }
    }

    /// `is_within_repeat_region_interval_tree`: (hit, start of the first interval in CPython
    /// set order).
    pub fn is_within_repeat_region_interval_tree(&self, chrm1: &str, pos: i64) -> (bool, i64) {
        let chrm = self.process_chrm_name(chrm1);
        let t = match self.m_interval_tree.get(&chrm) {
            Some(t) => t,
            None => return (false, -1),
        };
        let hits = t.at(pos);
        match hits.first() {
            Some((b, _)) => (true, *b),
            None => (false, -1),
        }
    }

    /// Boolean-only version (cheaper when the start is not needed).
    pub fn in_repeat(&self, chrm1: &str, pos: i64) -> bool {
        self.is_within_repeat_region_interval_tree(chrm1, pos).0
    }

    /// `is_within_repeat_region_interval_tree2`: all hits (start positions) in set order.
    pub fn is_within_repeat_region_interval_tree2(&self, chrm1: &str, pos: i64) -> Vec<(bool, i64)> {
        let chrm = self.process_chrm_name(chrm1);
        match self.m_interval_tree.get(&chrm) {
            Some(t) => t.at(pos).into_iter().map(|(b, _)| (true, b)).collect(),
            None => Vec::new(),
        }
    }

    /// `get_div_subfamily` -> (div_rate, sub_family, family, pos, pos_end); None when absent.
    pub fn get_div_subfamily(&self, chrm1: &str, pos: i64) -> Option<(f64, String, String, i64, i64)> {
        let chrm = self.process_chrm_name(chrm1);
        let m = self.m_rmsk_annotation.get(&chrm)?;
        let recs = m.get(&pos)?;
        let r = &recs[0];
        Some((r.div, r.sub_type.to_string(), r.family.to_string(), pos, r.end))
    }

    /// `get_div_subfamily_with_min_div`
    pub fn get_div_subfamily_with_min_div(&self, chrm1: &str, pos: i64) -> Option<(f64, String, String, i64, i64)> {
        let chrm = self.process_chrm_name(chrm1);
        let m = self.m_rmsk_annotation.get(&chrm)?;
        let recs = m.get(&pos)?;
        let mut div = recs[0].div;
        let mut sub = recs[0].sub_type.clone();
        let mut fam = recs[0].family.clone();
        let mut end = recs[0].end;
        for r in recs.iter().skip(1) {
            if r.div < div {
                div = r.div;
                sub = r.sub_type.clone();
                fam = r.family.clone();
                end = r.end;
            }
        }
        Some((div, sub.to_string(), fam.to_string(), pos, end))
    }

    /// `get_div_subfamily_with_min_div2`: all records at the start position.
    pub fn get_div_subfamily_with_min_div2(&self, chrm1: &str, pos: i64) -> Vec<(f64, String, String, i64, i64)> {
        let chrm = self.process_chrm_name(chrm1);
        let mut v = Vec::new();
        if let Some(m) = self.m_rmsk_annotation.get(&chrm) {
            if let Some(recs) = m.get(&pos) {
                for r in recs {
                    v.push((r.div, r.sub_type.to_string(), r.family.to_string(), pos, r.end));
                }
            }
        }
        v
    }

    pub fn get_annotation_te_region_lenth(&self, chrm: &str, start_pos: i64) -> i64 {
        self.m_rmsk_annotation[chrm][&start_pos][0].end - start_pos
    }

    /// `collect_seqs_of_TE_from_ref`
    pub fn collect_seqs_of_te_from_ref(&self, sf_ref: &str, flank_lth: i64, b_with_flank: bool, sf_out_fa: &str) -> Result<()> {
        let fa = FastaFile::open(sf_ref)?;
        let refs: std::collections::HashSet<&String> = fa.references().iter().collect();
        let mut fout = utils::create(sf_out_fa)?;
        for (chrm, m) in self.m_rmsk_annotation.iter() {
            if !refs.contains(chrm) {
                continue;
            }
            for (start_pos0, recs) in m.iter() {
                let rec = &recs[0];
                let mut start_pos = *start_pos0;
                let end_pos = rec.end;
                let mut pstart = start_pos;
                let mut pend = end_pos;
                if b_with_flank {
                    pstart = start_pos - flank_lth;
                    pend = end_pos + flank_lth;
                }
                if start_pos < 0 {
                    start_pos = 0;
                }
                if pstart < 0 {
                    continue;
                }
                let mut s_reg = fa.fetch(chrm, pstart, pend)?;
                if rec.rc {
                    s_reg = reverse_complement_xtea(&s_reg);
                }
                writeln!(fout, ">{}{}{}{}{}{}{}", chrm, gv::S_DELIM, start_pos, gv::S_DELIM, end_pos, gv::S_DELIM, rec.sub_type)?;
                writeln!(fout, "{}", s_reg)?;
            }
        }
        fout.flush()?;
        Ok(())
    }

    /// `collect_seqs_of_TE_from_ref_bed_fmt`
    pub fn collect_seqs_of_te_from_ref_bed_fmt(&self, sf_ref: &str, sf_out_fa: &str, flank_lth: i64) -> Result<()> {
        let fa = FastaFile::open(sf_ref)?;
        let refs: std::collections::HashSet<&String> = fa.references().iter().collect();
        let mut fout = utils::create(sf_out_fa)?;
        for line in utils::read_lines(&self.sf_annotation)? {
            let fields = split_ws(&line);
            if fields.len() < 3 {
                continue;
            }
            let chrm = fields[0].to_string();
            let pstart = parse_i(fields[1])?;
            let pend = parse_i(fields[2])?;
            let mut s_type = fields[fields.len() - 1].to_string();
            if !s_type.contains("SVA") {
                s_type = String::new();
            }
            if !refs.contains(&chrm) {
                continue;
            }
            if pend <= 0 {
                continue;
            }
            let s_reg = fa.fetch(&chrm, pstart - flank_lth, pend + flank_lth)?;
            writeln!(fout, ">{}{}{}{}{}{}{}", chrm, gv::S_DELIM, pstart, gv::S_DELIM, pend, gv::S_DELIM, s_type)?;
            writeln!(fout, "{}", s_reg)?;
        }
        fout.flush()?;
        Ok(())
    }

    /// `collect_flank_regions_of_TE_from_ref`
    pub fn collect_flank_regions_of_te_from_ref(&self, sf_ref: &str, flank_lth: i64, sf_out_fa: &str) -> Result<()> {
        let fa = FastaFile::open(sf_ref)?;
        let refs: std::collections::HashSet<&String> = fa.references().iter().collect();
        let mut fout = utils::create(sf_out_fa)?;
        for (chrm, m) in self.m_rmsk_annotation.iter() {
            if !refs.contains(chrm) {
                continue;
            }
            for (start_pos, recs) in m.iter() {
                let rec = &recs[0];
                let start_pos = *start_pos;
                let end_pos = rec.end;
                let bi_rc = if rec.rc { 1 } else { 0 };
                let lstart = start_pos - flank_lth;
                if lstart < 0 {
                    continue;
                }
                let s_lreg = fa.fetch(chrm, lstart, start_pos)?;
                writeln!(fout, ">{}{}{}{}{}{}{}{}{}L", chrm, gv::S_DELIM, start_pos, gv::S_DELIM, end_pos, gv::S_DELIM, rec.sub_type, gv::S_DELIM, bi_rc)?;
                writeln!(fout, "{}", s_lreg)?;
                let rstart = end_pos;
                if rstart <= 0 {
                    continue;
                }
                let s_rreg = fa.fetch(chrm, rstart, end_pos + flank_lth)?;
                writeln!(fout, ">{}{}{}{}{}{}{}{}{}R", chrm, gv::S_DELIM, start_pos, gv::S_DELIM, end_pos, gv::S_DELIM, rec.sub_type, gv::S_DELIM, bi_rc)?;
                writeln!(fout, "{}", s_rreg)?;
            }
        }
        fout.flush()?;
        Ok(())
    }

    /// `bwa_index_TE_seqs`
    pub fn bwa_index_te_seqs(&self, sf_fa: &str) {
        let cmd = format!("{} index {}", gv::BWA_PATH, sf_fa);
        utils::run_cmd_small_output(&cmd);
    }

    /// `select_L1_annotation_from_rmsk`
    pub fn select_l1_annotation_from_rmsk(&self, lth_cutoff: i64, sf_out: &str) -> Result<()> {
        let mut fout = utils::create(sf_out)?;
        for line in utils::read_raw_lines(&self.sf_annotation)? {
            let fields = split_ws(&line);
            if Self::skip_rmsk_line(&line, &fields) {
                continue;
            }
            let b_rc = fields[8] == "C";
            let sub_type = fields[9];
            if fields[10] != "LINE/L1" {
                continue;
            }
            let (cs, ce) = if b_rc { (parse_i(fields[13])?, parse_i(fields[12])?) } else { (parse_i(fields[11])?, parse_i(fields[12])?) };
            if (ce - cs).abs() >= lth_cutoff || sub_type == "L1HS" {
                fout.write_all(line.as_bytes())?;
            }
        }
        fout.flush()?;
        Ok(())
    }
}
