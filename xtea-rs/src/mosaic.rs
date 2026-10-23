//! Mosaic insertion calling from high-coverage bulk data (port of `x_mosaic_calling.py`) and
//! the joint filtering of mosaic calls over several samples (port of `x_joint_calling.py`).

use crate::intermediate_sites::SiteMap;
use crate::post_filter::{AFConflictFilter, XPostFilter, XTEARsltParser, XTPRTFilter, ONE_HALF_SIDE, ONE_SIDE, TWO_SIDE};
use crate::py::text::py_int;
use crate::somatic::SomaticMeiCaller;
use crate::utils::{self, split_ws};
use anyhow::{anyhow, Result};
use indexmap::IndexMap;
use std::collections::HashMap;
use std::io::Write;

/// `MosaicCaller` (derived from `XPostFilter`)
pub struct MosaicCaller {
    pub post: XPostFilter,
    pub brkpnt_slack: i64,
}

/// Which repeat-specific loop of `run_call_mosaic` to use.
#[derive(Clone, Copy, PartialEq)]
enum MosaicRep {
    L1,
    Alu,
    Sva,
}

impl MosaicCaller {
    pub fn new(swfolder: &str, n_jobs: usize) -> Self {
        MosaicCaller { post: XPostFilter::new(swfolder, n_jobs), brkpnt_slack: 100 }
    }

    /// `get_AF_cutoff`: (upper, lower) allele-fraction bounds per support type.
    pub fn get_af_cutoff(l_types: &[&str]) -> HashMap<String, (f64, f64)> {
        l_types.iter().map(|t| (t.to_string(), (0.35, 0.02))).collect()
    }

    /// `run_call_mosaic`
    pub fn run_call_mosaic(&self, sf_xtea_rslt: &str, sf_rmsk: &str, i_min_copy_len: i64, i_rep_type: i64, sf_black_list: &str, sf_new_out: &str) -> Result<()> {
        let l_old_rcd = XTEARsltParser::load_in_xtea_rslt(sf_xtea_rslt)?;
        // b_rslt_with_chr is True, and the default rmsk boundary extension (100) is used
        let xannotation = XPostFilter::construct_interval_tree(sf_rmsk, i_min_copy_len, true, 100)?;
        let xtprt_filter = XTPRTFilter::default();
        let l_types = AFConflictFilter::get_rep_type();
        let m_cutoff = Self::get_af_cutoff(&l_types);
        let sf_before = format!("{}.before_blacklist.txt", sf_new_out);
        let rep = if i_rep_type & 1 != 0 {
            Some(MosaicRep::L1)
        } else if i_rep_type & 2 != 0 {
            Some(MosaicRep::Alu)
        } else if i_rep_type & 4 != 0 {
            Some(MosaicRep::Sva)
        } else {
            None
        };
        if let Some(rep) = rep {
            let mut fout = utils::create(&sf_before)?;
            for rcd in l_old_rcd.iter() {
                let s_type = XTEARsltParser::get_ins_sub_type(rcd)?;
                let ins_chrm = rcd.first().ok_or_else(|| anyhow!("IndexError: list index out of range"))?;
                let ins_pos = py_int(rcd.get(1).ok_or_else(|| anyhow!("IndexError: list index out of range"))?)?;
                let (b_in_rep, _i_pos) = xannotation.is_within_repeat_region_interval_tree(ins_chrm, ins_pos);
                if !xtprt_filter.has_polya_signal(rcd)? {
                    continue;
                }
                let s_line = rcd.last().ok_or_else(|| anyhow!("IndexError: list index out of range"))?;
                if !AFConflictFilter::is_qualified_mosaic_rcd(s_line, &m_cutoff)? {
                    continue;
                }
                if s_type == TWO_SIDE {
                    if xtprt_filter.is_polya_dominant_two_side(rcd)? {
                        continue;
                    }
                } else if s_type == ONE_HALF_SIDE {
                    if b_in_rep {
                        continue;
                    }
                    // the SVA loop has this check commented out
                    if rep != MosaicRep::Sva && xtprt_filter.is_polya_dominant_two_side(rcd)? {
                        continue;
                    }
                } else if s_type == ONE_SIDE {
                    if b_in_rep {
                        continue;
                    }
                } else {
                    continue;
                }
                let name = match rep {
                    MosaicRep::L1 => "LINE1",
                    MosaicRep::Alu => "Alu",
                    MosaicRep::Sva => "SVA",
                };
                let s_in_rep = if b_in_rep { format!("Fall_in_{}_copy", name) } else { format!("not_in_{}_copy", name) };
                writeln!(fout, "{}\t{}", utils::rstrip(s_line), s_in_rep)?;
            }
            fout.flush()?;
        }
        self.filter_by_blacklist(&sf_before, sf_black_list, self.brkpnt_slack, sf_new_out)
    }

    /// `filter_by_blacklist`
    pub fn filter_by_blacklist(&self, sf_in: &str, sf_black_list: &str, i_slack: i64, sf_out: &str) -> Result<()> {
        let m_meis = SomaticMeiCaller::load_sites(sf_in)?;
        let m_pass = SomaticMeiCaller::filter_by_given_list(sf_black_list, &m_meis, i_slack)?;
        SomaticMeiCaller::export_slcted_candidates(&m_pass, sf_out)
    }
}

/// `MosaicJointCalling`: keep the calls seen in a single sample.
pub struct MosaicJointCalling {
    pub swfolder: String,
    pub n_jobs: usize,
}

impl MosaicJointCalling {
    pub fn new(s_working_folder: &str, n_jobs: usize) -> Result<Self> {
        let mut wf = s_working_folder.to_string();
        if !wf.ends_with('/') {
            if wf.is_empty() {
                return Err(anyhow!("IndexError: string index out of range"));
            }
            wf.push('/');
        }
        Ok(MosaicJointCalling { swfolder: wf, n_jobs })
    }

    /// `_is_within_nearby_region`
    fn is_within_nearby_region<T>(m: &SiteMap<T>, chrm: &str, i_pos: i64, i_slack: i64) -> (bool, i64) {
        let mm = match m.get(chrm) {
            None => return (false, -1),
            Some(mm) => mm,
        };
        for i in (i_pos - i_slack)..(i_pos + i_slack) {
            if mm.contains_key(&i) {
                return (true, i);
            }
        }
        (false, -1)
    }

    /// `_load_rslts` (the last line of a position wins)
    fn load_rslts(sf: &str) -> Result<SiteMap<String>> {
        let mut m: SiteMap<String> = IndexMap::new();
        for line in utils::read_raw_lines(sf)? {
            let fields = split_ws(&line);
            let chrm = fields.first().ok_or_else(|| anyhow!("IndexError: list index out of range"))?;
            let pos = py_int(fields.get(1).ok_or_else(|| anyhow!("IndexError: list index out of range"))?)?;
            m.entry(chrm.to_string()).or_default().insert(pos, line.clone());
        }
        Ok(m)
    }

    /// `_get_diff_set`: in the original the "not within the selected set" test compares a
    /// tuple with False and is never true, so only chromosomes absent from the selected set
    /// contribute.
    fn get_diff_set(sf_all: &str, sf_slct: &str, m_diff: &mut SiteMap<()>) -> Result<()> {
        let m_all = Self::load_rslts(sf_all)?;
        let m_slct = Self::load_rslts(sf_slct)?;
        for (chrm, m) in m_all.iter() {
            for &pos in m.keys() {
                if !m_slct.contains_key(chrm) {
                    m_diff.entry(chrm.clone()).or_default().insert(pos, ());
                }
            }
        }
        Ok(())
    }

    /// `call_mosaic_from_multi_samples`
    pub fn call_mosaic_from_multi_samples(&self, sf_rslt_list: &str, islack: i64, sf_out: &str) -> Result<()> {
        let mut l_rslts: Vec<(String, String, String)> = Vec::new();
        for line in utils::read_raw_lines(sf_rslt_list)? {
            let fields = split_ws(&line);
            let g = |i: usize| -> Result<String> {
                fields.get(i).map(|s| s.to_string()).ok_or_else(|| anyhow!("IndexError: list index out of range"))
            };
            let (sid, sf_rslt, sf_rslt_all) = (g(0)?, g(1)?, g(2)?);
            if !utils::is_file(&sf_rslt) || !utils::is_file(&sf_rslt_all) {
                continue;
            }
            l_rslts.push((sid, sf_rslt, sf_rslt_all));
        }
        let mut m_candidates: SiteMap<Vec<(String, String)>> = IndexMap::new();
        let mut m_not_slcted: SiteMap<()> = IndexMap::new();
        for (sid, sf_rslt, sf_rslt_all) in l_rslts.iter() {
            Self::get_diff_set(sf_rslt_all, sf_rslt, &mut m_not_slcted)?;
            for line in utils::read_raw_lines(sf_rslt)? {
                let fields = split_ws(&line);
                if fields.len() < 2 {
                    continue;
                }
                let chrm = fields[0];
                let pos = py_int(fields[1])?;
                if !m_candidates.contains_key(chrm) {
                    let mut m = IndexMap::new();
                    m.insert(pos, vec![(sid.clone(), line.clone())]);
                    m_candidates.insert(chrm.to_string(), m);
                } else {
                    let (b_hit, i_hit_pos) = Self::is_within_nearby_region(&m_candidates, chrm, pos, islack);
                    let m = m_candidates.get_mut(chrm).unwrap();
                    if b_hit {
                        m.get_mut(&i_hit_pos).unwrap().push((sid.clone(), line.clone()));
                    } else {
                        m.insert(pos, vec![(sid.clone(), line.clone())]);
                    }
                }
            }
        }
        // `(chrm in m_not_slcted) and _is_within_nearby_region(...)` is a non-empty tuple, i.e.
        // true whenever the chromosome has any not-selected site
        {
            let mut fout_singleton = utils::create(&format!("{}.singleton_for_screenshot", sf_out))?;
            let mut fout_non = utils::create(&format!("{}.non_singleton_for_screenshot", sf_out))?;
            for (chrm, m) in m_candidates.iter() {
                for v in m.values() {
                    if v.len() == 1 {
                        if m_not_slcted.contains_key(chrm) {
                            continue;
                        }
                        writeln!(fout_singleton, "{}\t{}", v[0].0, utils::rstrip(&v[0].1))?;
                    } else {
                        let s_tmp_line = utils::rstrip(&v[0].1);
                        for rcd in v.iter() {
                            writeln!(fout_non, "{}\t{}", rcd.0, s_tmp_line)?;
                        }
                    }
                }
            }
            fout_singleton.flush()?;
            fout_non.flush()?;
        }
        {
            let mut fout_singleton2 = utils::create(&format!("{}.singleton2", sf_out))?;
            let mut fout_non2 = utils::create(&format!("{}.non_singleton2", sf_out))?;
            for (chrm, m) in m_candidates.iter() {
                for v in m.values() {
                    if v.len() == 1 {
                        if m_not_slcted.contains_key(chrm) {
                            continue;
                        }
                        writeln!(fout_singleton2, "{}\t{}", utils::rstrip(&v[0].1), v[0].0)?;
                    } else {
                        let mut s_tmp_line = utils::rstrip(&v[0].1).to_string();
                        for rcd in v.iter() {
                            s_tmp_line.push('\t');
                            s_tmp_line.push_str(&rcd.0);
                        }
                        writeln!(fout_non2, "{}", s_tmp_line)?;
                    }
                }
            }
            fout_singleton2.flush()?;
            fout_non2.flush()?;
        }
        Ok(())
    }
}
