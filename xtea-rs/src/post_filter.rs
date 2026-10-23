//! Port of `x_post_filter.py`: `XPostFilter`, `XTEARsltParser`, `XTPRTFilter` and
//! `AFConflictFilter`.
//!
//! The records are kept as the whitespace separated fields of the result lines (strings), plus
//! - for `load_in_xtea_rslt` - the original line as last element, exactly like the Python lists.

use crate::annotation::XAnnotation;
use crate::black_list::XBlackList;
use crate::globals as gv;
use crate::num::PyNum;
use crate::utils::{self, split_ws};
use anyhow::{anyhow, Result};
use indexmap::IndexMap;
use std::collections::HashMap;
use std::io::Write;

/// Python `int(float(s))`
pub fn int_float(s: &str) -> Result<i64> {
    let f = crate::py::fmt::parse_py_float(s).ok_or_else(|| anyhow!("ValueError: could not convert string to float: '{}'", s))?;
    Ok(f.trunc() as i64)
}

fn pint(s: &str) -> Result<i64> {
    utils::pint(s)
}

fn pfloat(s: &str) -> Result<f64> {
    crate::py::fmt::parse_py_float(s).ok_or_else(|| anyhow!("ValueError: could not convert string to float: '{}'", s))
}

fn fld<'a>(rcd: &'a [String], i: usize) -> Result<&'a str> {
    rcd.get(i).map(|s| s.as_str()).ok_or_else(|| anyhow!("IndexError: list index {} out of range", i))
}

/// (start, end) of a "start:end" cluster field, with `int(float())` conversion.
fn cluster(s: &str) -> Result<(i64, i64)> {
    let f: Vec<&str> = s.split(':').collect();
    if f.len() < 2 {
        return Err(anyhow!("IndexError: malformed cluster field {}", s));
    }
    Ok((int_float(f[0])?, int_float(f[1])?))
}

// ------------------------------------------------------------------------------------------------
// XTEARsltParser

pub const TWO_SIDE: &str = "two_side";
pub const ONE_SIDE: &str = "one_side";
pub const ONE_HALF_SIDE: &str = "one_half";
pub const OTHER: &str = "other";

pub struct XTEARsltParser;

impl XTEARsltParser {
    /// `load_in_xTEA_rslt`: fields + the raw line (with its newline) as last element
    pub fn load_in_xtea_rslt(sf_rslt: &str) -> Result<Vec<Vec<String>>> {
        let mut l = Vec::new();
        for line in utils::read_raw_lines(sf_rslt)? {
            let mut rcd: Vec<String> = split_ws(utils::rstrip(&line)).into_iter().map(|s| s.to_string()).collect();
            rcd.push(line);
            l.push(rcd);
        }
        Ok(l)
    }

    /// `load_in_xTEA_rslt_ori`: fields only
    pub fn load_in_xtea_rslt_ori(sf_rslt: &str) -> Result<Vec<Vec<String>>> {
        let mut l = Vec::new();
        for line in utils::read_raw_lines(sf_rslt)? {
            l.push(split_ws(utils::rstrip(&line)).into_iter().map(|s| s.to_string()).collect());
        }
        Ok(l)
    }

    /// `dump_xTEA_rslt_to_file`
    pub fn dump_xtea_rslt_to_file(l_rcd: &[Vec<String>], sf_rslt: &str) -> Result<()> {
        let mut f = utils::create(sf_rslt)?;
        for r in l_rcd {
            if let Some(last) = r.last() {
                f.write_all(last.as_bytes())?;
            }
        }
        Ok(())
    }

    fn sub_type_of(s_type: &str) -> &'static str {
        if s_type.contains("two_side") || s_type.contains("both-side") || s_type.contains("one_side_and_half_transduction") {
            TWO_SIDE
        } else if s_type.contains("one_side") {
            ONE_SIDE
        } else if s_type.contains("one_half") || s_type.contains("one-half") {
            ONE_HALF_SIDE
        } else {
            OTHER
        }
    }

    /// `get_ins_sub_type`
    pub fn get_ins_sub_type(rcd: &[String]) -> Result<&'static str> {
        Ok(Self::sub_type_of(fld(rcd, 32)?))
    }

    /// `get_ins_sub_type_alu`
    pub fn get_ins_sub_type_alu(rcd: &[String]) -> Result<&'static str> {
        let n_ef_lclip = pint(fld(rcd, 5)?)?;
        let n_ef_rclip = pint(fld(rcd, 6)?)?;
        let t = Self::sub_type_of(fld(rcd, 32)?);
        if t == ONE_HALF_SIDE && n_ef_lclip > 0 && n_ef_rclip > 0 {
            return Ok(TWO_SIDE);
        }
        Ok(t)
    }

    /// `get_ins_sub_type_sva`
    pub fn get_ins_sub_type_sva(rcd: &[String]) -> Result<&'static str> {
        Self::get_ins_sub_type(rcd)
    }

    /// `replace_ins_length`: drop the last field (and the line), append the re-estimated
    /// length and the new line.
    pub fn replace_ins_length(old_rcd: &[String], icns_lth: i64) -> Result<Vec<String>> {
        let n = old_rcd.len().saturating_sub(2);
        let mut rcd: Vec<String> = old_rcd[..n].to_vec();
        let new_len = Self::estimate_insertion_length(&rcd, icns_lth)?;
        rcd.push(new_len.to_string());
        let s_rcd = rcd.join("\t");
        rcd.push(s_rcd);
        Ok(rcd)
    }

    /// `estimate_insertion_length`
    pub fn estimate_insertion_length(rcd: &[String], icns_lth: i64) -> Result<i64> {
        let mut i_tei_len: i64 = 0;
        let (i_lclip_start, i_lclip_end) = cluster(fld(rcd, 19)?)?;
        let (i_rclip_start, i_rclip_end) = cluster(fld(rcd, 20)?)?;
        let (i_ldisc_start, i_ldisc_end) = cluster(fld(rcd, 21)?)?;
        let (i_rdisc_start, i_rdisc_end) = cluster(fld(rcd, 22)?)?;
        if !Self::is_transduction(rcd)? {
            if i_lclip_start > 0 && i_rclip_start > 0 {
                i_tei_len = i_rclip_end - i_lclip_start;
                if i_tei_len < 0 {
                    i_tei_len = i_lclip_end - i_rclip_start;
                }
            } else if i_ldisc_start > 0 && i_rdisc_start > 0 {
                i_tei_len = i_rdisc_end - i_ldisc_start;
                if i_tei_len < 0 {
                    i_tei_len = i_ldisc_end - i_rdisc_start;
                }
            } else if i_lclip_start > 0 && i_ldisc_start > 0 {
                i_tei_len = i_ldisc_end - i_lclip_start;
                if i_tei_len < 0 {
                    i_tei_len = i_lclip_end - i_ldisc_start;
                }
            } else if i_rclip_start > 0 && i_rdisc_start > 0 {
                i_tei_len = i_rdisc_end - i_rclip_start;
                if i_tei_len < 0 {
                    i_tei_len = i_rclip_end - i_rdisc_start;
                }
            }
        } else if i_rdisc_start > 0 && i_ldisc_start > 0 {
            if i_rdisc_start < i_ldisc_start {
                i_tei_len = icns_lth - i_rdisc_start;
            } else {
                i_tei_len = icns_lth - i_ldisc_start;
            }
        } else if i_rdisc_start < 0 && i_ldisc_start > 0 {
            i_tei_len = icns_lth - i_ldisc_start;
        } else if i_rdisc_start > 0 && i_ldisc_start < 0 {
            i_tei_len = icns_lth - i_rdisc_start;
        } else {
            i_tei_len = 0;
        }
        Ok(i_tei_len.abs())
    }

    /// `is_sibling_transduction_related_ins`
    pub fn is_sibling_transduction_related_ins(rcd: &[String]) -> Result<bool> {
        let s = fld(rcd, 32)?;
        Ok(s.contains("one_side") || s.contains("one_half"))
    }

    /// `get_two_side_coverage_diff`
    pub fn get_two_side_coverage_diff(rcd: &[String]) -> Result<f64> {
        Ok((pfloat(fld(rcd, 11)?)? - pfloat(fld(rcd, 12)?)?).abs())
    }

    /// `is_transduction`
    pub fn is_transduction(rcd: &[String]) -> Result<bool> {
        let s_src = fld(rcd, 23)?;
        Ok(!(s_src == "not_transduction" || s_src.contains("may")))
    }

    /// `is_orphan_transduction`
    pub fn is_orphan_transduction(rcd: &[String]) -> Result<bool> {
        Ok(fld(rcd, 32)?.contains("orphan"))
    }

    /// `get_transduction_source`
    pub fn get_transduction_source(rcd: &[String]) -> Result<String> {
        Ok(fld(rcd, 23)?.to_string())
    }

    fn set(rcd: &mut [String], i: usize, v: String) -> Result<()> {
        match rcd.get_mut(i) {
            Some(x) => {
                *x = v;
                Ok(())
            }
            None => Err(anyhow!("IndexError: list assignment index {} out of range", i)),
        }
    }

    /// `update_td_info`
    pub fn update_td_info(rcd: &mut [String], s_td_info: &str, n_clip: i64, n_disc: i64, n_polya: i64) -> Result<()> {
        Self::set(rcd, 23, s_td_info.to_string())?;
        Self::set(rcd, 24, n_clip.to_string())?;
        Self::set(rcd, 25, n_disc.to_string())?;
        Self::set(rcd, 26, n_polya.to_string())
    }

    /// `update_td_source_only`
    pub fn update_td_source_only(rcd: &mut [String], s_td_info: &str) -> Result<()> {
        Self::set(rcd, 23, s_td_info.to_string())
    }
}

// ------------------------------------------------------------------------------------------------
// XTPRTFilter

pub struct XTPRTFilter {
    pub f_side_polya_cutoff: f64,
}

impl Default for XTPRTFilter {
    fn default() -> Self {
        XTPRTFilter { f_side_polya_cutoff: gv::ONE_SIDE_POLYA_CUTOFF }
    }
}

impl XTPRTFilter {
    /// `is_polyA_dominant_two_side`
    pub fn is_polya_dominant_two_side(&self, rcd: &[String]) -> Result<bool> {
        let n_lpolya = pint(fld(rcd, 9)?)?;
        let n_rpolya = pint(fld(rcd, 10)?)?;
        let n_ef_lclip = pint(fld(rcd, 5)?)?;
        let n_ef_rclip = pint(fld(rcd, 6)?)?;
        let s_type = fld(rcd, 32)?;
        if n_ef_lclip <= 0 || n_ef_rclip <= 0 {
            return Ok(false);
        }
        if self.is_two_clip_form_different_cluster(rcd)? {
            return Ok(false);
        }
        if s_type.contains(gv::TWO_SIDE_TPRT_BOTH) && self.is_two_disc_form_different_cluster(rcd)? {
            return Ok(false);
        }
        let b_l = (n_lpolya as f64) / (n_ef_lclip as f64) > self.f_side_polya_cutoff;
        let b_r = (n_rpolya as f64) / (n_ef_rclip as f64) > self.f_side_polya_cutoff;
        Ok(b_l && b_r)
    }

    /// `is_two_side_polyA_same_orientation`
    pub fn is_two_side_polya_same_orientation(&self, rcd: &[String]) -> Result<bool> {
        let a = pint(fld(rcd, 17)?)?;
        let b = pint(fld(rcd, 18)?)?;
        Ok(a * b > 0)
    }

    /// `is_polyA_dominant_td`
    pub fn is_polya_dominant_td(&self, rcd: &[String]) -> Result<bool> {
        let f_cutoff = gv::MAX_POLYA_RATIO;
        let n_lpolya = pint(fld(rcd, 9)?)?;
        let n_rpolya = pint(fld(rcd, 10)?)?;
        let n_ef_lclip = pint(fld(rcd, 5)?)?;
        let n_ef_rclip = pint(fld(rcd, 6)?)?;
        let s_type = fld(rcd, 32)?;
        if s_type.contains(gv::TWO_SIDE_TPRT_BOTH) {
            return Ok(false);
        }
        let mut b_l = true;
        if n_ef_lclip != 0 {
            b_l = (n_lpolya as f64) / (n_ef_lclip as f64) > f_cutoff;
        }
        let mut b_r = true;
        if n_ef_rclip != 0 {
            b_r = (n_rpolya as f64) / (n_ef_rclip as f64) > f_cutoff;
        }
        Ok(b_l && b_r)
    }

    /// `is_two_side_consist_with_enough_ins_size`
    pub fn is_two_side_consist_with_enough_ins_size(&self, rcd: &[String], min_ins_len: i64) -> Result<bool> {
        let s_consist = fld(rcd, 33)?;
        let (ls, le) = cluster(fld(rcd, 19)?)?;
        let (rs, re) = cluster(fld(rcd, 20)?)?;
        Ok(s_consist == gv::BOTH_END_CONSISTNT && ((ls - re).abs() > min_ins_len || (rs - le).abs() > min_ins_len))
    }

    /// `is_two_side_tprt_both_and_both_consistent`
    pub fn is_two_side_tprt_both_and_both_consistent(&self, rcd: &[String], i_cns_tail: i64, min_ins_len: i64) -> Result<bool> {
        let s_type = fld(rcd, 32)?;
        let s_consist = fld(rcd, 33)?;
        let (ls, le) = cluster(fld(rcd, 19)?)?;
        let (rs, re) = cluster(fld(rcd, 20)?)?;
        let (_lds, lde) = cluster(fld(rcd, 21)?)?;
        let (_rds, rde) = cluster(fld(rcd, 22)?)?;
        if (le < i_cns_tail && re < i_cns_tail) || (lde < i_cns_tail && rde < i_cns_tail) {
            return Ok(false);
        }
        if (ls - re).abs() < min_ins_len && (rs - le).abs() < min_ins_len {
            return Ok(false);
        }
        Ok(s_type.contains(gv::TWO_SIDE_TPRT_BOTH) && s_consist == gv::BOTH_END_CONSISTNT)
    }

    /// `is_polyA_dominant_two_side_sva`
    pub fn is_polya_dominant_two_side_sva(&self, rcd: &[String], i_pos_head: i64) -> Result<bool> {
        let n_lpolya = pint(fld(rcd, 9)?)?;
        let n_rpolya = pint(fld(rcd, 10)?)?;
        let n_ef_lclip = pint(fld(rcd, 5)?)?;
        let n_ef_rclip = pint(fld(rcd, 6)?)?;
        let s_type = fld(rcd, 32)?;
        if n_ef_lclip <= 0 || n_ef_rclip <= 0 {
            return Ok(false);
        }
        if s_type.contains(gv::TWO_SIDE_TPRT_BOTH) && self.disc_cluster_hit_cns_head(rcd, i_pos_head)? {
            return Ok(false);
        }
        let b_l = (n_lpolya as f64) / (n_ef_lclip as f64) > self.f_side_polya_cutoff;
        let b_r = (n_rpolya as f64) / (n_ef_rclip as f64) > self.f_side_polya_cutoff;
        Ok(b_l && b_r)
    }

    fn clip_clusters(rcd: &[String]) -> Result<Option<((i64, i64), (i64, i64))>> {
        let l = fld(rcd, 19)?;
        let r = fld(rcd, 20)?;
        if l.contains("-1") || r.contains("-1") {
            return Ok(None);
        }
        Ok(Some((cluster(l)?, cluster(r)?)))
    }

    /// `is_two_side_clip_both_polyA_sva`
    pub fn is_two_side_clip_both_polya_sva(&self, rcd: &[String], i_cns_tail: i64) -> Result<bool> {
        Ok(match Self::clip_clusters(rcd)? {
            None => false,
            Some(((ls, _), (rs, _))) => ls >= i_cns_tail && rs >= i_cns_tail,
        })
    }

    /// `is_two_side_clip_both_non_polyA_sva`
    pub fn is_two_side_clip_both_non_polya_sva(&self, rcd: &[String], i_cns_tail: i64) -> Result<bool> {
        Ok(match Self::clip_clusters(rcd)? {
            None => false,
            Some(((ls, _), (rs, _))) => ls < i_cns_tail && rs < i_cns_tail,
        })
    }

    /// `_is_two_clip_form_different_cluster`
    fn is_two_clip_form_different_cluster(&self, rcd: &[String]) -> Result<bool> {
        Ok(match Self::clip_clusters(rcd)? {
            None => true,
            Some(((ls, le), (rs, re))) => {
                let c = gv::two_clip_cluster_diff_cutoff();
                !((le - rs).abs() < c || (re - ls).abs() < c)
            }
        })
    }

    /// `_is_two_disc_form_different_cluster`
    fn is_two_disc_form_different_cluster(&self, rcd: &[String]) -> Result<bool> {
        let l = fld(rcd, 21)?;
        let r = fld(rcd, 22)?;
        if l.contains("-1") || r.contains("-1") {
            return Ok(true);
        }
        let (ls, le) = cluster(l)?;
        let (rs, re) = cluster(r)?;
        let c = gv::two_clip_cluster_diff_cutoff();
        Ok(!((le - rs).abs() < c || (re - ls).abs() < c))
    }

    /// `_disc_cluster_hit_cns_head`
    fn disc_cluster_hit_cns_head(&self, rcd: &[String], i_pos_head: i64) -> Result<bool> {
        let l = fld(rcd, 21)?;
        let r = fld(rcd, 22)?;
        if l.contains("-1") || r.contains("-1") {
            return Ok(false);
        }
        let (ls, _) = cluster(l)?;
        let (rs, _) = cluster(r)?;
        Ok(ls < i_pos_head || rs < i_pos_head)
    }

    /// `has_polyA_signal`
    pub fn has_polya_signal(&self, rcd: &[String]) -> Result<bool> {
        let l = pint(fld(rcd, 9)?)?;
        let r = pint(fld(rcd, 10)?)?;
        Ok(!(l <= 0 && r <= 0))
    }

    /// `is_polyA_dominant_one_side`
    pub fn is_polya_dominant_one_side(&self, rcd: &[String], nclip_half_cutoff: i64) -> Result<bool> {
        let f_cutoff = self.f_side_polya_cutoff;
        let n_lpolya = pint(fld(rcd, 9)?)?;
        let n_rpolya = pint(fld(rcd, 10)?)?;
        let n_ef_lclip = pint(fld(rcd, 5)?)?;
        let n_ef_rclip = pint(fld(rcd, 6)?)?;
        let n_lr_clip = pint(fld(rcd, 35)?)?;
        if n_ef_lclip <= 0 && n_ef_rclip <= 0 {
            return Ok(true);
        }
        if n_ef_lclip <= 0 {
            let b_r = (n_rpolya as f64) / (n_ef_rclip as f64) > f_cutoff;
            if b_r && (n_lr_clip - n_ef_rclip) >= nclip_half_cutoff {
                return Ok(true);
            }
        }
        if n_ef_rclip <= 0 {
            let b_l = (n_lpolya as f64) / (n_ef_lclip as f64) > f_cutoff;
            if b_l && (n_lr_clip - n_ef_lclip) >= nclip_half_cutoff {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// `is_two_side_tprt_and_with_polyA`
    pub fn is_two_side_tprt_and_with_polya(&self, rcd: &[String], i_cns_end: i64) -> Result<bool> {
        let s_type = fld(rcd, 32)?;
        if s_type.contains(gv::TWO_SIDE_TPRT) {
            if let Some(((_, le), (_, re))) = Self::clip_clusters(rcd)? {
                if le > i_cns_end || re > i_cns_end {
                    return Ok(true);
                }
            } else {
                return Ok(false);
            }
        }
        Ok(false)
    }

    /// `is_two_side_clip_and_both_hit_end`
    pub fn is_two_side_clip_and_both_hit_end(&self, rcd: &[String], i_cns_end: i64) -> Result<bool> {
        Ok(match Self::clip_clusters(rcd)? {
            None => false,
            Some(((ls, _), (rs, _))) => ls > i_cns_end && rs > i_cns_end,
        })
    }

    /// `hit_consensus_tail`
    pub fn hit_consensus_tail(&self, rcd: &[String], i_cns_end: i64) -> Result<bool> {
        Ok(match Self::clip_clusters(rcd)? {
            None => true,
            Some(((_, le), (_, re))) => le >= i_cns_end || re >= i_cns_end,
        })
    }

    /// `hit_end_of_cns`
    pub fn hit_end_of_cns(&self, rcd: &[String]) -> Result<bool> {
        Ok(fld(rcd, 34)? == gv::HIT_END_OF_CNS)
    }

    /// `has_enough_indel_reads`
    pub fn has_enough_indel_reads(&self, rcd: &[String], ncutoff: i64) -> Result<bool> {
        Ok(pint(fld(rcd, 41)?)? > ncutoff)
    }

    /// `cov_is_abnormal`
    pub fn cov_is_abnormal(&self, rcd: &[String], f_cutoff: f64) -> Result<bool> {
        let l = pfloat(fld(rcd, 11)?)?;
        let r = pfloat(fld(rcd, 12)?)?;
        Ok(l > f_cutoff || r > f_cutoff)
    }
}

// ------------------------------------------------------------------------------------------------
// AFConflictFilter

pub type AfCutoff = HashMap<String, (f64, f64, f64, f64)>;

pub struct AFConflictFilter;

impl AFConflictFilter {
    /// `get_rep_type`
    pub fn get_rep_type() -> Vec<&'static str> {
        vec![
            gv::ONE_SIDE_FLANKING,
            gv::TWO_SIDE,
            gv::TWO_SIDE_TPRT_BOTH,
            gv::TWO_SIDE_TPRT,
            gv::ONE_HALF_SIDE,
            gv::ONE_HALF_SIDE_TRPT_BOTH,
            gv::ONE_HALF_SIDE_TRPT,
            gv::ONE_HALF_SIDE_POLYA_DOMINANT,
            gv::ONE_SIDE,
            gv::ONE_SIDE_COVERAGE_CONFLICT,
            gv::ONE_SIDE_TRSDCT,
            gv::ONE_SIDE_WEAK,
            gv::ONE_SIDE_OTHER,
            gv::ONE_SIDE_SV,
            gv::ONE_SIDE_POLYA_DOMINANT,
            gv::TWO_SIDE_POLYA_DOMINANT,
            gv::HIGH_COV_ISD,
            gv::OTHER_TYPE,
        ]
    }

    /// `get_cutoff_by_type` (all categories use the same values in the original)
    pub fn get_cutoff_by_type(l_types: &[&str], b_tumor: bool) -> AfCutoff {
        let mut m = HashMap::new();
        for t in l_types {
            let v = if b_tumor { (0.001, 0.001, 0.001, 0.001) } else { (0.075, 0.075, 0.075, 0.075) };
            m.insert(t.to_string(), v);
        }
        m
    }

    /// `is_qualified_td`
    pub fn is_qualified_td(rcd: &[String]) -> Result<bool> {
        let n_disc = pint(fld(rcd, 39)?)?;
        let n_concod = pint(fld(rcd, 40)?)?;
        let mut f = 0.0;
        if n_disc + n_concod != 0 {
            f = n_disc as f64 / (n_disc + n_concod) as f64;
        }
        Ok(!(f < 0.08))
    }

    fn ratios(s_line: &str) -> Result<(f64, f64, f64, f64, String)> {
        let fields: Vec<String> = split_ws(s_line).into_iter().map(|s| s.to_string()).collect();
        let _ = pint(fld(&fields, 9)?)?;
        let _ = pint(fld(&fields, 10)?)?;
        let n_ef_clip = pint(fld(&fields, 5)?)? + pint(fld(&fields, 6)?)?;
        let n_ef_disc = pint(fld(&fields, 7)?)? + pint(fld(&fields, 8)?)?;
        let n_clip = pint(fld(&fields, 35)?)?;
        let n_full_map = pint(fld(&fields, 36)?)?;
        let n_disc = pint(fld(&fields, 39)?)?;
        let n_concod = pint(fld(&fields, 40)?)?;
        let mut f_ef_clip = 0.0;
        if n_clip != 0 {
            f_ef_clip = n_ef_clip as f64 / n_clip as f64;
        }
        let mut f_ef_disc = 0.0;
        if n_disc != 0 {
            f_ef_disc = n_ef_disc as f64 / n_disc as f64;
        }
        let mut f_clip_full_map = 0.0;
        if n_clip + n_full_map != 0 {
            f_clip_full_map = n_clip as f64 / (n_clip + n_full_map) as f64;
        }
        let mut f_disc_concod = 0.0;
        if n_disc + n_concod != 0 {
            f_disc_concod = n_disc as f64 / (n_disc + n_concod) as f64;
        }
        Ok((f_ef_clip, f_ef_disc, f_clip_full_map, f_disc_concod, fld(&fields, 32)?.to_string()))
    }

    /// `is_qualified_rcd`
    pub fn is_qualified_rcd(s_line: &str, m_cutoff: &AfCutoff) -> Result<bool> {
        let (a, b, c, d, s_type) = Self::ratios(s_line)?;
        Self::is_ins_pass_cutoff(m_cutoff, &s_type, a, b, c, d)
    }

    /// `is_ins_pass_cutoff`
    pub fn is_ins_pass_cutoff(m_cutoff: &AfCutoff, s_type: &str, f_ef_clip: f64, f_ef_disc: f64, f_clip_full_map: f64, f_disc_concod: f64) -> Result<bool> {
        let (c1, c2, c3, c4) = *m_cutoff.get(s_type).ok_or_else(|| anyhow!("KeyError: '{}'", s_type))?;
        Ok(f_ef_clip > c1 && f_ef_disc > c2 && f_clip_full_map > c3 && f_disc_concod > c4)
    }

    /// `is_qualified_mosaic_rcd` (the cutoff map holds (upper, lower) allele fractions)
    pub fn is_qualified_mosaic_rcd(s_line: &str, m_cutoff: &HashMap<String, (f64, f64)>) -> Result<bool> {
        let (_a, _b, c, d, s_type) = Self::ratios(s_line)?;
        let (f_upper_af, f_lower_af) = *m_cutoff.get(&s_type).ok_or_else(|| anyhow!("KeyError: '{}'", s_type))?;
        let b_clip = c >= f_lower_af && c <= f_upper_af;
        let b_disc = !(d > f_upper_af);
        Ok(b_clip && b_disc)
    }
}

// ------------------------------------------------------------------------------------------------
// XPostFilter

pub struct XPostFilter {
    pub swfolder: String,
    pub n_jobs: usize,
}

const REP_ALU_MIN_LEN: i64 = 200;
const REP_ALU_POLYA_START: i64 = 255;
const REP_LINE_POLYA_START: i64 = 5950;
const REP_SVA_CNS_HEAD: i64 = 400;
const REP_SVA_POLYA_START: i64 = 1900;
const NCLIP_HALF_CUTOFF: i64 = 5;
const SVA_CLIP_CLUSTER_DIFF_CUTOFF: i64 = 200;
const INDEL_READS_MAX_RATIO: f64 = 0.3;
const ABNORMAL_COV_TIMES: f64 = 2.0;
const L1_BOUNDARY_EXTND: i64 = 400;
const L1_MIN_REF_COPY_LEN: i64 = 125;
const BOUNDARY_EXTND: i64 = 100;

/// Python `str(div_rate)` where div_rate is the parsed float of the RepeatMasker column (or
/// the int -1 fallback).
fn div_str(d: PyNum) -> String {
    d.py_str()
}

impl XPostFilter {
    pub fn new(swfolder: &str, n_jobs: usize) -> Self {
        XPostFilter { swfolder: utils::with_slash(swfolder), n_jobs }
    }

    /// `construct_interval_tree`
    pub fn construct_interval_tree(sf_rmsk: &str, i_min_copy_len: i64, b_with_chr: bool, boundary_extnd: i64) -> Result<XAnnotation> {
        let mut xannotation = XAnnotation::new(sf_rmsk);
        xannotation.set_with_chr(b_with_chr);
        xannotation.load_rmsk_annotation_with_extnd_div_with_lenth_cutoff(boundary_extnd, i_min_copy_len)?;
        xannotation.index_rmsk_annotation_interval_tree();
        Ok(xannotation)
    }

    /// `fall_in_low_div_same_type_rep` -> (hit, lowest div (int -1 when none), copy_start, copy_end)
    pub fn fall_in_low_div_same_type_rep(xannotation: &XAnnotation, ins_chrm: &str, ins_pos: i64, i_cutoff: i64) -> (bool, PyNum, i64, i64) {
        let l_hits = xannotation.is_within_repeat_region_interval_tree2(ins_chrm, ins_pos);
        if l_hits.is_empty() {
            return (false, PyNum::Int(-1), -1, -1);
        }
        let mut f_lowest_div = PyNum::Int(-1);
        let mut b_hit = false;
        let mut copy_start = -1i64;
        let mut copy_end = -1i64;
        for (b_in_rep, i_pos) in l_hits {
            for (div_rate, _sub, _fam, pos_start, pos_end) in xannotation.get_div_subfamily_with_min_div2(ins_chrm, i_pos) {
                if b_in_rep && (f_lowest_div.f() == -1.0 || f_lowest_div.f() > div_rate) {
                    f_lowest_div = PyNum::Float(div_rate);
                    b_hit = true;
                    copy_start = pos_start;
                    copy_end = pos_end;
                }
            }
        }
        if b_hit && f_lowest_div.f() <= i_cutoff as f64 && f_lowest_div.f() != -1.0 {
            return (true, f_lowest_div, copy_start, copy_end);
        }
        (false, f_lowest_div, copy_start, copy_end)
    }

    /// `run_post_filtering`
    #[allow(clippy::too_many_arguments)]
    pub fn run_post_filtering(&self, sf_xtea_rslt: &str, sf_rmsk: &str, i_min_copy_len: i64, i_rep_type: i64, f_cov: f64, sf_black_list: &str, sf_new_out: &str, b_tumor: bool) -> Result<()> {
        let l_old_rcd = XTEARsltParser::load_in_xtea_rslt(sf_xtea_rslt)?;
        let xtprt = XTPRTFilter::default();
        let l_types = AFConflictFilter::get_rep_type();
        let m_cutoff = AFConflictFilter::get_cutoff_by_type(&l_types, b_tumor);
        let mut x_blklist = XBlackList::new();
        x_blklist.load_index_regions(sf_black_list)?;
        if i_rep_type & 1 != 0 {
            let sf_td_new_sites = format!("{}{}{}", sf_xtea_rslt, gv::TD_NON_SIBLING_SUFFIX, gv::TD_NEW_SITES_SUFFIX);
            let xannotation = Self::construct_interval_tree(sf_rmsk, L1_MIN_REF_COPY_LEN, true, L1_BOUNDARY_EXTND)?;
            self.post_processing_l1(&l_old_rcd, &xtprt, &xannotation, &m_cutoff, f_cov, &x_blklist, b_tumor, sf_new_out, Some(&sf_td_new_sites))
        } else if i_rep_type & 4 != 0 {
            let xannotation = Self::construct_interval_tree(sf_rmsk, i_min_copy_len, true, BOUNDARY_EXTND)?;
            self.post_processing_sva(&l_old_rcd, &xtprt, &xannotation, &m_cutoff, &x_blklist, b_tumor, sf_new_out)
        } else {
            let xannotation = Self::construct_interval_tree(sf_rmsk, i_min_copy_len, true, BOUNDARY_EXTND)?;
            self.post_processing_alu(&l_old_rcd, &xtprt, &xannotation, &m_cutoff, f_cov, &x_blklist, sf_new_out)
        }
    }

    fn div_of(xannotation: &XAnnotation, ins_chrm: &str, ins_pos: i64) -> (bool, PyNum) {
        let (b_in_rep, i_pos) = xannotation.is_within_repeat_region_interval_tree(ins_chrm, ins_pos);
        let div = match xannotation.get_div_subfamily(ins_chrm, i_pos) {
            Some((d, _, _, _, _)) => PyNum::Float(d),
            None => PyNum::Int(-1),
        };
        (b_in_rep, div)
    }

    /// `post_processing_SVA`
    #[allow(clippy::too_many_arguments)]
    fn post_processing_sva(&self, l_old_rcd: &[Vec<String>], xtprt: &XTPRTFilter, xannotation: &XAnnotation, m_cutoff: &AfCutoff, x_blklist: &XBlackList, b_tumor: bool, sf_new_out: &str) -> Result<()> {
        let i_dft_cluster_diff_cutoff = gv::two_clip_cluster_diff_cutoff();
        gv::set_two_clip_cluster_diff_cutoff(SVA_CLIP_CLUSTER_DIFF_CUTOFF);
        let mut fout = utils::create(sf_new_out)?;
        for old_rcd in l_old_rcd {
            let rcd = XTEARsltParser::replace_ins_length(old_rcd, REP_SVA_POLYA_START + 100)?;
            let s_rep_supt_type = XTEARsltParser::get_ins_sub_type(&rcd)?;
            let ins_chrm = rcd[0].as_str();
            let ins_pos = pint(&rcd[1])?;
            if x_blklist.fall_in_region(ins_chrm, ins_pos).0 {
                continue;
            }
            let (b_in_rep, div_rate) = Self::div_of(xannotation, ins_chrm, ins_pos);
            let b_with_polya = xtprt.has_polya_signal(&rcd)?;
            let rep_div = gv::REP_DIVERGENT_CUTOFF as f64;
            if s_rep_supt_type == TWO_SIDE {
                if xtprt.is_polya_dominant_two_side_sva(&rcd, REP_SVA_CNS_HEAD)? {
                    continue;
                }
                if !AFConflictFilter::is_qualified_rcd(rcd.last().unwrap(), m_cutoff)? {
                    continue;
                }
                if !b_with_polya {
                    continue;
                }
            } else if s_rep_supt_type == ONE_HALF_SIDE {
                if b_in_rep && div_rate.f() < rep_div {
                    continue;
                }
                if !b_with_polya {
                    continue;
                }
                if xtprt.is_polya_dominant_two_side_sva(&rcd, REP_SVA_CNS_HEAD)? {
                    continue;
                }
            } else if s_rep_supt_type == ONE_SIDE {
                if b_in_rep && div_rate.f() < rep_div {
                    continue;
                }
                if !b_with_polya {
                    continue;
                }
                if (XTEARsltParser::is_orphan_transduction(&rcd)? || XTEARsltParser::is_transduction(&rcd)?) && !b_tumor {
                    continue;
                }
            } else {
                continue;
            }
            if b_in_rep && xtprt.is_two_side_clip_both_polya_sva(&rcd, REP_SVA_POLYA_START)? {
                continue;
            }
            if b_in_rep && xtprt.is_two_side_clip_both_non_polya_sva(&rcd, REP_SVA_POLYA_START)? {
                continue;
            }
            let s_in_rep = if b_in_rep { format!("Fall_in_SVA_copy_{}", div_str(div_rate)) } else { "not_in_SVA_copy".to_string() };
            writeln!(fout, "{}\t{}", utils::rstrip(rcd.last().unwrap()), s_in_rep)?;
        }
        fout.flush()?;
        gv::set_two_clip_cluster_diff_cutoff(i_dft_cluster_diff_cutoff);
        Ok(())
    }

    /// `post_processing_Alu`
    #[allow(clippy::too_many_arguments)]
    fn post_processing_alu(&self, l_old_rcd: &[Vec<String>], xtprt: &XTPRTFilter, xannotation: &XAnnotation, m_cutoff: &AfCutoff, f_cov: f64, x_blklist: &XBlackList, sf_new_out: &str) -> Result<()> {
        let mut fout = utils::create(sf_new_out)?;
        let rep_div = gv::REP_DIVERGENT_CUTOFF as f64;
        for old_rcd in l_old_rcd {
            let rcd = XTEARsltParser::replace_ins_length(old_rcd, REP_ALU_POLYA_START + 15)?;
            let s_rep_supt_type = XTEARsltParser::get_ins_sub_type_alu(&rcd)?;
            let ins_chrm = rcd[0].as_str();
            let ins_pos = pint(&rcd[1])?;
            if x_blklist.fall_in_region(ins_chrm, ins_pos).0 {
                continue;
            }
            let (b_in_rep, div_rate) = Self::div_of(xannotation, ins_chrm, ins_pos);
            let b_with_polya = xtprt.has_polya_signal(&rcd)?;
            if s_rep_supt_type == TWO_SIDE {
                if xtprt.is_polya_dominant_two_side(&rcd)? && xtprt.is_two_side_polya_same_orientation(&rcd)? {
                    continue;
                }
                if !AFConflictFilter::is_qualified_rcd(rcd.last().unwrap(), m_cutoff)? {
                    continue;
                }
                if !b_with_polya && !xtprt.is_two_side_tprt_and_with_polya(&rcd, REP_ALU_MIN_LEN)? {
                    continue;
                }
            } else if s_rep_supt_type == ONE_HALF_SIDE || s_rep_supt_type == ONE_SIDE {
                if b_in_rep && div_rate.f() < rep_div {
                    continue;
                }
                if xtprt.is_polya_dominant_one_side(&rcd, NCLIP_HALF_CUTOFF)? {
                    continue;
                }
                if b_in_rep && !b_with_polya {
                    continue;
                }
            } else {
                continue;
            }
            if b_in_rep && xtprt.is_two_side_clip_and_both_hit_end(&rcd, REP_ALU_POLYA_START)? {
                continue;
            }
            if b_in_rep && !xtprt.hit_consensus_tail(&rcd, REP_ALU_POLYA_START)? {
                continue;
            }
            let b_cov_abnormal = xtprt.cov_is_abnormal(&rcd, f_cov * ABNORMAL_COV_TIMES)?;
            let n_indel_cutoff = crate::py::fmt::py_int(f_cov * INDEL_READS_MAX_RATIO);
            if xtprt.has_enough_indel_reads(&rcd, n_indel_cutoff)? && n_indel_cutoff > 0 && b_cov_abnormal {
                continue;
            }
            let s_in_rep = if b_in_rep { format!("Fall_in_Alu_copy_{}", div_str(div_rate)) } else { "not_in_Alu_copy".to_string() };
            writeln!(fout, "{}\t{}", utils::rstrip(rcd.last().unwrap()), s_in_rep)?;
        }
        fout.flush()?;
        Ok(())
    }

    /// `post_processing_L1`
    #[allow(clippy::too_many_arguments)]
    fn post_processing_l1(
        &self,
        l_old_rcd: &[Vec<String>],
        xtprt: &XTPRTFilter,
        xannotation: &XAnnotation,
        m_cutoff: &AfCutoff,
        f_cov: f64,
        x_blklist: &XBlackList,
        b_tumor: bool,
        sf_new_out: &str,
        sf_td_new: Option<&str>,
    ) -> Result<()> {
        let mut m_td_new: IndexMap<String, IndexMap<String, ()>> = IndexMap::new();
        if let Some(sf) = sf_td_new {
            if utils::is_file(sf) {
                for line in utils::read_lines(sf)? {
                    let f = split_ws(&line);
                    if f.len() < 2 {
                        return Err(anyhow!("IndexError: malformed line in {}", sf));
                    }
                    m_td_new.entry(f[0].to_string()).or_default().insert(f[1].to_string(), ());
                }
            }
        }
        let mut fout = utils::create(sf_new_out)?;
        let mut flog = utils::create(&format!("{}.post_filtering.log", sf_new_out))?;
        let rep_div = gv::REP_DIVERGENT_CUTOFF as f64;
        macro_rules! log {
            ($($arg:tt)*) => {{
                let s = format!($($arg)*);
                println!("{}", s);
                writeln!(flog, "{}", s)?;
            }};
        }
        for old_rcd in l_old_rcd {
            let rcd = XTEARsltParser::replace_ins_length(old_rcd, REP_LINE_POLYA_START + 100)?;
            let s_rep_supt_type = XTEARsltParser::get_ins_sub_type(&rcd)?;
            let ins_chrm = rcd[0].clone();
            let s_ins_pos = rcd[1].clone();
            let ins_pos = pint(&s_ins_pos)?;
            if x_blklist.fall_in_region(&ins_chrm, ins_pos).0 {
                log!("{}:{} is filtered out, because fall in centromere region", ins_chrm, s_ins_pos);
                continue;
            }
            let (b_in_low_div_rep, div_rate, copy_start, copy_end) =
                Self::fall_in_low_div_same_type_rep(xannotation, &ins_chrm, ins_pos, gv::REP_LOW_DIVERGENT_CUTOFF);
            let ref_copy_start = copy_start + L1_BOUNDARY_EXTND;
            let ref_copy_end = copy_end - L1_BOUNDARY_EXTND;
            let mut b_in_rep = false;
            let mut s_in_rep = "not_in_LINE1_copy".to_string();
            let mut b_within_copy = false;
            if div_rate.f() >= 0.0 {
                b_in_rep = true;
                let tmp_pos = ins_pos;
                if tmp_pos >= ref_copy_start && tmp_pos <= ref_copy_end {
                    b_within_copy = true;
                    s_in_rep = format!("Fall_in_LINE1_copy_{}", div_str(div_rate));
                } else if ref_copy_end > 0 && tmp_pos > ref_copy_end {
                    s_in_rep = format!("Left_is_LINE1_copy_{}", div_str(div_rate));
                } else if ref_copy_start > tmp_pos {
                    s_in_rep = format!("Right_is_LINE1_copy_{}", div_str(div_rate));
                } else {
                    s_in_rep = format!("Some_where_close_is_LINE1_copy_{}", div_str(div_rate));
                }
            }
            let s_pass_info = format!("{}\t{}\n", utils::rstrip(rcd.last().unwrap()), s_in_rep);
            if XTEARsltParser::is_orphan_transduction(&rcd)? {
                if s_rep_supt_type != ONE_SIDE || b_tumor {
                    fout.write_all(s_pass_info.as_bytes())?;
                }
                continue;
            } else if XTEARsltParser::is_transduction(&rcd)? {
                if b_in_low_div_rep && b_within_copy && !xtprt.is_two_side_tprt_both_and_both_consistent(&rcd, 5950, 100)? {
                    continue;
                }
                if xtprt.is_polya_dominant_td(&rcd)? || !AFConflictFilter::is_qualified_td(&rcd)? {
                    continue;
                }
                if m_td_new.get(&ins_chrm).map_or(false, |m| m.contains_key(&s_ins_pos)) {
                    if s_rep_supt_type != ONE_SIDE || b_tumor {
                        fout.write_all(s_pass_info.as_bytes())?;
                    }
                    continue;
                }
            }
            let b_with_polya = xtprt.has_polya_signal(&rcd)?;
            let min_ins_len = 200;
            if !b_with_polya && !xtprt.hit_end_of_cns(&rcd)? {
                log!("{}:{} is filtered out, because no any polyA signal detected!", ins_chrm, s_ins_pos);
                continue;
            }
            if s_rep_supt_type == TWO_SIDE {
                if xtprt.is_polya_dominant_two_side(&rcd)?
                    && xtprt.is_two_side_polya_same_orientation(&rcd)?
                    && !xtprt.is_two_side_consist_with_enough_ins_size(&rcd, min_ins_len)?
                {
                    log!("{}:{} is filtered out, because two side polyA dominant!", ins_chrm, s_ins_pos);
                    continue;
                }
                if !xtprt.hit_end_of_cns(&rcd)? {
                    log!("{}:{} is filtered out, because doesn't hit the end of consensus!", ins_chrm, s_ins_pos);
                    continue;
                }
                if !AFConflictFilter::is_qualified_rcd(rcd.last().unwrap(), m_cutoff)? {
                    log!("{}:{} is filtered out, because AF is not qualified!", ins_chrm, s_ins_pos);
                    continue;
                }
                if div_rate.f() >= 0.0 && div_rate.f() < rep_div && !xtprt.is_two_side_tprt_both_and_both_consistent(&rcd, 6000, 200)? {
                    log!("{}:{} is filtered out, because fall in low div L1 region", ins_chrm, s_ins_pos);
                    continue;
                }
            } else if s_rep_supt_type == ONE_HALF_SIDE {
                if div_rate.f() >= 0.0 && div_rate.f() < rep_div {
                    log!("{}:{} is filtered out, because fall in repeat region whose divergent rate is low!", ins_chrm, s_ins_pos);
                    continue;
                }
                if xtprt.is_polya_dominant_two_side(&rcd)? && xtprt.is_two_side_polya_same_orientation(&rcd)? {
                    log!("{}:{} is filtered out, because two side polyA dominant!", ins_chrm, s_ins_pos);
                    continue;
                }
            } else if s_rep_supt_type == ONE_SIDE {
                if div_rate.f() >= 0.0 && div_rate.f() < rep_div {
                    log!("{}:{} is filtered out, because fall in repeat region whose divergent rate is low!", ins_chrm, s_ins_pos);
                    continue;
                }
            } else {
                continue;
            }
            let b_cov_abnormal = xtprt.cov_is_abnormal(&rcd, f_cov * ABNORMAL_COV_TIMES)?;
            if (b_in_rep || b_cov_abnormal)
                && xtprt.is_two_side_clip_and_both_hit_end(&rcd, REP_LINE_POLYA_START)?
                && xtprt.is_two_side_polya_same_orientation(&rcd)?
            {
                log!(
                    "{}:{} is filtered out, because coverage is abnormal and also two side clip both hit end of consensus!",
                    ins_chrm,
                    s_ins_pos
                );
                continue;
            }
            if b_in_rep && !xtprt.hit_consensus_tail(&rcd, REP_LINE_POLYA_START)? {
                log!("{}:{} is filtered out, because fall in repeat region, and none of clip side fall in repeat!", ins_chrm, s_ins_pos);
                continue;
            }
            fout.write_all(s_pass_info.as_bytes())?;
        }
        fout.flush()?;
        flog.flush()?;
        Ok(())
    }
}
