//! Port of `x_population_AF.py` (stand-alone tool, `x_population_AF.py <vcf> <out_csv>`):
//! cumulative allele-frequency counts of the population-specific sites of a population VCF.
//!
//! A record is counted only when exactly one `<POP>_AF` INFO term of the eight populations is
//! present (the scan stops as soon as a second one is seen).  Its AF increments the count of
//! every level `i` with `i * 0.00001 <= AF` (levels 0..n_level-1, n_level =
//! `int(1.0 / 0.00001) + 1` = 100000); records with AF > 0.1 are copied to
//! `<out_csv>.high_af.sites`.  Like the script, an INFO term without "=" (before the scan
//! stops) or a record with fewer than 8 columns raises an error.

use super::{create_file, open_text, py_float, py_rstrip};
use crate::py::fmt::float_repr;
use anyhow::{anyhow, Result};
use std::ffi::{OsStr, OsString};
use std::io::Write;

/// `l_pop_code`
pub const L_POP_CODE: [&str; 8] = ["AFR", "AMR", "EAS", "EUR", "SAS", "WAS", "OCN", "CAS"];
/// `f_base`
pub const F_BASE: f64 = 0.00001;
/// `n_high_af`
pub const N_HIGH_AF: f64 = 0.1;

/// `n_level = int(1.0 / f_base) + 1` (1.0 / 0.00001 = 99999.99999999999, so 100000).
pub fn n_level() -> usize {
    crate::py::fmt::py_int(1.0 / F_BASE) as usize + 1
}

/// Number of leading levels `add_a_site` increments for `f_af`: its loop
/// `for i in range(n_level): if float(i) * f_base > f_af: break` stops at the first level whose
/// value exceeds the AF.  `float(i) * f_base` is non-decreasing in `i` (exact integers times a
/// positive constant, correctly rounded), so the break point is found by bisection; NaN never
/// breaks.
pub fn n_levels_hit(f_af: f64, n_level: usize) -> usize {
    let (mut lo, mut hi) = (0usize, n_level);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if (mid as f64) * F_BASE > f_af {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    lo
}

/// The population code and AF value of a record when exactly one `<POP>_AF` term is present.
fn single_population_af(sinfo: &str) -> Result<Option<(&'static str, &str)>> {
    let mut n_cnt = 0;
    let mut s_tmp_pop_code = "";
    let mut s_tmp_af = "";
    for sterm in sinfo.split(';') {
        let mut info_sub_fields = sterm.split('=');
        let s_term_id = info_sub_fields.next().unwrap_or("");
        let s_term_value = info_sub_fields.next().ok_or_else(|| anyhow!("IndexError: list index out of range (INFO term {:?} has no '=')", sterm))?;
        for s_pop_code in L_POP_CODE {
            // s_pop_code + "_AF" == s_term_id
            if s_term_id.strip_suffix("_AF") == Some(s_pop_code) {
                n_cnt += 1;
                s_tmp_pop_code = s_pop_code;
                s_tmp_af = s_term_value;
            }
        }
        if n_cnt > 1 {
            break;
        }
    }
    Ok(if n_cnt == 1 { Some((s_tmp_pop_code, s_tmp_af)) } else { None })
}

/// `cnt_population_AF(sf_vcf, sf_out_csv)`
pub fn cnt_population_af(sf_vcf: &OsStr, sf_out_csv: &OsStr) -> Result<()> {
    let n_level = n_level();
    // histogram of the number of levels hit, per population (in l_pop_code order)
    let mut m_hist: Vec<Vec<u64>> = vec![vec![0u64; n_level + 1]; L_POP_CODE.len()];

    let mut sf_sites = OsString::from(sf_out_csv);
    sf_sites.push(".high_af.sites");
    let text = open_text(sf_vcf)?;
    {
        let mut fout_sites = create_file(&sf_sites)?;
        for line in &text.lines {
            if line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            let sinfo = fields.get(7).ok_or_else(|| anyhow!("IndexError: list index out of range (record with fewer than 8 columns)"))?;
            if let Some((s_pop, s_af)) = single_population_af(sinfo)? {
                let f_tmp_af = py_float(s_af)?;
                let i_pop = L_POP_CODE.iter().position(|p| *p == s_pop).unwrap_or(0);
                m_hist[i_pop][n_levels_hit(f_tmp_af, n_level)] += 1;
                if f_tmp_af > N_HIGH_AF {
                    fout_sites.write_all(py_rstrip(line).as_bytes())?;
                    fout_sites.write_all(b"\n")?;
                }
            }
        }
        fout_sites.flush()?;
    }
    text.check()?;

    // level strings: str(float(i) * f_base)
    let l_af: Vec<String> = (0..n_level).map(|i| float_repr((i as f64) * F_BASE)).collect();
    let mut fout_csv = create_file(sf_out_csv)?;
    fout_csv.write_all(b"AF,Population,Count\n")?;
    for (i_pop, pop_code) in L_POP_CODE.iter().enumerate() {
        // count of level i = number of sites hitting more than i levels
        let hist = &m_hist[i_pop];
        let mut n_remaining: u64 = hist.iter().sum();
        for (i, s_af) in l_af.iter().enumerate() {
            n_remaining -= hist[i];
            writeln!(fout_csv, "{},{},{}", s_af, pop_code, n_remaining)?;
        }
    }
    fout_csv.flush()?;
    Ok(())
}

/// The script's `__main__` (`sys.argv[1]`, `sys.argv[2]`; further arguments are ignored).
pub fn main_with_args(args: &[OsString]) -> Result<()> {
    let sf_vcf = args.get(1).ok_or_else(|| anyhow!("IndexError: list index out of range (usage: xtea_population_af <vcf> <out_csv>)"))?;
    let sf_out = args.get(2).ok_or_else(|| anyhow!("IndexError: list index out of range (usage: xtea_population_af <vcf> <out_csv>)"))?;
    cnt_population_af(sf_vcf, sf_out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The loop of the script's `add_a_site`, counting the increments.
    fn naive(f_af: f64, n_level: usize) -> usize {
        let mut n = 0;
        for i in 0..n_level {
            if (i as f64) * F_BASE > f_af {
                break;
            }
            n += 1;
        }
        n
    }

    #[test]
    fn levels() {
        assert_eq!(n_level(), 100000);
        let n = n_level();
        let mut vals = vec![-1.0, -0.0, 0.0, 1e-6, 1e-5, 2e-5, 3e-5, 3.0000000000000004e-05, 0.1, 0.10000000000000002, 0.5, 0.99999, 0.999995, 1.0, 2.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN];
        for i in 0..2000 {
            // values exactly on and next to the level grid
            let x = (i as f64) * 0.00037;
            vals.push(x);
            vals.push(((i * 7) as f64) * F_BASE);
            vals.push(f64::from_bits((((i * 7) as f64) * F_BASE).to_bits() + 1));
            vals.push(f64::from_bits((((i * 7) as f64) * F_BASE).to_bits().saturating_sub(1)));
        }
        for v in vals {
            assert_eq!(n_levels_hit(v, n), naive(v, n), "af={:?}", v);
        }
    }

    #[test]
    fn af_strings() {
        // Python: [str(float(i) * 0.00001) for i in (0, 1, 3, 7, 10, 99999)]
        let got: Vec<String> = [0usize, 1, 3, 7, 10, 99999].iter().map(|&i| float_repr((i as f64) * F_BASE)).collect();
        assert_eq!(got, vec!["0.0", "1e-05", "3.0000000000000004e-05", "7.000000000000001e-05", "0.0001", "0.99999"]);
    }

    #[test]
    fn info_scan() {
        assert_eq!(single_population_af("AN=10;EUR_AF=0.5;GENE=x").unwrap(), Some(("EUR", "0.5")));
        assert_eq!(single_population_af("AF=0.1;EUR_AF=0.5;AFR_AF=0.1;BAD").unwrap(), None);
        assert!(single_population_af("AF=0.1;BAD;EUR_AF=0.5").is_err());
        assert_eq!(single_population_af("EUR_AF=0.5=7").unwrap(), Some(("EUR", "0.5")));
        assert_eq!(single_population_af("EUR_AF=0.5;EUR_AF=0.2;x").unwrap(), None);
        assert_eq!(single_population_af("EUR_AFX=0.5;xEUR_AF=1").unwrap(), None);
    }
}
