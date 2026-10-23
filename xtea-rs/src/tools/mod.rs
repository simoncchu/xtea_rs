//! Ports of the standalone xTea tools, each of which has its own `__main__` in Python:
//!
//! * [`vcf_merger`] - `x_vcf_merger.py` (merge the per-repeat VCFs of a sample; population VCF)
//! * [`cmp`] - `x_cmp.py` (compare calls against a benchmark; F1 / sensitivity tables)
//! * [`population_af`] - `x_population_AF.py` (population allele-frequency counting)
//!
//! This file holds the small Python-semantics helpers the three tools share (text-mode file
//! reading, `str.split()`, `int()`, ...).  They are deliberately stricter than the generic
//! helpers in `crate::utils` so that the tools behave exactly like the scripts:
//! `str.split()`/`str.rstrip()` treat U+001C..U+001F and the non-ASCII Unicode spaces as
//! whitespace, and text files are decoded as strict UTF-8 with universal newlines.

pub mod cmp;
pub mod population_af;
pub mod vcf_merger;

use anyhow::{anyhow, Result};
use std::collections::BTreeSet;
use std::fs::File;
use std::io::BufWriter;
use std::ops::Bound;
use std::path::Path;

pub use crate::py::text::{
    is_py_space, open_text, py_float, py_int, py_repr_str, py_rstrip, py_split_ws, read_text_lines, split_universal_newlines, TextFile,
};

/// Python `s[start:]` (slicing counts code points, not bytes).
pub fn py_slice_from(s: &str, start: usize) -> &str {
    match s.char_indices().nth(start) {
        Some((i, _)) => &s[i..],
        None => "",
    }
}

/// Python `s[1:-1]`.
pub fn py_slice_inner(s: &str) -> &str {
    let n = s.chars().count();
    if n <= 2 {
        return "";
    }
    let start = s.char_indices().nth(1).map(|(i, _)| i).unwrap_or(0);
    let end = s.char_indices().nth(n - 1).map(|(i, _)| i).unwrap_or(s.len());
    &s[start..end]
}

/// Python `v[i]` for a non-negative index (IndexError -> Err).
#[inline]
pub fn at<'a>(v: &[&'a str], i: usize) -> Result<&'a str> {
    v.get(i).copied().ok_or_else(|| anyhow!("IndexError: list index out of range (index {})", i))
}

/// Python `v[-k]` (k >= 1).
#[inline]
pub fn at_back<'a>(v: &[&'a str], k: usize) -> Result<&'a str> {
    if k == 0 || k > v.len() {
        return Err(anyhow!("IndexError: list index out of range (index -{})", k));
    }
    Ok(v[v.len() - k])
}

/// A Python value that may be `None` where the script would raise a `TypeError` when using it
/// (`None + "..."`, `open(None)`, ...).
#[inline]
pub fn need<'a>(v: Option<&'a str>, what: &str) -> Result<&'a str> {
    v.ok_or_else(|| anyhow!("TypeError: {} is None (option not given)", what))
}

/// Python `a / b` for floats (ZeroDivisionError -> Err).
#[inline]
pub fn py_fdiv(a: f64, b: f64) -> Result<f64> {
    if b == 0.0 {
        return Err(anyhow!("ZeroDivisionError: float division by zero"));
    }
    Ok(a / b)
}

/// The smallest key `k` with `center - ext <= k < center + ext`, i.e. the first hit of the
/// scripts' `for i in range(-ext, ext): if center + i in keys` loops (no hit when `ext <= 0`).
pub fn first_key_in_window(keys: &BTreeSet<i64>, center: i64, ext: i64) -> Option<i64> {
    let lo = center as i128 - ext as i128;
    let hi = center as i128 + ext as i128;
    if hi <= lo || lo > i64::MAX as i128 || hi <= i64::MIN as i128 {
        return None;
    }
    let lower = Bound::Included(lo.max(i64::MIN as i128) as i64);
    let upper = if hi > i64::MAX as i128 { Bound::Unbounded } else { Bound::Excluded(hi as i64) };
    keys.range((lower, upper)).next().copied()
}

/// Python `open(path, "w")`.
pub fn create_file<P: AsRef<Path>>(path: P) -> Result<BufWriter<File>> {
    let p = path.as_ref();
    let f = File::create(p).map_err(|e| anyhow!("cannot create {}: {}", p.display(), e))?;
    Ok(BufWriter::new(f))
}

/// Python `os.remove(path)` (raises when the file is missing).
pub fn remove_file<P: AsRef<Path>>(path: P) -> Result<()> {
    let p = path.as_ref();
    std::fs::remove_file(p).map_err(|e| anyhow!("cannot remove {}: {}", p.display(), e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slicing_counts_code_points() {
        assert_eq!(py_slice_from("TD_SRC=xyz", 7), "xyz");
        assert_eq!(py_slice_from("TD_SRC", 7), "");
        assert_eq!(py_slice_from("TD_SRÇ=é1", 7), "é1");
        assert_eq!(py_slice_inner("<INS:ME:ALU>"), "INS:ME:ALU");
        assert_eq!(py_slice_inner("ab"), "");
        assert_eq!(py_slice_inner("a"), "");
        assert_eq!(py_slice_inner(""), "");
        assert_eq!(py_slice_inner("<é>"), "é");
    }

    #[test]
    fn float_division_by_zero() {
        assert!(py_fdiv(1.0, 0.0).is_err());
        assert!(py_fdiv(0.0, -0.0).is_err());
    }

    #[test]
    fn window_lookup_matches_range_loop() {
        let keys: BTreeSet<i64> = [5, 10, 20, 21, 50].into_iter().collect();
        for center in -10..70 {
            for ext in -2..25 {
                let mut exp = None;
                for i in -ext..ext {
                    if keys.contains(&(center + i)) {
                        exp = Some(center + i);
                        break;
                    }
                }
                assert_eq!(first_key_in_window(&keys, center, ext), exp, "center={} ext={}", center, ext);
            }
        }
        let big: BTreeSet<i64> = [i64::MAX, i64::MIN].into_iter().collect();
        assert_eq!(first_key_in_window(&big, i64::MAX, 1), Some(i64::MAX));
        assert_eq!(first_key_in_window(&big, i64::MIN, 1), Some(i64::MIN));
        // [0 - MAX, 0 + MAX) = [MIN + 1, MAX): neither extreme is inside
        assert_eq!(first_key_in_window(&big, 0, i64::MAX), None);
        assert_eq!(first_key_in_window(&big, -1, i64::MAX), Some(i64::MIN));
    }


}
