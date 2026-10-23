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

/// Rewrite argv so that clap sees the values optparse sees.  optparse takes the rest of a
/// short-option cluster verbatim as the option's value ("-p=out" gives "=out", "-Sy=7" gives
/// -S and "=7"), while clap drops a leading "="; such attached values are moved into their own
/// argument, which clap then takes unchanged (the options allow hyphen values).  Arguments
/// that optparse consumes as option values, and everything after "--", are left alone.  The
/// option table is taken from the clap definition of `T`, and long options are resolved by
/// unique prefix like optparse does.
pub fn optparse_compatible_argv<T: clap::CommandFactory>(args: Vec<std::ffi::OsString>) -> Vec<std::ffi::OsString> {
    use std::ffi::OsString;
    let cmd = T::command();
    // (short, long, takes a value); optparse's -h/--help is a flag
    let mut specs: Vec<(Option<char>, Option<String>, bool)> = vec![(Some('h'), Some("help".to_string()), false)];
    for a in cmd.get_arguments() {
        if a.get_short().is_none() && a.get_long().is_none() {
            continue; // positional
        }
        specs.push((a.get_short(), a.get_long().map(|s| s.to_string()), a.get_action().takes_values()));
    }
    let long_takes_value = |body: &str| -> bool {
        if let Some(s) = specs.iter().find(|s| s.1.as_deref() == Some(body)) {
            return s.2;
        }
        let cands: Vec<&(Option<char>, Option<String>, bool)> = specs.iter().filter(|s| s.1.as_deref().is_some_and(|l| l.starts_with(body))).collect();
        cands.len() == 1 && cands[0].2
    };

    let mut out: Vec<OsString> = Vec::with_capacity(args.len() + 2);
    let mut iter = args.into_iter();
    if let Some(prog) = iter.next() {
        out.push(prog);
    }
    let mut value_next = false;
    let mut positional_only = false;
    for a in iter {
        if positional_only || value_next {
            value_next = false;
            out.push(a);
            continue;
        }
        let s = match a.to_str() {
            Some(s) => s.to_string(),
            None => {
                out.push(a);
                continue;
            }
        };
        if s == "--" {
            positional_only = true;
        } else if let Some(body) = s.strip_prefix("--") {
            value_next = !body.contains('=') && long_takes_value(body);
        } else if s.len() > 1 && s.starts_with('-') {
            let body = &s[1..];
            for (i, c) in body.char_indices() {
                match specs.iter().find(|sp| sp.0 == Some(c)) {
                    Some(sp) if sp.2 => {
                        let k = 1 + i + c.len_utf8();
                        if k == s.len() {
                            value_next = true;
                        } else if s[k..].starts_with('=') {
                            out.push(OsString::from(&s[..k]));
                            out.push(OsString::from(&s[k..]));
                            break;
                        }
                        out.push(a.clone());
                        break;
                    }
                    Some(_) if i + c.len_utf8() < body.len() => continue, // a flag inside the cluster
                    _ => {
                        out.push(a.clone()); // last flag, or an unknown option (clap reports it)
                        break;
                    }
                }
            }
            continue;
        }
        out.push(a);
    }
    out
}

/// `T::parse_from` on the process arguments made optparse-compatible
/// ([`optparse_compatible_argv`]).
pub fn parse_optparse_args<T: clap::Parser>() -> T {
    T::parse_from(optparse_compatible_argv::<T>(std::env::args_os().collect()))
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

    #[test]
    fn optparse_argv_normalisation() {
        use clap::Parser;
        use vcf_merger::Opts;
        let norm = |v: &[&str]| -> Vec<String> {
            let a: Vec<std::ffi::OsString> = std::iter::once("prog").chain(v.iter().copied()).map(Into::into).collect();
            optparse_compatible_argv::<Opts>(a).into_iter().skip(1).map(|s| s.into_string().unwrap()).collect()
        };
        assert_eq!(norm(&["-p=out"]), ["-p", "=out"]);
        assert_eq!(norm(&["-Sy=7", "-x", "l"]), ["-Sy", "=7", "-x", "l"]);
        assert_eq!(norm(&["-SPp=x"]), ["-SPp", "=x"]);
        assert_eq!(norm(&["-pout", "-Sw5"]), ["-pout", "-Sw5"]);
        // values consumed by the previous option, unknown options and "--" are left alone
        assert_eq!(norm(&["-i", "-p=x"]), ["-i", "-p=x"]);
        assert_eq!(norm(&["--rep", "-y=3"]), ["--rep", "-y=3"]);
        assert_eq!(norm(&["--path=-p=x", "-w=1"]), ["--path=-p=x", "-w", "=1"]);
        assert_eq!(norm(&["--", "-p=x"]), ["--", "-p=x"]);
        assert_eq!(norm(&["-Zp=x"]), ["-Zp=x"]);
        assert_eq!(norm(&["--sample", "-p=x"]), ["--sample", "-p", "=x"]);
        assert_eq!(norm(&["--re", "-p=x"]), ["--re", "-p", "=x"]); // ambiguous prefix: not a value-taking option
        assert_eq!(norm(&["-p"]), ["-p"]);
        // and clap then sees optparse's values
        let o = Opts::try_parse_from(optparse_compatible_argv::<Opts>(["prog", "-S", "-p=out", "-i", "-x"].iter().map(Into::into).collect())).unwrap();
        assert_eq!(o.wfolder.as_deref(), Some("=out"));
        assert_eq!(o.sample_info.as_deref(), Some("-x"));
        assert!(Opts::try_parse_from(optparse_compatible_argv::<Opts>(["prog", "-y=7"].iter().map(Into::into).collect())).is_err());
    }

}
