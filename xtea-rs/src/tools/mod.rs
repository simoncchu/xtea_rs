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

/// Python's whitespace test used by `str.split()` / `str.strip()` (`Py_UNICODE_ISSPACE`):
/// Unicode White_Space plus the ASCII separators U+001C..U+001F.
#[inline]
pub fn is_py_space(c: char) -> bool {
    c.is_whitespace() || ('\u{1c}'..='\u{1f}').contains(&c)
}

/// Python `s.split()` (no separator): split on runs of whitespace, dropping empty fields.
pub fn py_split_ws(s: &str) -> Vec<&str> {
    s.split(is_py_space).filter(|t| !t.is_empty()).collect()
}

/// Python `s.rstrip()`.
#[inline]
pub fn py_rstrip(s: &str) -> &str {
    s.trim_end_matches(is_py_space)
}

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

/// Decimal value of a Unicode decimal digit (category Nd), like `Py_UNICODE_TODECIMAL`.
fn unicode_decimal(c: char) -> Option<u8> {
    let u = c as u32;
    let i = ND_ZEROS.partition_point(|&z| z <= u);
    if i == 0 {
        return None;
    }
    let off = u - ND_ZEROS[i - 1];
    if off < 10 {
        Some(off as u8)
    } else {
        None
    }
}

/// The text `int()` / `float()` actually parse: CPython's
/// `_PyUnicode_TransformDecimalAndSpaceToASCII` (non-ASCII whitespace becomes ' ', Unicode
/// decimal digits become ASCII digits, any other non-ASCII character is invalid) followed by
/// the strip of ASCII whitespace (`Py_ISSPACE`: note that U+001C..U+001F are *not* stripped
/// here, unlike `str.strip()`).
fn py_numeric_ascii(s: &str) -> Option<String> {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        if (c as u32) < 127 {
            out.push(c);
        } else if is_py_space(c) {
            out.push(' ');
        } else {
            // any other non-ASCII character makes the literal invalid
            let d = unicode_decimal(c)?;
            out.push((b'0' + d) as char);
        }
    }
    let t = out.trim_matches(|c: char| matches!(c, ' ' | '\t' | '\n' | '\x0b' | '\x0c' | '\r'));
    Some(t.to_string())
}

/// Python `int(s)` for a `str` in base 10: surrounding whitespace, an optional sign, single
/// underscores between digits and Unicode decimal digits are accepted (values must fit in an
/// i64, Python's ints are unbounded).
pub fn py_int(s: &str) -> Result<i64> {
    let err = || anyhow!("ValueError: invalid literal for int() with base 10: {}", py_repr_str(s));
    let t = py_numeric_ascii(s).ok_or_else(err)?;
    let t = t.as_str();
    let (neg, body) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    let b = body.as_bytes();
    if b.is_empty() || b[0] == b'_' || b[b.len() - 1] == b'_' {
        return Err(err());
    }
    let mut v: i128 = 0;
    let mut prev_us = false;
    for &c in b {
        if c == b'_' {
            if prev_us {
                return Err(err());
            }
            prev_us = true;
            continue;
        }
        prev_us = false;
        if !c.is_ascii_digit() {
            return Err(err());
        }
        v = v * 10 + (c - b'0') as i128;
        if v > (i64::MAX as i128) + 1 {
            return Err(anyhow!("integer {} does not fit in 64 bits", s));
        }
    }
    let v = if neg { -v } else { v };
    i64::try_from(v).map_err(|_| anyhow!("integer {} does not fit in 64 bits", s))
}

/// Python `float(s)` for a `str` (whitespace, `inf`/`nan`, underscores between digits, Unicode
/// decimal digits).
pub fn py_float(s: &str) -> Result<f64> {
    let err = || anyhow!("ValueError: could not convert string to float: {}", py_repr_str(s));
    let t = py_numeric_ascii(s).ok_or_else(err)?;
    crate::py::optparse::parse_float(&t).map_err(|_| err())
}

/// Python `str.isprintable()` for one character.
fn is_py_printable(c: char) -> bool {
    let u = c as u32;
    let i = NON_PRINTABLE.partition_point(|&(lo, _)| lo <= u);
    !(i > 0 && u <= NON_PRINTABLE[i - 1].1)
}

/// Python `repr(str)` (CPython `unicode_repr`: quote choice, `\t \n \r` escapes, `\xhh`,
/// `\uhhhh`, `\Uhhhhhhhh` for every non-printable character).
pub fn py_repr_str(s: &str) -> String {
    use std::fmt::Write as _;
    let quote = if s.contains('\'') && !s.contains('"') { '"' } else { '\'' };
    let mut out = String::with_capacity(s.len() + 2);
    out.push(quote);
    for c in s.chars() {
        let u = c as u32;
        if c == quote || c == '\\' {
            out.push('\\');
            out.push(c);
        } else if c == '\t' {
            out.push_str("\\t");
        } else if c == '\n' {
            out.push_str("\\n");
        } else if c == '\r' {
            out.push_str("\\r");
        } else if u < 0x20 || u == 0x7f {
            let _ = write!(out, "\\x{:02x}", u);
        } else if u < 0x7f || is_py_printable(c) {
            out.push(c);
        } else if u <= 0xff {
            let _ = write!(out, "\\x{:02x}", u);
        } else if u <= 0xffff {
            let _ = write!(out, "\\u{:04x}", u);
        } else {
            let _ = write!(out, "\\U{:08x}", u);
        }
    }
    out.push(quote);
    out
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

/// Lines of a text file as Python 3 yields them when iterating over `open(path)`: the file is
/// decoded as strict UTF-8 (the scripts run in UTF-8 mode; invalid bytes raise
/// UnicodeDecodeError) and universal newlines apply ("\r\n" and a lone "\r" become "\n").
/// Every line keeps its "\n", except possibly the last one.
pub fn read_text_lines<P: AsRef<Path>>(path: P) -> Result<Vec<String>> {
    open_text(path)?.into_lines()
}

/// A text file opened like Python's `open(path)` (see [`open_text`]).
pub struct TextFile {
    /// The lines iteration yields (all of them, or those before a decoding error).
    pub lines: Vec<String>,
    decode_error: Option<String>,
}

impl TextFile {
    /// The UnicodeDecodeError the iteration raises after the yielded lines, if any.
    pub fn check(&self) -> Result<()> {
        match &self.decode_error {
            Some(e) => Err(anyhow!("{}", e)),
            None => Ok(()),
        }
    }

    /// All lines, or the decoding error.
    pub fn into_lines(self) -> Result<Vec<String>> {
        self.check()?;
        Ok(self.lines)
    }
}

/// Size of the chunks `TextIOWrapper` decodes (`_CHUNK_SIZE`, read with `read1`).
const PY_TEXT_CHUNK: usize = 8192;

/// Open a text file like Python's `open(path)`: a missing file fails at once (the scripts open
/// their input before creating their outputs), while a UTF-8 decoding error only surfaces
/// during iteration.  CPython decodes 8192-byte chunks and raises when it decodes the chunk in
/// which the invalid byte is detected, so the lines completed by the previous chunks are still
/// yielded first ([`TextFile::lines`]; a trailing incomplete character and a trailing "\r" of
/// the decoded text are held back by the incremental decoders), then [`TextFile::check`] fails.
pub fn open_text<P: AsRef<Path>>(path: P) -> Result<TextFile> {
    let p = path.as_ref();
    let data = std::fs::read(p).map_err(|e| anyhow!("cannot open {}: {}", p.display(), e))?;
    match String::from_utf8(data) {
        Ok(text) => Ok(TextFile { lines: split_universal_newlines(&text), decode_error: None }),
        Err(e) => {
            let data = e.as_bytes();
            let ue = e.utf8_error();
            let valid = ue.valid_up_to();
            let decoded_end = match ue.error_len() {
                // incomplete sequence at the end of the file: raised by the final decode
                None => valid,
                Some(n) => {
                    // CPython's stateful decoder holds back a truncated surrogate "\xED[\xA0-\xBF]"
                    // that ends the decoded data ("Truncated surrogate code in range D800-DFFF")
                    let lead = data[valid];
                    let held_surrogate = lead == 0xED
                        && n == 1
                        && valid + 1 < data.len()
                        && (0xA0..=0xBF).contains(&data[valid + 1])
                        && ((valid + 2) % PY_TEXT_CHUNK == 0 || valid + 2 == data.len());
                    // byte at which the decoder notices the error: the bad lead byte itself, or
                    // the first byte that cannot continue a valid lead byte
                    let trigger = if held_surrogate {
                        valid + 2
                    } else if (0xC2..=0xF4).contains(&lead) {
                        valid + n
                    } else {
                        valid
                    };
                    let chunk_start = if trigger == data.len() { trigger } else { trigger / PY_TEXT_CHUNK * PY_TEXT_CHUNK };
                    if chunk_start >= valid {
                        valid
                    } else {
                        let mut b = chunk_start;
                        while b > 0 && (data[b] & 0xC0) == 0x80 {
                            b -= 1;
                        }
                        b
                    }
                }
            };
            let text = std::str::from_utf8(&data[..decoded_end]).unwrap_or("");
            let text = text.strip_suffix('\r').unwrap_or(text);
            let mut lines = split_universal_newlines(text);
            if lines.last().is_some_and(|l| !l.ends_with('\n')) {
                lines.pop(); // the incomplete last line is never yielded
            }
            Ok(TextFile { lines, decode_error: Some(format!("UnicodeDecodeError in {}: {}", p.display(), ue)) })
        }
    }
}

/// Universal-newline line splitting (see [`read_text_lines`]).
pub fn split_universal_newlines(text: &str) -> Vec<String> {
    let b = text.as_bytes();
    let mut lines = Vec::new();
    let mut start = 0;
    let mut i = 0;
    while i < b.len() {
        match b[i] {
            b'\n' | b'\r' => {
                let mut line = String::with_capacity(i - start + 1);
                line.push_str(&text[start..i]);
                line.push('\n');
                lines.push(line);
                if b[i] == b'\r' && i + 1 < b.len() && b[i + 1] == b'\n' {
                    i += 1;
                }
                i += 1;
                start = i;
            }
            _ => i += 1,
        }
    }
    if start < b.len() {
        lines.push(text[start..].to_string());
    }
    lines
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

/// Code point ranges for which CPython 3.8.20's `str.isprintable()` is false (Unicode 12.1.0), generated
/// with `unicodedata`; used by [`py_repr_str`].
#[rustfmt::skip]
const NON_PRINTABLE: [(u32, u32); 668] = [
    (0x0, 0x1F), (0x7F, 0xA0), (0xAD, 0xAD), (0x378, 0x379), (0x380, 0x383), (0x38B, 0x38B),
    (0x38D, 0x38D), (0x3A2, 0x3A2), (0x530, 0x530), (0x557, 0x558), (0x58B, 0x58C), (0x590, 0x590),
    (0x5C8, 0x5CF), (0x5EB, 0x5EE), (0x5F5, 0x605), (0x61C, 0x61D), (0x6DD, 0x6DD), (0x70E, 0x70F),
    (0x74B, 0x74C), (0x7B2, 0x7BF), (0x7FB, 0x7FC), (0x82E, 0x82F), (0x83F, 0x83F), (0x85C, 0x85D),
    (0x85F, 0x85F), (0x86B, 0x89F), (0x8B5, 0x8B5), (0x8BE, 0x8D2), (0x8E2, 0x8E2), (0x984, 0x984),
    (0x98D, 0x98E), (0x991, 0x992), (0x9A9, 0x9A9), (0x9B1, 0x9B1), (0x9B3, 0x9B5), (0x9BA, 0x9BB),
    (0x9C5, 0x9C6), (0x9C9, 0x9CA), (0x9CF, 0x9D6), (0x9D8, 0x9DB), (0x9DE, 0x9DE), (0x9E4, 0x9E5),
    (0x9FF, 0xA00), (0xA04, 0xA04), (0xA0B, 0xA0E), (0xA11, 0xA12), (0xA29, 0xA29), (0xA31, 0xA31),
    (0xA34, 0xA34), (0xA37, 0xA37), (0xA3A, 0xA3B), (0xA3D, 0xA3D), (0xA43, 0xA46), (0xA49, 0xA4A),
    (0xA4E, 0xA50), (0xA52, 0xA58), (0xA5D, 0xA5D), (0xA5F, 0xA65), (0xA77, 0xA80), (0xA84, 0xA84),
    (0xA8E, 0xA8E), (0xA92, 0xA92), (0xAA9, 0xAA9), (0xAB1, 0xAB1), (0xAB4, 0xAB4), (0xABA, 0xABB),
    (0xAC6, 0xAC6), (0xACA, 0xACA), (0xACE, 0xACF), (0xAD1, 0xADF), (0xAE4, 0xAE5), (0xAF2, 0xAF8),
    (0xB00, 0xB00), (0xB04, 0xB04), (0xB0D, 0xB0E), (0xB11, 0xB12), (0xB29, 0xB29), (0xB31, 0xB31),
    (0xB34, 0xB34), (0xB3A, 0xB3B), (0xB45, 0xB46), (0xB49, 0xB4A), (0xB4E, 0xB55), (0xB58, 0xB5B),
    (0xB5E, 0xB5E), (0xB64, 0xB65), (0xB78, 0xB81), (0xB84, 0xB84), (0xB8B, 0xB8D), (0xB91, 0xB91),
    (0xB96, 0xB98), (0xB9B, 0xB9B), (0xB9D, 0xB9D), (0xBA0, 0xBA2), (0xBA5, 0xBA7), (0xBAB, 0xBAD),
    (0xBBA, 0xBBD), (0xBC3, 0xBC5), (0xBC9, 0xBC9), (0xBCE, 0xBCF), (0xBD1, 0xBD6), (0xBD8, 0xBE5),
    (0xBFB, 0xBFF), (0xC0D, 0xC0D), (0xC11, 0xC11), (0xC29, 0xC29), (0xC3A, 0xC3C), (0xC45, 0xC45),
    (0xC49, 0xC49), (0xC4E, 0xC54), (0xC57, 0xC57), (0xC5B, 0xC5F), (0xC64, 0xC65), (0xC70, 0xC76),
    (0xC8D, 0xC8D), (0xC91, 0xC91), (0xCA9, 0xCA9), (0xCB4, 0xCB4), (0xCBA, 0xCBB), (0xCC5, 0xCC5),
    (0xCC9, 0xCC9), (0xCCE, 0xCD4), (0xCD7, 0xCDD), (0xCDF, 0xCDF), (0xCE4, 0xCE5), (0xCF0, 0xCF0),
    (0xCF3, 0xCFF), (0xD04, 0xD04), (0xD0D, 0xD0D), (0xD11, 0xD11), (0xD45, 0xD45), (0xD49, 0xD49),
    (0xD50, 0xD53), (0xD64, 0xD65), (0xD80, 0xD81), (0xD84, 0xD84), (0xD97, 0xD99), (0xDB2, 0xDB2),
    (0xDBC, 0xDBC), (0xDBE, 0xDBF), (0xDC7, 0xDC9), (0xDCB, 0xDCE), (0xDD5, 0xDD5), (0xDD7, 0xDD7),
    (0xDE0, 0xDE5), (0xDF0, 0xDF1), (0xDF5, 0xE00), (0xE3B, 0xE3E), (0xE5C, 0xE80), (0xE83, 0xE83),
    (0xE85, 0xE85), (0xE8B, 0xE8B), (0xEA4, 0xEA4), (0xEA6, 0xEA6), (0xEBE, 0xEBF), (0xEC5, 0xEC5),
    (0xEC7, 0xEC7), (0xECE, 0xECF), (0xEDA, 0xEDB), (0xEE0, 0xEFF), (0xF48, 0xF48), (0xF6D, 0xF70),
    (0xF98, 0xF98), (0xFBD, 0xFBD), (0xFCD, 0xFCD), (0xFDB, 0xFFF), (0x10C6, 0x10C6), (0x10C8, 0x10CC),
    (0x10CE, 0x10CF), (0x1249, 0x1249), (0x124E, 0x124F), (0x1257, 0x1257), (0x1259, 0x1259), (0x125E, 0x125F),
    (0x1289, 0x1289), (0x128E, 0x128F), (0x12B1, 0x12B1), (0x12B6, 0x12B7), (0x12BF, 0x12BF), (0x12C1, 0x12C1),
    (0x12C6, 0x12C7), (0x12D7, 0x12D7), (0x1311, 0x1311), (0x1316, 0x1317), (0x135B, 0x135C), (0x137D, 0x137F),
    (0x139A, 0x139F), (0x13F6, 0x13F7), (0x13FE, 0x13FF), (0x1680, 0x1680), (0x169D, 0x169F), (0x16F9, 0x16FF),
    (0x170D, 0x170D), (0x1715, 0x171F), (0x1737, 0x173F), (0x1754, 0x175F), (0x176D, 0x176D), (0x1771, 0x1771),
    (0x1774, 0x177F), (0x17DE, 0x17DF), (0x17EA, 0x17EF), (0x17FA, 0x17FF), (0x180E, 0x180F), (0x181A, 0x181F),
    (0x1879, 0x187F), (0x18AB, 0x18AF), (0x18F6, 0x18FF), (0x191F, 0x191F), (0x192C, 0x192F), (0x193C, 0x193F),
    (0x1941, 0x1943), (0x196E, 0x196F), (0x1975, 0x197F), (0x19AC, 0x19AF), (0x19CA, 0x19CF), (0x19DB, 0x19DD),
    (0x1A1C, 0x1A1D), (0x1A5F, 0x1A5F), (0x1A7D, 0x1A7E), (0x1A8A, 0x1A8F), (0x1A9A, 0x1A9F), (0x1AAE, 0x1AAF),
    (0x1ABF, 0x1AFF), (0x1B4C, 0x1B4F), (0x1B7D, 0x1B7F), (0x1BF4, 0x1BFB), (0x1C38, 0x1C3A), (0x1C4A, 0x1C4C),
    (0x1C89, 0x1C8F), (0x1CBB, 0x1CBC), (0x1CC8, 0x1CCF), (0x1CFB, 0x1CFF), (0x1DFA, 0x1DFA), (0x1F16, 0x1F17),
    (0x1F1E, 0x1F1F), (0x1F46, 0x1F47), (0x1F4E, 0x1F4F), (0x1F58, 0x1F58), (0x1F5A, 0x1F5A), (0x1F5C, 0x1F5C),
    (0x1F5E, 0x1F5E), (0x1F7E, 0x1F7F), (0x1FB5, 0x1FB5), (0x1FC5, 0x1FC5), (0x1FD4, 0x1FD5), (0x1FDC, 0x1FDC),
    (0x1FF0, 0x1FF1), (0x1FF5, 0x1FF5), (0x1FFF, 0x200F), (0x2028, 0x202F), (0x205F, 0x206F), (0x2072, 0x2073),
    (0x208F, 0x208F), (0x209D, 0x209F), (0x20C0, 0x20CF), (0x20F1, 0x20FF), (0x218C, 0x218F), (0x2427, 0x243F),
    (0x244B, 0x245F), (0x2B74, 0x2B75), (0x2B96, 0x2B97), (0x2C2F, 0x2C2F), (0x2C5F, 0x2C5F), (0x2CF4, 0x2CF8),
    (0x2D26, 0x2D26), (0x2D28, 0x2D2C), (0x2D2E, 0x2D2F), (0x2D68, 0x2D6E), (0x2D71, 0x2D7E), (0x2D97, 0x2D9F),
    (0x2DA7, 0x2DA7), (0x2DAF, 0x2DAF), (0x2DB7, 0x2DB7), (0x2DBF, 0x2DBF), (0x2DC7, 0x2DC7), (0x2DCF, 0x2DCF),
    (0x2DD7, 0x2DD7), (0x2DDF, 0x2DDF), (0x2E50, 0x2E7F), (0x2E9A, 0x2E9A), (0x2EF4, 0x2EFF), (0x2FD6, 0x2FEF),
    (0x2FFC, 0x3000), (0x3040, 0x3040), (0x3097, 0x3098), (0x3100, 0x3104), (0x3130, 0x3130), (0x318F, 0x318F),
    (0x31BB, 0x31BF), (0x31E4, 0x31EF), (0x321F, 0x321F), (0x4DB6, 0x4DBF), (0x9FF0, 0x9FFF), (0xA48D, 0xA48F),
    (0xA4C7, 0xA4CF), (0xA62C, 0xA63F), (0xA6F8, 0xA6FF), (0xA7C0, 0xA7C1), (0xA7C7, 0xA7F6), (0xA82C, 0xA82F),
    (0xA83A, 0xA83F), (0xA878, 0xA87F), (0xA8C6, 0xA8CD), (0xA8DA, 0xA8DF), (0xA954, 0xA95E), (0xA97D, 0xA97F),
    (0xA9CE, 0xA9CE), (0xA9DA, 0xA9DD), (0xA9FF, 0xA9FF), (0xAA37, 0xAA3F), (0xAA4E, 0xAA4F), (0xAA5A, 0xAA5B),
    (0xAAC3, 0xAADA), (0xAAF7, 0xAB00), (0xAB07, 0xAB08), (0xAB0F, 0xAB10), (0xAB17, 0xAB1F), (0xAB27, 0xAB27),
    (0xAB2F, 0xAB2F), (0xAB68, 0xAB6F), (0xABEE, 0xABEF), (0xABFA, 0xABFF), (0xD7A4, 0xD7AF), (0xD7C7, 0xD7CA),
    (0xD7FC, 0xF8FF), (0xFA6E, 0xFA6F), (0xFADA, 0xFAFF), (0xFB07, 0xFB12), (0xFB18, 0xFB1C), (0xFB37, 0xFB37),
    (0xFB3D, 0xFB3D), (0xFB3F, 0xFB3F), (0xFB42, 0xFB42), (0xFB45, 0xFB45), (0xFBC2, 0xFBD2), (0xFD40, 0xFD4F),
    (0xFD90, 0xFD91), (0xFDC8, 0xFDEF), (0xFDFE, 0xFDFF), (0xFE1A, 0xFE1F), (0xFE53, 0xFE53), (0xFE67, 0xFE67),
    (0xFE6C, 0xFE6F), (0xFE75, 0xFE75), (0xFEFD, 0xFF00), (0xFFBF, 0xFFC1), (0xFFC8, 0xFFC9), (0xFFD0, 0xFFD1),
    (0xFFD8, 0xFFD9), (0xFFDD, 0xFFDF), (0xFFE7, 0xFFE7), (0xFFEF, 0xFFFB), (0xFFFE, 0xFFFF), (0x1000C, 0x1000C),
    (0x10027, 0x10027), (0x1003B, 0x1003B), (0x1003E, 0x1003E), (0x1004E, 0x1004F), (0x1005E, 0x1007F), (0x100FB, 0x100FF),
    (0x10103, 0x10106), (0x10134, 0x10136), (0x1018F, 0x1018F), (0x1019C, 0x1019F), (0x101A1, 0x101CF), (0x101FE, 0x1027F),
    (0x1029D, 0x1029F), (0x102D1, 0x102DF), (0x102FC, 0x102FF), (0x10324, 0x1032C), (0x1034B, 0x1034F), (0x1037B, 0x1037F),
    (0x1039E, 0x1039E), (0x103C4, 0x103C7), (0x103D6, 0x103FF), (0x1049E, 0x1049F), (0x104AA, 0x104AF), (0x104D4, 0x104D7),
    (0x104FC, 0x104FF), (0x10528, 0x1052F), (0x10564, 0x1056E), (0x10570, 0x105FF), (0x10737, 0x1073F), (0x10756, 0x1075F),
    (0x10768, 0x107FF), (0x10806, 0x10807), (0x10809, 0x10809), (0x10836, 0x10836), (0x10839, 0x1083B), (0x1083D, 0x1083E),
    (0x10856, 0x10856), (0x1089F, 0x108A6), (0x108B0, 0x108DF), (0x108F3, 0x108F3), (0x108F6, 0x108FA), (0x1091C, 0x1091E),
    (0x1093A, 0x1093E), (0x10940, 0x1097F), (0x109B8, 0x109BB), (0x109D0, 0x109D1), (0x10A04, 0x10A04), (0x10A07, 0x10A0B),
    (0x10A14, 0x10A14), (0x10A18, 0x10A18), (0x10A36, 0x10A37), (0x10A3B, 0x10A3E), (0x10A49, 0x10A4F), (0x10A59, 0x10A5F),
    (0x10AA0, 0x10ABF), (0x10AE7, 0x10AEA), (0x10AF7, 0x10AFF), (0x10B36, 0x10B38), (0x10B56, 0x10B57), (0x10B73, 0x10B77),
    (0x10B92, 0x10B98), (0x10B9D, 0x10BA8), (0x10BB0, 0x10BFF), (0x10C49, 0x10C7F), (0x10CB3, 0x10CBF), (0x10CF3, 0x10CF9),
    (0x10D28, 0x10D2F), (0x10D3A, 0x10E5F), (0x10E7F, 0x10EFF), (0x10F28, 0x10F2F), (0x10F5A, 0x10FDF), (0x10FF7, 0x10FFF),
    (0x1104E, 0x11051), (0x11070, 0x1107E), (0x110BD, 0x110BD), (0x110C2, 0x110CF), (0x110E9, 0x110EF), (0x110FA, 0x110FF),
    (0x11135, 0x11135), (0x11147, 0x1114F), (0x11177, 0x1117F), (0x111CE, 0x111CF), (0x111E0, 0x111E0), (0x111F5, 0x111FF),
    (0x11212, 0x11212), (0x1123F, 0x1127F), (0x11287, 0x11287), (0x11289, 0x11289), (0x1128E, 0x1128E), (0x1129E, 0x1129E),
    (0x112AA, 0x112AF), (0x112EB, 0x112EF), (0x112FA, 0x112FF), (0x11304, 0x11304), (0x1130D, 0x1130E), (0x11311, 0x11312),
    (0x11329, 0x11329), (0x11331, 0x11331), (0x11334, 0x11334), (0x1133A, 0x1133A), (0x11345, 0x11346), (0x11349, 0x1134A),
    (0x1134E, 0x1134F), (0x11351, 0x11356), (0x11358, 0x1135C), (0x11364, 0x11365), (0x1136D, 0x1136F), (0x11375, 0x113FF),
    (0x1145A, 0x1145A), (0x1145C, 0x1145C), (0x11460, 0x1147F), (0x114C8, 0x114CF), (0x114DA, 0x1157F), (0x115B6, 0x115B7),
    (0x115DE, 0x115FF), (0x11645, 0x1164F), (0x1165A, 0x1165F), (0x1166D, 0x1167F), (0x116B9, 0x116BF), (0x116CA, 0x116FF),
    (0x1171B, 0x1171C), (0x1172C, 0x1172F), (0x11740, 0x117FF), (0x1183C, 0x1189F), (0x118F3, 0x118FE), (0x11900, 0x1199F),
    (0x119A8, 0x119A9), (0x119D8, 0x119D9), (0x119E5, 0x119FF), (0x11A48, 0x11A4F), (0x11AA3, 0x11ABF), (0x11AF9, 0x11BFF),
    (0x11C09, 0x11C09), (0x11C37, 0x11C37), (0x11C46, 0x11C4F), (0x11C6D, 0x11C6F), (0x11C90, 0x11C91), (0x11CA8, 0x11CA8),
    (0x11CB7, 0x11CFF), (0x11D07, 0x11D07), (0x11D0A, 0x11D0A), (0x11D37, 0x11D39), (0x11D3B, 0x11D3B), (0x11D3E, 0x11D3E),
    (0x11D48, 0x11D4F), (0x11D5A, 0x11D5F), (0x11D66, 0x11D66), (0x11D69, 0x11D69), (0x11D8F, 0x11D8F), (0x11D92, 0x11D92),
    (0x11D99, 0x11D9F), (0x11DAA, 0x11EDF), (0x11EF9, 0x11FBF), (0x11FF2, 0x11FFE), (0x1239A, 0x123FF), (0x1246F, 0x1246F),
    (0x12475, 0x1247F), (0x12544, 0x12FFF), (0x1342F, 0x143FF), (0x14647, 0x167FF), (0x16A39, 0x16A3F), (0x16A5F, 0x16A5F),
    (0x16A6A, 0x16A6D), (0x16A70, 0x16ACF), (0x16AEE, 0x16AEF), (0x16AF6, 0x16AFF), (0x16B46, 0x16B4F), (0x16B5A, 0x16B5A),
    (0x16B62, 0x16B62), (0x16B78, 0x16B7C), (0x16B90, 0x16E3F), (0x16E9B, 0x16EFF), (0x16F4B, 0x16F4E), (0x16F88, 0x16F8E),
    (0x16FA0, 0x16FDF), (0x16FE4, 0x16FFF), (0x187F8, 0x187FF), (0x18AF3, 0x1AFFF), (0x1B11F, 0x1B14F), (0x1B153, 0x1B163),
    (0x1B168, 0x1B16F), (0x1B2FC, 0x1BBFF), (0x1BC6B, 0x1BC6F), (0x1BC7D, 0x1BC7F), (0x1BC89, 0x1BC8F), (0x1BC9A, 0x1BC9B),
    (0x1BCA0, 0x1CFFF), (0x1D0F6, 0x1D0FF), (0x1D127, 0x1D128), (0x1D173, 0x1D17A), (0x1D1E9, 0x1D1FF), (0x1D246, 0x1D2DF),
    (0x1D2F4, 0x1D2FF), (0x1D357, 0x1D35F), (0x1D379, 0x1D3FF), (0x1D455, 0x1D455), (0x1D49D, 0x1D49D), (0x1D4A0, 0x1D4A1),
    (0x1D4A3, 0x1D4A4), (0x1D4A7, 0x1D4A8), (0x1D4AD, 0x1D4AD), (0x1D4BA, 0x1D4BA), (0x1D4BC, 0x1D4BC), (0x1D4C4, 0x1D4C4),
    (0x1D506, 0x1D506), (0x1D50B, 0x1D50C), (0x1D515, 0x1D515), (0x1D51D, 0x1D51D), (0x1D53A, 0x1D53A), (0x1D53F, 0x1D53F),
    (0x1D545, 0x1D545), (0x1D547, 0x1D549), (0x1D551, 0x1D551), (0x1D6A6, 0x1D6A7), (0x1D7CC, 0x1D7CD), (0x1DA8C, 0x1DA9A),
    (0x1DAA0, 0x1DAA0), (0x1DAB0, 0x1DFFF), (0x1E007, 0x1E007), (0x1E019, 0x1E01A), (0x1E022, 0x1E022), (0x1E025, 0x1E025),
    (0x1E02B, 0x1E0FF), (0x1E12D, 0x1E12F), (0x1E13E, 0x1E13F), (0x1E14A, 0x1E14D), (0x1E150, 0x1E2BF), (0x1E2FA, 0x1E2FE),
    (0x1E300, 0x1E7FF), (0x1E8C5, 0x1E8C6), (0x1E8D7, 0x1E8FF), (0x1E94C, 0x1E94F), (0x1E95A, 0x1E95D), (0x1E960, 0x1EC70),
    (0x1ECB5, 0x1ED00), (0x1ED3E, 0x1EDFF), (0x1EE04, 0x1EE04), (0x1EE20, 0x1EE20), (0x1EE23, 0x1EE23), (0x1EE25, 0x1EE26),
    (0x1EE28, 0x1EE28), (0x1EE33, 0x1EE33), (0x1EE38, 0x1EE38), (0x1EE3A, 0x1EE3A), (0x1EE3C, 0x1EE41), (0x1EE43, 0x1EE46),
    (0x1EE48, 0x1EE48), (0x1EE4A, 0x1EE4A), (0x1EE4C, 0x1EE4C), (0x1EE50, 0x1EE50), (0x1EE53, 0x1EE53), (0x1EE55, 0x1EE56),
    (0x1EE58, 0x1EE58), (0x1EE5A, 0x1EE5A), (0x1EE5C, 0x1EE5C), (0x1EE5E, 0x1EE5E), (0x1EE60, 0x1EE60), (0x1EE63, 0x1EE63),
    (0x1EE65, 0x1EE66), (0x1EE6B, 0x1EE6B), (0x1EE73, 0x1EE73), (0x1EE78, 0x1EE78), (0x1EE7D, 0x1EE7D), (0x1EE7F, 0x1EE7F),
    (0x1EE8A, 0x1EE8A), (0x1EE9C, 0x1EEA0), (0x1EEA4, 0x1EEA4), (0x1EEAA, 0x1EEAA), (0x1EEBC, 0x1EEEF), (0x1EEF2, 0x1EFFF),
    (0x1F02C, 0x1F02F), (0x1F094, 0x1F09F), (0x1F0AF, 0x1F0B0), (0x1F0C0, 0x1F0C0), (0x1F0D0, 0x1F0D0), (0x1F0F6, 0x1F0FF),
    (0x1F10D, 0x1F10F), (0x1F16D, 0x1F16F), (0x1F1AD, 0x1F1E5), (0x1F203, 0x1F20F), (0x1F23C, 0x1F23F), (0x1F249, 0x1F24F),
    (0x1F252, 0x1F25F), (0x1F266, 0x1F2FF), (0x1F6D6, 0x1F6DF), (0x1F6ED, 0x1F6EF), (0x1F6FB, 0x1F6FF), (0x1F774, 0x1F77F),
    (0x1F7D9, 0x1F7DF), (0x1F7EC, 0x1F7FF), (0x1F80C, 0x1F80F), (0x1F848, 0x1F84F), (0x1F85A, 0x1F85F), (0x1F888, 0x1F88F),
    (0x1F8AE, 0x1F8FF), (0x1F90C, 0x1F90C), (0x1F972, 0x1F972), (0x1F977, 0x1F979), (0x1F9A3, 0x1F9A4), (0x1F9AB, 0x1F9AD),
    (0x1F9CB, 0x1F9CC), (0x1FA54, 0x1FA5F), (0x1FA6E, 0x1FA6F), (0x1FA74, 0x1FA77), (0x1FA7B, 0x1FA7F), (0x1FA83, 0x1FA8F),
    (0x1FA96, 0x1FFFF), (0x2A6D7, 0x2A6FF), (0x2B735, 0x2B73F), (0x2B81E, 0x2B81F), (0x2CEA2, 0x2CEAF), (0x2EBE1, 0x2F7FF),
    (0x2FA1E, 0xE00FF), (0xE01F0, 0x10FFFF),
];

/// First code point (digit zero) of every block of ten Unicode decimal digits (category Nd,
/// Unicode 12.1.0), which `int()` and `float()` accept like ASCII digits.
#[rustfmt::skip]
const ND_ZEROS: [u32; 63] = [
    0x30, 0x660, 0x6F0, 0x7C0, 0x966, 0x9E6, 0xA66, 0xAE6, 0xB66, 0xBE6,
    0xC66, 0xCE6, 0xD66, 0xDE6, 0xE50, 0xED0, 0xF20, 0x1040, 0x1090, 0x17E0,
    0x1810, 0x1946, 0x19D0, 0x1A80, 0x1A90, 0x1B50, 0x1BB0, 0x1C40, 0x1C50, 0xA620,
    0xA8D0, 0xA900, 0xA9D0, 0xA9F0, 0xAA50, 0xABF0, 0xFF10, 0x104A0, 0x10D30, 0x11066,
    0x110F0, 0x11136, 0x111D0, 0x112F0, 0x11450, 0x114D0, 0x11650, 0x116C0, 0x11730, 0x118E0,
    0x11C50, 0x11D50, 0x11DA0, 0x16A60, 0x16B50, 0x1D7CE, 0x1D7D8, 0x1D7E2, 0x1D7EC, 0x1D7F6,
    0x1E140, 0x1E2F0, 0x1E950,
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn whitespace_like_python() {
        // Python: 'a b\x1c c\x85d e'.split() == ['a', 'b', 'c', 'd', 'e']
        assert_eq!(py_split_ws("a b\u{1c} c\u{85}d\u{2028}e"), vec!["a", "b", "c", "d", "e"]);
        // Python: ' x\x1c\x85　'.rstrip() == ' x'
        assert_eq!(py_rstrip(" x\u{1c}\u{85}\u{3000}"), " x");
        assert_eq!(py_split_ws("  "), Vec::<&str>::new());
        assert_eq!(py_split_ws("chr1\t100\t101\n"), vec!["chr1", "100", "101"]);
    }

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
    fn whitespace_set_is_pythons() {
        // CPython 3.8: [c for c in range(0x110000) if chr(c).isspace()]
        let mut exp: Vec<u32> = vec![0x9, 0xa, 0xb, 0xc, 0xd, 0x1c, 0x1d, 0x1e, 0x1f, 0x20, 0x85, 0xa0, 0x1680];
        exp.extend(0x2000..=0x200a);
        exp.extend([0x2028, 0x2029, 0x202f, 0x205f, 0x3000]);
        let got: Vec<u32> = (0..=0x10ffffu32).filter_map(char::from_u32).filter(|&c| is_py_space(c)).map(|c| c as u32).collect();
        assert_eq!(got, exp);
    }

    #[test]
    fn int_like_python() {
        // values observed with CPython 3.8 int(s)
        assert_eq!(py_int("100").unwrap(), 100);
        assert_eq!(py_int(" 7 \n").unwrap(), 7);
        assert_eq!(py_int("\t5\n").unwrap(), 5);
        assert_eq!(py_int("\x0b5\x0c").unwrap(), 5);
        assert_eq!(py_int("+5").unwrap(), 5);
        assert_eq!(py_int("-12").unwrap(), -12);
        assert_eq!(py_int("007").unwrap(), 7);
        assert_eq!(py_int("0_7").unwrap(), 7);
        assert_eq!(py_int("1_000").unwrap(), 1000);
        assert_eq!(py_int("\u{85}5").unwrap(), 5);
        assert_eq!(py_int("5\u{a0}").unwrap(), 5);
        assert_eq!(py_int("\u{3000}5").unwrap(), 5);
        assert_eq!(py_int("\u{ff15}").unwrap(), 5);
        assert_eq!(py_int("+\u{ff15}").unwrap(), 5);
        assert_eq!(py_int("\u{661}\u{662}").unwrap(), 12);
        assert_eq!(py_int("\u{1d7d8}\u{1d7e1}").unwrap(), 9); // mathematical double-struck 0 and 9
        assert_eq!(py_int("-9223372036854775808").unwrap(), i64::MIN);
        for bad in [
            "", " ", "1__0", "_1", "1_", "5.0", "0x10", "- 5", "+-5", "abc", "\u{1c}5", "5\u{1c}", "5\u{7f}", "1\u{200b}", "\u{b2}",
            "9223372036854775808",
        ] {
            assert!(py_int(bad).is_err(), "{:?}", bad);
        }
    }

    #[test]
    fn float_like_python() {
        // values observed with CPython 3.8 float(s)
        assert_eq!(py_float("0.5\n").unwrap(), 0.5);
        assert_eq!(py_float("1_0.5").unwrap(), 10.5);
        assert_eq!(py_float("1e5_0").unwrap(), 1e50);
        assert_eq!(py_float("+.5").unwrap(), 0.5);
        assert_eq!(py_float("5.").unwrap(), 5.0);
        assert_eq!(py_float("\u{ff11}.\u{ff15}").unwrap(), 1.5);
        assert_eq!(py_float("1\u{a0}").unwrap(), 1.0);
        assert_eq!(py_float(" 1e-3 ").unwrap(), 0.001);
        assert_eq!(py_float("-iNF").unwrap(), f64::NEG_INFINITY);
        assert!(py_float("infinity").unwrap().is_infinite());
        assert!(py_float("+nan").unwrap().is_nan());
        for bad in ["", "0.5x", "1._5", "1_e5", "1e_5", "1e", "e5", "0x1p3", "--1", "in f", "\u{1c}5", "5\u{1c}", "1\u{200b}"] {
            assert!(py_float(bad).is_err(), "{:?}", bad);
        }
        assert!(py_fdiv(1.0, 0.0).is_err());
        assert!(py_fdiv(0.0, -0.0).is_err());
    }

    #[test]
    fn repr_like_python() {
        // CPython 3.8: repr('\xad​\U000e0001͸a\xa0\x85 ')
        assert_eq!(
            py_repr_str("\u{ad}\u{200b}\u{e000}\u{e0001}\u{378}a\u{a0}\u{85}\u{2028}"),
            "'\\xad\\u200b\\ue000\\U000e0001\\u0378a\\xa0\\x85\\u2028'"
        );
        assert_eq!(py_repr_str("abc"), "'abc'");
        assert_eq!(py_repr_str("a'b"), "\"a'b\"");
        assert_eq!(py_repr_str("a'b\""), "'a\\'b\"'");
        assert_eq!(py_repr_str("x\ty\n\r\\"), "'x\\ty\\n\\r\\\\'");
        assert_eq!(py_repr_str("\x01\x7f é€😀"), "'\\x01\\x7f é€😀'");
        assert!(is_py_printable(' ') && !is_py_printable('\u{3000}') && is_py_printable('\u{3001}'));
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

    #[test]
    fn universal_newlines() {
        // Python: [l for l in open(f)] for b"a\tb\r\nc\rd\n\r\re\r"
        assert_eq!(split_universal_newlines("a\tb\r\nc\rd\n\r\re\r"), vec!["a\tb\n", "c\n", "d\n", "\n", "\n", "e\n"]);
        assert_eq!(split_universal_newlines("x\ny"), vec!["x\n", "y"]);
        assert_eq!(split_universal_newlines(""), Vec::<String>::new());
    }
}
