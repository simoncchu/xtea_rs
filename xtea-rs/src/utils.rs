//! Small helpers replacing Python standard-library behaviour (os.path, subprocess, str.split).

use anyhow::{anyhow, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::process::{Command, Stdio};

/// `CMD_RUNNER.run_cmd_small_output`: print and run a shell command, capturing stdout.
pub fn run_cmd_small_output(cmd: &str) {
    println!("Running command: {}\n", cmd);
    let _ = Command::new("sh").arg("-c").arg(cmd).stdout(Stdio::piped()).output();
}

/// `CMD_RUNNER.run_cmd_to_file`: run a shell command writing stdout into `sf_out`.
pub fn run_cmd_to_file(cmd: &str, sf_out: &str) {
    println!("Running command with output: {}\n", cmd);
    if let Ok(f) = File::create(sf_out) {
        let _ = Command::new("sh").arg("-c").arg(cmd).stdout(Stdio::from(f)).status();
    }
    let sf_err = format!("{}.err", sf_out);
    if Path::new(&sf_err).is_file() {
        let _ = std::fs::remove_file(&sf_err);
    }
}

/// `Popen(cmd, shell=True, stdout=PIPE).communicate()` without printing.
pub fn run_cmd_quiet(cmd: &str) -> String {
    match Command::new("sh").arg("-c").arg(cmd).stdout(Stdio::piped()).output() {
        Ok(o) => String::from_utf8_lossy(&o.stdout).into_owned(),
        Err(_) => String::new(),
    }
}

/// `print(cmd); Popen(cmd, shell=True, stdout=PIPE).communicate()` (pipeline generators).
pub fn run_cmd_print(cmd: &str) {
    println!("{}", cmd);
    let _ = Command::new("sh").arg("-c").arg(cmd).stdout(Stdio::piped()).output();
}

/// os.path.isfile
pub fn is_file(p: &str) -> bool {
    Path::new(p).is_file()
}

/// os.path.exists
pub fn exists(p: &str) -> bool {
    Path::new(p).exists()
}

/// os.path.islink
pub fn is_link(p: &str) -> bool {
    std::fs::symlink_metadata(p).map(|m| m.file_type().is_symlink()).unwrap_or(false)
}

/// os.remove (ignoring errors)
pub fn remove_file(p: &str) {
    let _ = std::fs::remove_file(p);
}

/// os.path.basename
pub fn basename(p: &str) -> String {
    match p.rfind('/') {
        Some(i) => p[i + 1..].to_string(),
        None => p.to_string(),
    }
}

/// Append "/" when missing.
pub fn with_slash(p: &str) -> String {
    if p.ends_with('/') {
        p.to_string()
    } else {
        format!("{}/", p)
    }
}

/// `mkdir -p`
pub fn mkdir_p(p: &str) {
    let _ = std::fs::create_dir_all(p);
}

/// Python `str.split()` (split on runs of whitespace, dropping empty fields).
#[inline]
pub fn split_ws(s: &str) -> Vec<&str> {
    crate::py::text::py_split_ws(s)
}

/// Lines of a text file without their line terminator (Python iteration over `open(path)`:
/// strict UTF-8, universal newlines, then the trailing "\n" removed).
pub fn read_lines(path: &str) -> Result<Vec<String>> {
    let mut v = crate::py::text::read_text_lines(path)?;
    for l in v.iter_mut() {
        if l.ends_with('\n') {
            l.pop();
        }
    }
    Ok(v)
}

/// Iterate over the raw lines of a file, each including its trailing "\n" if present (Python's
/// iteration over `open(path)`: strict UTF-8 and universal newlines).
pub fn read_raw_lines(path: &str) -> Result<Vec<String>> {
    crate::py::text::read_text_lines(path)
}

/// Python `line.rstrip()` (strip all trailing whitespace).
#[inline]
pub fn rstrip(s: &str) -> &str {
    crate::py::text::py_rstrip(s)
}

/// Buffered writer creation (Python `open(path, "w")`).
pub fn create(path: &str) -> Result<BufWriter<File>> {
    let f = File::create(path).map_err(|e| anyhow!("cannot create {}: {}", path, e))?;
    Ok(BufWriter::new(f))
}

/// Buffered appender (Python `open(path, "a")`).
pub fn append(path: &str) -> Result<BufWriter<File>> {
    let f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .map_err(|e| anyhow!("cannot open {}: {}", path, e))?;
    Ok(BufWriter::new(f))
}

/// Write a whole string to a file.
pub fn write_string(path: &str, s: &str) -> Result<()> {
    let mut w = create(path)?;
    w.write_all(s.as_bytes())?;
    w.flush()?;
    Ok(())
}

/// Python's `str(int)` for usize/i64 is the plain decimal; for floats use py::fmt::float_repr.
pub use crate::py::fmt::float_repr as fstr;

/// Parse an int field like Python `int(s)`; returns error with context.
pub fn pint(s: &str) -> Result<i64> {
    crate::py::fmt::parse_py_int(s).ok_or_else(|| anyhow!("invalid literal for int() with base 10: '{}'", s))
}

/// Parse a float field like Python `float(s)`.
pub fn pfloat(s: &str) -> Result<f64> {
    crate::py::fmt::parse_py_float(s).ok_or_else(|| anyhow!("could not convert string to float: '{}'", s))
}

/// Python `int(float(s))`-style parse used when files contain "12.0" in integer columns.
pub fn pint_lenient(s: &str) -> Result<i64> {
    if let Some(v) = crate::py::fmt::parse_py_int(s) {
        return Ok(v);
    }
    Err(anyhow!("invalid literal for int() with base 10: '{}'", s))
}

/// Create a symbolic link like `ln -s target link` (link may be a directory: then the link is
/// created inside it with the target's basename).
pub fn ln_s(target: &str, link: &str) {
    let lp = Path::new(link);
    let dest = if lp.is_dir() {
        lp.join(basename(target))
    } else {
        lp.to_path_buf()
    };
    let _ = std::os::unix::fs::symlink(target, dest);
}

/// Copy a file (shutil.copyfile).
pub fn copyfile(src: &str, dst: &str) -> Result<()> {
    std::fs::copy(src, dst).map_err(|e| anyhow!("copy {} -> {} failed: {}", src, dst, e))?;
    Ok(())
}

/// Python `max(list)` semantics for picking the index of the first maximal element.
pub fn first_max_index<T: PartialOrd>(v: &[T]) -> Option<usize> {
    if v.is_empty() {
        return None;
    }
    let mut best = 0;
    for i in 1..v.len() {
        if v[i] > v[best] {
            best = i;
        }
    }
    Some(best)
}

/// Run `f(i)` for i in 0..n using the global rayon pool, returning results in order
/// (multiprocessing.Pool.map replacement).
pub fn par_map<T: Send, F: Fn(usize) -> T + Sync + Send>(n: usize, f: F) -> Vec<T> {
    use rayon::prelude::*;
    (0..n).into_par_iter().map(|i| f(i)).collect()
}
