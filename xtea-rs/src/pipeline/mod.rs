//! Ports of xTea's pipeline-generator scripts:
//!
//! | Python                                            | Rust binary           | module            |
//! |---------------------------------------------------|-----------------------|-------------------|
//! | `gnrt_pipeline_local_v38.py` (= `bin/xtea`)       | `xtea`                | [`local_v38`]     |
//! | `gnrt_pipeline_local.py` (= `bin/xtea_hg19`)      | `xtea_hg19`           | [`local_hg19`]    |
//! | `gnrt_pipeline_local_chm13.py`                    | `xtea_chm13`          | [`local_chm13`]   |
//! | `gnrt_pipeline_germline_cloud_v38.py`             | `xtea_germline_cloud` | [`germline_cloud`]|
//!
//! Generated files and printed output are byte-identical to the Python scripts, with one
//! deliberate change of the invocation convention: the generated commands call `${XTEA_MAIN}`
//! (the Rust `xtea_main` / `xtea_main_germline_cloud` binaries) instead of
//! `python ${XTEA_PATH}"x_TEA_main*.py"`, and every `<rep>.config` gets one extra line
//! `XTEA_MAIN <path>` right after the `XTEA_PATH` line, so that it reaches the parameter block
//! of each generated script as `XTEA_MAIN=<path>` (see [`xtea_main_path`]).
//!
//! This module holds what the generators share: Python-compatible helpers (text files, `str.split`,
//! `Popen(..., shell=True)`, optparse command lines), the configuration/script writers, the
//! driver of the three local generators (`gnrt_running_shell`, `split_bam_list`,
//! `gnrt_pipelines`, the `__main__` block) and the command builders whose text is identical in
//! several scripts.  The script-specific command strings live in the per-generator modules.
//!
//! Python quirks that are reproduced on purpose (see the generator modules for more):
//! * `split_bam_list` decides whether to read the 10X list from the module-global `sf_bams_10X`
//!   instead of its own parameter;
//! * HERV configurations of the hg38/CHM13/cloud scripts are written with `REF` before
//!   `ANNOTATION1` (argument order of the `write_to_config` call);
//! * `load_par_config` keeps only the first two whitespace-separated fields of a line;
//! * options left at `None` are formatted as `None` (e.g. `#SBATCH -p None`, `tar -zxvf None`)
//!   or raise where Python raises (string concatenation, `os.path.isfile(None)`, `None & 1`).

pub mod germline_cloud;
pub mod local_chm13;
pub mod local_hg19;
pub mod local_v38;

use crate::py::fmt::{float_repr, py_str_repr};
use anyhow::{anyhow, Result};
use indexmap::{IndexMap, IndexSet};
use std::ffi::OsString;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

pub const S_VERSION: &str = "v0.1";
pub const ILLUMINA: &str = "illumina";
pub const X10: &str = "10X";
pub const PUB_CLIP: &str = "pub_clip";
pub const CASE_LIST_SUFFIX: &str = ".case_list";
pub const CTRL_LIST_SUFFIX: &str = ".ctrl_list";

/// Name of the Rust port of `x_TEA_main.py` (called by the local generators).
pub const XTEA_MAIN_EXE: &str = "xtea_main";
/// Name of the Rust port of `x_TEA_main_germline_cloud.py` (called by the cloud generator).
pub const XTEA_MAIN_GERMLINE_CLOUD_EXE: &str = "xtea_main_germline_cloud";

// ------------------------------------------------------------------------------------------
// Python error emulation
// ------------------------------------------------------------------------------------------

fn type_error(msg: &str) -> anyhow::Error {
    anyhow!("TypeError: {}", msg)
}

fn index_error(msg: &str) -> anyhow::Error {
    anyhow!("IndexError: {}", msg)
}

/// `str + None`
fn concat_none_error() -> anyhow::Error {
    type_error("can only concatenate str (not \"NoneType\") to str")
}

/// `iflag & 1` with `iflag` None (the `-f` option has no default in most scripts).
pub fn need_flag(iflag: Option<i64>) -> Result<i64> {
    iflag.ok_or_else(|| type_error("unsupported operand type(s) for &: 'NoneType' and 'int'"))
}

/// Value of a string option used in a concatenation (`"REF " + sf_ref + "\n"`).
fn need_str(o: Option<&str>) -> Result<&str> {
    o.ok_or_else(concat_none_error)
}

/// `if s[-1] != "/": s += "/"` (TypeError for None, IndexError for "").
pub fn py_add_slash(o: Option<&str>) -> Result<String> {
    let s = o.ok_or_else(|| type_error("'NoneType' object is not subscriptable"))?;
    let last = s.chars().last().ok_or_else(|| index_error("string index out of range"))?;
    Ok(if last != '/' { format!("{}/", s) } else { s.to_string() })
}

/// `fields[i]` of a `str.split()` result.
fn field<'a>(fields: &[&'a str], i: usize) -> Result<&'a str> {
    fields.get(i).copied().ok_or_else(|| index_error("list index out of range"))
}

/// Python's message for an `OSError` raised by `open()`.
fn os_error(e: &std::io::Error, path: &str) -> anyhow::Error {
    match e.raw_os_error() {
        Some(code) => {
            // SAFETY: strerror returns a pointer to a NUL-terminated static/thread-local string.
            let msg = unsafe { std::ffi::CStr::from_ptr(libc::strerror(code)) }.to_string_lossy().into_owned();
            let class = match code {
                libc::ENOENT => "FileNotFoundError",
                libc::EACCES | libc::EPERM => "PermissionError",
                libc::EISDIR => "IsADirectoryError",
                libc::ENOTDIR => "NotADirectoryError",
                _ => "OSError",
            };
            anyhow!("{}: [Errno {}] {}: {}", class, code, msg, py_str_repr(path))
        }
        None => anyhow!("OSError: {}: {}", e, py_str_repr(path)),
    }
}

/// `str(x)` of an optional string value (`"{0}".format(None)` gives `None`).
pub fn py_str_or_none(o: Option<&str>) -> &str {
    o.unwrap_or("None")
}

/// `str(x)` of an optional int value.
pub fn py_int_or_none(o: Option<i64>) -> String {
    match o {
        Some(v) => v.to_string(),
        None => "None".to_string(),
    }
}

/// `repr(list_of_str)` (as printed by `print(l_sh)`).
pub fn py_list_repr(items: &[String]) -> String {
    let parts: Vec<String> = items.iter().map(|s| py_str_repr(s)).collect();
    format!("[{}]", parts.join(", "))
}

/// Python 3 true division of two ints (`int(a) / int(b)`): the correctly rounded quotient.
pub fn py_int_truediv(a: i64, b: i64) -> Result<f64> {
    if b == 0 {
        return Err(anyhow!("ZeroDivisionError: division by zero"));
    }
    let neg = (a < 0) != (b < 0);
    let (ua, ub) = (a.unsigned_abs() as u128, b.unsigned_abs() as u128);
    if ua == 0 {
        return Ok(if neg { -0.0 } else { 0.0 });
    }
    let bits = |x: u128| 128 - x.leading_zeros() as i32;
    // scale the dividend so that the integer quotient has at least 55 significant bits; the
    // remainder then only acts as a sticky bit for the final (round-half-even) conversion
    let shift = (55 + bits(ub) - bits(ua)).max(0);
    let num = ua << shift;
    let (q, r) = (num / ub, num % ub);
    let q = if r != 0 { q | 1 } else { q };
    let v = (q as f64) * 2f64.powi(-shift);
    Ok(if neg { -v } else { v })
}

// ------------------------------------------------------------------------------------------
// Python text I/O and string helpers
// ------------------------------------------------------------------------------------------

/// `str.isspace()` for one character: the separators of `str.split()` / `str.rstrip()`.
pub fn py_isspace(c: char) -> bool {
    matches!(
        c,
        '\t' | '\n' | '\x0b' | '\x0c' | '\r' | '\x1c' | '\x1d' | '\x1e' | '\x1f' | ' ' | '\u{85}' | '\u{a0}' | '\u{1680}' | '\u{2000}'
            ..='\u{200a}' | '\u{2028}' | '\u{2029}' | '\u{202f}' | '\u{205f}' | '\u{3000}'
    )
}

/// `str.split()` without arguments.
pub fn py_split(s: &str) -> Vec<&str> {
    s.split(py_isspace).filter(|t| !t.is_empty()).collect()
}

/// `str.rstrip()` without arguments.
pub fn py_rstrip(s: &str) -> &str {
    s.trim_end_matches(py_isspace)
}

/// Split text into lines the way iterating a text-mode file does (universal newlines: `\r\n`
/// and `\r` are read as `\n`; every line keeps its `\n`, the last one may have none).
pub fn split_universal_lines(text: &str) -> Vec<String> {
    let mut lines = Vec::new();
    let mut cur = String::new();
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '\r' | '\n' => {
                if c == '\r' && chars.peek() == Some(&'\n') {
                    chars.next();
                }
                cur.push('\n');
                lines.push(std::mem::take(&mut cur));
            }
            _ => cur.push(c),
        }
    }
    if !cur.is_empty() {
        lines.push(cur);
    }
    lines
}

/// Decode file content as UTF-8 (strict, like Python's default text decoding).
fn decode_utf8(data: Vec<u8>) -> Result<String> {
    String::from_utf8(data).map_err(|e| {
        let pos = e.utf8_error().valid_up_to();
        let byte = e.as_bytes()[pos];
        anyhow!("UnicodeDecodeError: 'utf-8' codec can't decode byte 0x{:02x} in position {}", byte, pos)
    })
}

/// `for line in open(path)`: all lines of a text file.
pub fn read_py_lines(path: &str) -> Result<Vec<String>> {
    let data = std::fs::read(path).map_err(|e| os_error(&e, path))?;
    Ok(split_universal_lines(&decode_utf8(data)?))
}

/// `with open(path, "w") as f: f.write(content)`
pub fn write_py_file(path: &str, content: &str) -> Result<()> {
    std::fs::write(path, content).map_err(|e| os_error(&e, path))
}

/// `os.path.exists` (follows symlinks; false on any error).
pub fn exists(p: &str) -> bool {
    Path::new(p).exists()
}

/// `os.path.isfile` (follows symlinks; false on any error).
pub fn is_file(p: &str) -> bool {
    Path::new(p).is_file()
}

/// Python `print(s)`.
pub fn py_print(s: &str) {
    let out = std::io::stdout();
    let mut lock = out.lock();
    let _ = lock.write_all(s.as_bytes());
    let _ = lock.write_all(b"\n");
}

// ------------------------------------------------------------------------------------------
// subprocess
// ------------------------------------------------------------------------------------------

/// `Popen(cmd, shell=True, stdout=PIPE).communicate()`: run the command with `/bin/sh -c`,
/// discard its standard output (stdin and stderr are inherited), ignore the exit status.
pub fn popen(cmd: &str) -> Result<()> {
    let _ = std::io::stdout().flush();
    Command::new("/bin/sh")
        .arg("-c")
        .arg(cmd)
        .stdin(Stdio::inherit())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .output()
        .map_err(|e| anyhow!("cannot run /bin/sh: {}", e))?;
    Ok(())
}

/// `run_cmd`: `print(cmd)` followed by `Popen(cmd, shell=True, stdout=PIPE).communicate()`.
pub fn run_cmd(cmd: &str) -> Result<()> {
    py_print(cmd);
    popen(cmd)
}

/// `cp_file`: copy with `cp` when the source is a file.
pub fn cp_file(sf_from: &str, sf_to: &str) -> Result<()> {
    let cmd = format!("cp {} {}", sf_from, sf_to);
    if !is_file(sf_from) {
        return Ok(());
    }
    run_cmd(&cmd)
}

/// `decompress`: `tar -zxvf <tar> -C <folder>`.
pub fn decompress(sf_in_tar: &str, sf_out: &str) -> Result<()> {
    run_cmd(&format!("tar -zxvf {} -C {}", sf_in_tar, sf_out))
}

// ------------------------------------------------------------------------------------------
// Command line (optparse emulation on top of clap)
// ------------------------------------------------------------------------------------------

/// The long option (name, takes a value) and the short options of a clap command, in definition
/// order.  The generators define their own `-h/--help` flag (clap's is disabled) so that clap's
/// prefix matching sees it, as optparse's does (`--h` is ambiguous with `--hard`).
fn option_specs(cmd: &clap::Command) -> (Vec<(String, bool)>, Vec<(char, usize)>) {
    let mut cmd = cmd.clone();
    cmd.build();
    let mut specs = Vec::new();
    let mut shorts = Vec::new();
    for a in cmd.get_arguments() {
        if a.is_positional() {
            continue;
        }
        if let Some(long) = a.get_long() {
            if let Some(s) = a.get_short() {
                shorts.push((s, specs.len()));
            }
            specs.push((long.to_string(), a.get_action().takes_values()));
        }
    }
    (specs, shorts)
}

/// optparse's `_match_abbrev`: exact match, else the unique option starting with `name`.
fn resolve_long<'a>(specs: &'a [(String, bool)], name: &str) -> Option<&'a (String, bool)> {
    if let Some(s) = specs.iter().find(|s| s.0 == name) {
        return Some(s);
    }
    let mut it = specs.iter().filter(|s| s.0.starts_with(name));
    match (it.next(), it.next()) {
        (Some(s), None) => Some(s),
        _ => None,
    }
}

/// Result of [`optparse_argv`].
pub struct OptparseArgv {
    /// the command line for clap
    pub argv: Vec<OsString>,
    /// `-h/--help` was reached: optparse prints the help and exits at that point, ignoring the
    /// rest of the command line (`argv` then holds the arguments before it)
    pub help: bool,
}

/// Re-spell a command line the way optparse reads it, so that clap produces the values optparse
/// would have produced.  optparse takes an option value from the rest of a short-option cluster
/// or from the next argument whatever it looks like (`-p=x` means `=x`, `-p --` means `--`),
/// while clap strips the `=` of `-p=x`; every option recognised here (long names resolved by
/// exact match or unique prefix, as optparse does) is therefore passed on in the unambiguous
/// `--long=value` / `--long` form.  From the first argument optparse would reject (unknown or
/// ambiguous option, missing value, value given to a flag) the rest is passed through unchanged
/// for clap to report.
pub fn optparse_argv(cmd: &clap::Command, argv: Vec<OsString>) -> OptparseArgv {
    let (specs, shorts) = option_specs(cmd);
    let mut it = argv.into_iter();
    let mut out: Vec<OsString> = it.next().into_iter().collect();
    let rest: Vec<OsString> = it.collect();
    let args: Vec<String> = match rest.iter().map(|a| a.to_str().map(str::to_string)).collect::<Option<Vec<_>>>() {
        Some(v) => v,
        None => {
            out.extend(rest);
            return OptparseArgv { argv: out, help: false };
        }
    };
    let done = |mut out: Vec<OsString>, from: &[String]| {
        out.extend(from.iter().map(OsString::from));
        OptparseArgv { argv: out, help: false }
    };
    let help = |out: Vec<OsString>| OptparseArgv { argv: out, help: true };
    let mut i = 0;
    while i < args.len() {
        let arg = &args[i];
        if arg == "--" {
            return done(out, &args[i..]);
        } else if let Some(body) = arg.strip_prefix("--") {
            let (name, explicit) = match body.find('=') {
                Some(p) => (&body[..p], Some(&body[p + 1..])),
                None => (body, None),
            };
            match resolve_long(&specs, name) {
                Some((long, true)) => {
                    let value = match explicit {
                        Some(v) => v.to_string(),
                        None if i + 1 < args.len() => {
                            i += 1;
                            args[i].clone()
                        }
                        None => return done(out, &args[i..]),
                    };
                    out.push(format!("--{}={}", long, value).into());
                }
                Some((long, false)) if explicit.is_none() => {
                    if long == "help" {
                        return help(out);
                    }
                    out.push(format!("--{}", long).into());
                }
                _ => return done(out, &args[i..]),
            }
        } else if arg.len() > 1 && arg.starts_with('-') {
            let chars: Vec<char> = arg[1..].chars().collect();
            let mut j = 0;
            while j < chars.len() {
                let (long, takes_value) = match shorts.iter().find(|(c, _)| *c == chars[j]) {
                    Some(&(_, idx)) => &specs[idx],
                    None => {
                        out.push(format!("-{}", chars[j..].iter().collect::<String>()).into());
                        return done(out, &args[i + 1..]);
                    }
                };
                if *takes_value {
                    let value = if j + 1 < chars.len() {
                        chars[j + 1..].iter().collect::<String>()
                    } else if i + 1 < args.len() {
                        i += 1;
                        args[i].clone()
                    } else {
                        out.push(format!("-{}", chars[j]).into());
                        return done(out, &[]);
                    };
                    out.push(format!("--{}={}", long, value).into());
                    break;
                }
                if long == "help" {
                    return help(out);
                }
                out.push(format!("--{}", long).into());
                j += 1;
            }
        } else {
            out.push(arg.into());
        }
        i += 1;
    }
    OptparseArgv { argv: out, help: false }
}

/// Parse the process arguments with optparse semantics (see [`optparse_argv`]): usage errors
/// exit with status 2; `-h/--help` prints optparse's help text ([`optparse_help`]) and exits
/// with status 0.
pub fn parse_opts<T: clap::Parser>(name: &'static str) -> T {
    let cmd = T::command().name(name).bin_name(name);
    let argv: Vec<OsString> = std::env::args_os().collect();
    let prog = argv
        .first()
        .map(|a| Path::new(a).file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_default())
        .unwrap_or_else(|| name.to_string());
    let parsed = optparse_argv(&cmd, argv);
    let matches = cmd.clone().try_get_matches_from(parsed.argv).unwrap_or_else(|e| e.exit());
    if parsed.help {
        let out = std::io::stdout();
        let mut lock = out.lock();
        let _ = lock.write_all(optparse_help(&cmd, &prog).as_bytes());
        let _ = lock.flush();
        std::process::exit(0);
    }
    T::from_arg_matches(&matches).unwrap_or_else(|e| e.exit())
}

/// The help text of optparse (`OptionParser.format_help` with the default
/// `IndentedHelpFormatter`, width taken from `$COLUMNS`), built from the clap definition: the
/// options in definition order, `-x METAVAR, --long=METAVAR` (the metavar is the value name,
/// i.e. `FILE` or the upper-cased destination), help texts wrapped with `textwrap.wrap`.
pub fn optparse_help(cmd: &clap::Command, prog: &str) -> String {
    let mut cmd = cmd.clone();
    cmd.build();
    let mut entries: Vec<(String, String)> = Vec::new();
    for a in cmd.get_arguments() {
        if a.is_positional() || a.is_hide_set() {
            continue;
        }
        let takes_value = a.get_action().takes_values();
        let metavar =
            a.get_value_names().and_then(|v| v.first()).map(|s| s.to_string()).unwrap_or_else(|| a.get_id().as_str().to_uppercase());
        let mut opts: Vec<String> = Vec::new();
        if let Some(s) = a.get_short() {
            opts.push(if takes_value { format!("-{} {}", s, metavar) } else { format!("-{}", s) });
        }
        if let Some(l) = a.get_long() {
            opts.push(if takes_value { format!("--{}={}", l, metavar) } else { format!("--{}", l) });
        }
        entries.push((opts.join(", "), a.get_help().map(|h| h.to_string()).unwrap_or_default()));
    }
    // HelpFormatter.__init__ (indent_increment=2, max_help_position=24)
    let width: i64 = std::env::var("COLUMNS").ok().and_then(|v| py_int_str(&v)).unwrap_or(80) - 2;
    let max_help_position = 24.min((width - 20).max(4));
    // store_option_strings (options are indented by 2)
    let max_len = entries.iter().map(|(o, _)| o.chars().count() as i64 + 2).max().unwrap_or(0);
    let help_position = (max_len + 2).min(max_help_position);
    let help_width = (width - help_position).max(11);
    let current_indent = 2;
    let mut s = format!("Usage: {} [options]\n\nOptions:\n", prog);
    for (opts, help) in &entries {
        // format_option
        let opt_width = help_position - current_indent - 2;
        let indent_first;
        if opts.chars().count() as i64 > opt_width {
            s += &format!("{:w$}{}\n", "", opts, w = current_indent as usize);
            indent_first = help_position;
        } else {
            s += &format!("{:w$}{:<ow$}  ", "", opts, w = current_indent as usize, ow = opt_width as usize);
            indent_first = 0;
        }
        if !help.is_empty() {
            let lines = py_textwrap(help, help_width as usize);
            let first = lines.first().map(String::as_str).unwrap_or("");
            s += &format!("{:w$}{}\n", "", first, w = indent_first as usize);
            for line in lines.iter().skip(1) {
                s += &format!("{:w$}{}\n", "", line, w = help_position as usize);
            }
        } else if !s.ends_with('\n') {
            s.push('\n');
        }
    }
    s
}

/// Python `int(str)` (surrounding whitespace, sign and single underscores between digits).
fn py_int_str(v: &str) -> Option<i64> {
    let t = v.trim_matches(py_isspace);
    let (neg, digits) = match t.strip_prefix('-') {
        Some(d) => (true, d),
        None => (false, t.strip_prefix('+').unwrap_or(t)),
    };
    if digits.is_empty() || digits.starts_with('_') || digits.ends_with('_') || digits.contains("__") {
        return None;
    }
    let n: i64 = digits.replace('_', "").parse().ok()?;
    Some(if neg { -n } else { n })
}

/// `textwrap.wrap(text, width)` of Python 3.8 with the default options (tabs expanded,
/// whitespace replaced, hyphenated words breakable, long words broken, whitespace dropped
/// at line ends).
pub fn py_textwrap(text: &str, width: usize) -> Vec<String> {
    let is_ws = |c: char| matches!(c, '\t' | '\n' | '\x0b' | '\x0c' | '\r' | ' ');
    // _munge_whitespace: expandtabs(8), then every whitespace character becomes a space
    let mut munged: Vec<char> = Vec::new();
    let mut col = 0;
    for c in text.chars() {
        match c {
            '\t' => {
                let n = 8 - col % 8;
                munged.extend(std::iter::repeat(' ').take(n));
                col += n;
            }
            '\n' | '\r' => {
                munged.push(' ');
                col = 0;
            }
            _ => {
                munged.push(if is_ws(c) { ' ' } else { c });
                col += 1;
            }
        }
    }
    // _split: wordsep_re.split
    let t = &munged;
    let n = t.len();
    let is_word = |c: char| c.is_alphanumeric() || c == '_';
    let is_letter = |c: char| is_word(c) && !c.is_numeric();
    let is_wp = |c: char| is_word(c) || matches!(c, '!' | '"' | '\'' | '&' | '.' | ',' | '?');
    // `-{2,}\w` at position p
    let em_dash_at = |p: usize| {
        let mut j = p;
        while j < n && t[j] == '-' {
            j += 1;
        }
        if j - p >= 2 && j < n && is_word(t[j]) {
            Some(j)
        } else {
            None
        }
    };
    let mut chunks: Vec<String> = Vec::new();
    let mut i = 0;
    while i < n {
        let end = if is_ws(t[i]) {
            let mut j = i;
            while j < n && is_ws(t[j]) {
                j += 1;
            }
            j
        } else if let Some(j) = (i > 0 && is_wp(t[i - 1])).then(|| em_dash_at(i)).flatten() {
            j
        } else {
            // word: shortest non-whitespace run followed by a breakable hyphen, the end of the
            // word or an em-dash
            let mut e = i + 1;
            loop {
                if e < n && t[e] == '-' {
                    let behind = (e >= 2 && is_letter(t[e - 2]) && is_letter(t[e - 1]))
                        || (e >= 3 && is_letter(t[e - 3]) && t[e - 2] == '-' && is_letter(t[e - 1]));
                    let ahead =
                        e + 2 < n && is_letter(t[e + 1]) && (is_letter(t[e + 2]) || (t[e + 2] == '-' && e + 3 < n && is_letter(t[e + 3])));
                    if behind && ahead {
                        break e + 1;
                    }
                }
                if e == n || is_ws(t[e]) {
                    break e;
                }
                if is_wp(t[e - 1]) && em_dash_at(e).is_some() {
                    break e;
                }
                e += 1;
            }
        };
        chunks.push(t[i..end].iter().collect());
        i = end;
    }
    // _wrap_chunks
    let blank = |c: &str| c.chars().all(py_isspace);
    let len = |c: &str| c.chars().count();
    chunks.reverse();
    let mut lines: Vec<String> = Vec::new();
    while !chunks.is_empty() {
        let mut cur_line: Vec<String> = Vec::new();
        let mut cur_len = 0;
        if blank(chunks.last().unwrap()) && !lines.is_empty() {
            chunks.pop();
        }
        while let Some(c) = chunks.last() {
            let l = len(c);
            if cur_len + l <= width {
                cur_len += l;
                cur_line.push(chunks.pop().unwrap());
            } else {
                break;
            }
        }
        if let Some(c) = chunks.last_mut() {
            if len(c) > width {
                // _handle_long_word
                let space_left = if width < 1 { 1 } else { width - cur_len };
                let head: String = c.chars().take(space_left).collect();
                let tail: String = c.chars().skip(space_left).collect();
                cur_line.push(head);
                *c = tail;
            }
        }
        if cur_line.last().map_or(false, |c| blank(c)) {
            cur_line.pop();
        }
        if !cur_line.is_empty() {
            lines.push(cur_line.concat());
        }
    }
    lines
}

/// Exit status of a generator run (Python exits with 1 on an uncaught exception).
pub fn finish(r: Result<()>) -> i32 {
    let _ = std::io::stdout().flush();
    match r {
        Ok(()) => 0,
        Err(e) => {
            eprintln!("Error: {:#}", e);
            1
        }
    }
}

/// Value of `XTEA_MAIN`: the absolute path of `exe_name` when it sits next to the running
/// generator executable, otherwise the bare name (looked up on `PATH` by the shell).
pub fn xtea_main_path(exe_name: &str) -> String {
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            let p = dir.join(exe_name);
            if p.is_file() {
                if let Some(s) = p.to_str() {
                    return s.to_string();
                }
            }
        }
    }
    exe_name.to_string()
}

// ------------------------------------------------------------------------------------------
// Repeat types and library configuration
// ------------------------------------------------------------------------------------------

/// One repeat type of a generator: its name (folder and `<name>.config`), its `-y` bit (also the
/// `--rtype` value returned by `get_flag_by_rep_type`) and its library files (relative to the
/// repeat library folder).
#[derive(Debug)]
pub struct RepType {
    pub name: &'static str,
    pub bit: i64,
    pub anno: &'static str,
    pub anno1: &'static str,
    pub copy_with_flank: &'static str,
    /// `None` writes `SF_FLANK null`
    pub flank: Option<&'static str>,
    pub cns: &'static str,
    /// The `write_to_config(sf_anno, sf_ref, sf_anno1, ...)` call of HERV in the hg38, CHM13 and
    /// cloud scripts (REF is written before ANNOTATION1).
    pub ref_before_anno1: bool,
}

/// Repeat types selected by `-y` (`if i_rep_type & bit != 0`), in bit order.
pub fn select_rep_types(table: &'static [RepType], i_rep_type: i64) -> Vec<&'static RepType> {
    table.iter().filter(|r| i_rep_type & r.bit != 0).collect()
}

/// `load_par_config`: (key, value) of every non-comment line.
pub fn load_par_config(sf_par_config: &str) -> Result<Vec<(String, String)>> {
    let mut l_pars = Vec::new();
    for line in read_py_lines(sf_par_config)? {
        if line.starts_with('#') {
            continue;
        }
        let fields = py_split(&line);
        l_pars.push((field(&fields, 0)?.to_string(), field(&fields, 1)?.to_string()));
    }
    Ok(l_pars)
}

/// `gnrt_parameters`: `KEY=value` lines.
pub fn gnrt_parameters(l_pars: &[(String, String)]) -> String {
    let mut s_pars = String::new();
    for (sid, svalue) in l_pars {
        s_pars.push_str(sid);
        s_pars.push('=');
        s_pars.push_str(svalue);
        s_pars.push('\n');
    }
    s_pars
}

/// `gnrt_lib_config`: write `<prefix>/<rep>.config` for every selected repeat type.  The only
/// addition to the original is the `XTEA_MAIN <path>` line written right after `XTEA_PATH`.
#[allow(clippy::too_many_arguments)]
pub fn gnrt_lib_config(
    l_rep_type: &[&RepType],
    sf_folder_rep: Option<&str>,
    sf_ref: Option<&str>,
    sf_gene: Option<&str>,
    sf_black_list: &str,
    sf_folder_xtea: Option<&str>,
    sf_config_prefix: &str,
    xtea_main: &str,
) -> Result<()> {
    let sf_folder_rep = py_add_slash(sf_folder_rep)?;
    let sf_folder_xtea = py_add_slash(sf_folder_xtea)?;
    let sf_config_prefix = py_add_slash(Some(sf_config_prefix))?;

    let s_bl = "BAM_LIST ${PREFIX}\"bam_list.txt\"\n";
    let s_bam1 = "BAM1 ${PREFIX}\"10X_phased_possorted_bam.bam\"\n";
    let s_bc_bam = "BARCODE_BAM ${PREFIX}\"10X_barcode_indexed.sorted.bam\"\n";
    let s_tmp = "TMP ${PREFIX}\"tmp/\"\n";
    let s_tmp_clip = "TMP_CLIP ${PREFIX}\"tmp/clip/\"\n";
    let s_tmp_cns = "TMP_CNS ${PREFIX}\"tmp/cns/\"\nTMP_TNSD ${PREFIX}\"tmp/transduction/\"\n";
    let sf_ref = format!("REF {}\n", need_str(sf_ref)?);
    let sf_xtea = format!("XTEA_PATH {}\n", sf_folder_xtea);
    let sf_xtea_main = format!("XTEA_MAIN {}\n", xtea_main);
    let sf_gene_anno = format!("GENE {}\n", need_str(sf_gene)?);
    let sf_black_list = format!("BLACK_LIST {}\n", sf_black_list);
    for rep in l_rep_type {
        let sf_config = format!("{}{}.config", sf_config_prefix, rep.name);
        let sf_anno = format!("ANNOTATION {}{}\n", sf_folder_rep, rep.anno);
        let sf_anno1 = format!("ANNOTATION1 {}{}\n", sf_folder_rep, rep.anno1);
        let sf_copy_with_flank = format!("L1_COPY_WITH_FLANK {}{}\n", sf_folder_rep, rep.copy_with_flank);
        let sf_flank = match rep.flank {
            Some(f) => format!("SF_FLANK {}{}\n", sf_folder_rep, f),
            None => "SF_FLANK null\n".to_string(),
        };
        let sf_cns = format!("L1_CNS {}{}\n", sf_folder_rep, rep.cns);
        let (second, third) = if rep.ref_before_anno1 { (&sf_ref, &sf_anno1) } else { (&sf_anno1, &sf_ref) };
        // write_to_config
        let content: String = [
            sf_anno.as_str(),
            second,
            third,
            &sf_gene_anno,
            &sf_black_list,
            &sf_copy_with_flank,
            &sf_flank,
            &sf_cns,
            &sf_xtea,
            &sf_xtea_main,
            s_bl,
            s_bam1,
            s_bc_bam,
            s_tmp,
            s_tmp_clip,
            s_tmp_cns,
        ]
        .concat();
        write_py_file(&sf_config, &content)?;
    }
    Ok(())
}

// ------------------------------------------------------------------------------------------
// Job script headers
// ------------------------------------------------------------------------------------------

/// `gnrt_script_head` (slurm).
pub fn gnrt_script_head(spartition: Option<&str>, ncores: i64, stime: &str, smemory: i64, s_id: &str) -> String {
    let mut s_head = String::from("#!/bin/bash\n\n");
    s_head += &format!("#SBATCH -c {}\n", ncores);
    s_head += &format!("#SBATCH -t {}\n", stime);
    s_head += &format!("#SBATCH --mem={}G\n", smemory);
    s_head += &format!("#SBATCH -p {}\n", py_str_or_none(spartition));
    s_head += &format!("#SBATCH -o {}_%j.out\n", s_id);
    s_head += "#SBATCH --mail-type=NONE\n";
    s_head += "#SBATCH --mail-user=chong.simonchu@gmail.com\n";
    if spartition == Some("park") || spartition == Some("priopark") {
        s_head += "#SBATCH --account=park_contrib\n\n";
    }
    s_head
}

/// `gnrt_script_head_lsf`; the memory per core is the Python 3 true division
/// `int(smemory)/int(ncores)` (a float, printed with `repr`).
pub fn gnrt_script_head_lsf(spartition: Option<&str>, ncores: i64, stime: &str, smemory: i64, s_id: &str) -> Result<String> {
    let mut s_head = String::from("#!/bin/bash\n\n");
    s_head += &format!("#BSUB -n {}\n", ncores);
    if !stime.is_empty() {
        s_head += &format!("#BSUB -W {}\n", stime);
    }
    s_head += &format!("#BSUB -M {}G\n", smemory);
    s_head += &format!("#BSUB -q {}\n", py_str_or_none(spartition));
    s_head += &format!("#BSUB -o {}_%J.out\n", s_id);
    s_head += &format!("#BSUB -e {}_%J.err\n", s_id);
    s_head += &format!("#BSUB -R \"rusage[mem={}G]\"\n", float_repr(py_int_truediv(smemory, ncores)?));
    Ok(s_head)
}

// ------------------------------------------------------------------------------------------
// Pipeline scripts
// ------------------------------------------------------------------------------------------

/// `gnrt_pipelines`: create the working folders of one repeat type, the 10X links and
/// `bam_list.txt`, and write `run_xTEA_pipeline.sh`; returns the scripts written.  The cloud
/// script's copy is the same function without the 10X part, which is never reached there
/// (its 10X list is always "null").
#[allow(clippy::too_many_arguments)]
pub fn gnrt_pipelines(
    s_head: &str,
    s_libs: &str,
    s_calling_cmd: &str,
    sf_id: &str,
    sf_bams: &str,
    sf_bams_10x: &str,
    sf_root_folder: &str,
    rep_type: &str,
) -> Result<Vec<String>> {
    let mut l_sbatch_files = Vec::new();
    let sf_working_folder = py_add_slash(Some(&format!("{}{}", sf_root_folder, rep_type)))?;

    let mut m_id: IndexSet<String> = IndexSet::new();
    for line in read_py_lines(sf_id)? {
        let sid = py_rstrip(&line).to_string();
        m_id.insert(sid);

        let sf_folder = &sf_working_folder;
        if !exists(sf_folder) {
            run_cmd(&format!("mkdir {}", sf_folder))?;
        }
        // create the temporary folders
        for sub in ["/tmp", "/tmp/clip", "/tmp/cns", "/tmp/cns/asm", "/tmp/transduction", "/tmp/igv"] {
            let sf_tmp = format!("{}{}", sf_folder, sub);
            if !exists(&sf_tmp) {
                run_cmd(&format!("mkdir {}", sf_tmp))?;
            }
        }
    }
    let mut m_bams: IndexMap<String, Vec<String>> = IndexMap::new();
    if sf_bams != "null" {
        for line in read_py_lines(sf_bams)? {
            let fields = py_split(&line);
            let sid = field(&fields, 0)?;
            let s_bam = field(&fields, 1)?;
            m_bams.entry(sid.to_string()).or_default().push(s_bam.to_string());
        }
    }

    let mut m_bams_10x: IndexMap<String, String> = IndexMap::new();
    if sf_bams_10x != "null" {
        for line in read_py_lines(sf_bams_10x)? {
            let fields = py_split(&line);
            let sid = field(&fields, 0)?;
            if !m_id.contains(sid) {
                continue;
            }
            let s_bam = field(&fields, 1)?;
            let s_barcode_bam = field(&fields, 2)?;
            m_bams_10x.insert(sid.to_string(), s_bam.to_string());

            // soft-link the bams and their indices
            let links = [
                (s_bam.to_string(), "10X_phased_possorted_bam.bam"),
                (s_barcode_bam.to_string(), "10X_barcode_indexed.sorted.bam"),
                (format!("{}.bai", s_bam), "10X_phased_possorted_bam.bam.bai"),
                (format!("{}.bai", s_barcode_bam), "10X_barcode_indexed.sorted.bam.bai"),
            ];
            for (target, name) in links {
                let sf_link = format!("{}{}", sf_working_folder, name);
                if !is_file(&sf_link) {
                    run_cmd(&format!("ln -s {} {}", target, sf_link))?;
                }
            }
        }
    }

    for sid in &m_id {
        let sf_folder = &sf_working_folder;
        if !exists(sf_folder) {
            continue;
        }
        // the bam list file
        let mut s_bam_list = String::new();
        if let Some(bams) = m_bams.get(sid) {
            for sf_tmp_bam in bams {
                s_bam_list += &format!("{}\t{}\n", sf_tmp_bam, ILLUMINA);
            }
        }
        if let Some(s_bam) = m_bams_10x.get(sid) {
            s_bam_list += &format!("{}\t{}\n", s_bam, X10);
        }
        write_py_file(&format!("{}bam_list.txt", sf_folder), &s_bam_list)?;
        // the pipeline file
        let sf_out_sh = format!("{}run_xTEA_pipeline.sh", sf_folder);
        let s_script =
            format!("{}PREFIX={}\n############\n############\n{}############\n############\n{}", s_head, sf_folder, s_libs, s_calling_cmd);
        write_py_file(&sf_out_sh, &s_script)?;
        l_sbatch_files.push(sf_out_sh);
    }
    Ok(l_sbatch_files)
}

// ------------------------------------------------------------------------------------------
// Command builders shared by several scripts
// ------------------------------------------------------------------------------------------

/// Arguments of the `gnrt_calling_command*` functions.
#[derive(Debug, Clone)]
pub struct CallArgs {
    pub iclip_c: i64,
    pub iclip_rp: i64,
    pub idisc_c: i64,
    pub iflt_clip: i64,
    pub iflt_disc: i64,
    /// `str(ncores)` (`None` in the cloud script when `-n` is not given)
    pub ncores: String,
    pub iflk_len: i64,
    pub min_tei_len: i64,
    pub iflag: Option<i64>,
    pub b_mosaic: bool,
    pub b_user_par: bool,
    pub b_force: bool,
    pub b_tumor: bool,
    pub b_hard_cut: bool,
    pub b_resume: bool,
    pub f_purity: f64,
    pub i_rep_type: i64,
    pub s_cfolder: String,
}

impl CallArgs {
    pub fn s_user(&self) -> &'static str {
        if self.b_user_par {
            "--user"
        } else {
            ""
        }
    }
    pub fn s_clean(&self) -> &'static str {
        if self.b_force {
            "--force"
        } else {
            ""
        }
    }
    pub fn s_tumor(&self) -> &'static str {
        if self.b_tumor {
            "--tumor"
        } else {
            ""
        }
    }
    /// `"--purity {0}".format(f_purity)` for tumor samples
    pub fn s_purity(&self) -> String {
        if self.b_tumor {
            format!("--purity {}", float_repr(self.f_purity))
        } else {
            String::new()
        }
    }
    pub fn s_hard_cut(&self) -> &'static str {
        if self.b_hard_cut {
            "--hard"
        } else {
            ""
        }
    }
    pub fn s_resume(&self) -> &'static str {
        if self.b_resume {
            "--resume"
        } else {
            ""
        }
    }
}

/// Step commands of `gnrt_calling_command` whose text is the same in the hg38, hg19, CHM13 and
/// cloud scripts (the cloud script only uses some of them and never sets `b_mosaic`).
pub struct RnaSteps {
    pub sclip_step: String,
    pub sdisc_step: String,
    pub sbarcode_step: String,
    pub sfilter_10x: String,
    pub s_filter: String,
    pub sf_trsdct: String,
    pub sf_sibling: String,
    pub sf_collect: String,
    pub sf_asm: String,
    pub sf_alg_ctg: String,
    pub sf_sr_asm: String,
}

/// The common part of `gnrt_calling_command`; `sr_asm_input` is the input of the `--sr_asm`
/// step for non-tumor samples (it differs in the CHM13 script).
pub fn rna_steps(a: &CallArgs, b_sva: bool, sr_asm_input: &str) -> RnaSteps {
    let (s_user, s_clean, s_tumor, s_purity, s_resume) = (a.s_user(), a.s_clean(), a.s_tumor(), a.s_purity(), a.s_resume());
    let (iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc) = (a.iclip_c, a.iclip_rp, a.idisc_c, a.iflt_clip, a.iflt_disc);
    let (ncores, iflk_len, i_rep_type, s_cfolder) = (&a.ncores, a.iflk_len, a.i_rep_type, &a.s_cfolder);

    let mut sclip_step = format!(
        concat!(
            "${{XTEA_MAIN}} -C -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  ",
            "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} ",
            "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5} {6} {7} {8} {9}\n"
        ),
        iclip_c, iclip_c, iclip_rp, ncores, s_cfolder, s_purity, s_user, s_clean, s_tumor, s_resume
    );
    if b_sva {
        sclip_step = format!(
            concat!(
                "${{XTEA_MAIN}} -C --sva -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  ",
                "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} ",
                "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5} {6} {7} {8} {9}\n"
            ),
            iclip_c, iclip_c, iclip_rp, ncores, s_cfolder, s_purity, s_user, s_clean, s_tumor, s_resume
        );
    }
    let mut sdisc_step = format!(
        concat!(
            "${{XTEA_MAIN}}  -D -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} ",
            "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} ",
            "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2} {3} {4} {5}\n"
        ),
        idisc_c, ncores, s_purity, s_user, s_tumor, s_resume
    );
    if b_sva {
        sdisc_step = format!(
            concat!(
                "${{XTEA_MAIN}}  -D --sva -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} ",
                "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} ",
                "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2} {3} {4} {5}\n"
            ),
            idisc_c, ncores, s_purity, s_user, s_tumor, s_resume
        );
    }
    if a.b_mosaic {
        sclip_step = format!(
            concat!(
                "${{XTEA_MAIN}} --mosaic -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  ",
                "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} ",
                "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5} {6}\n"
            ),
            iclip_c, iclip_c, iclip_rp, ncores, s_cfolder, s_user, s_clean
        );
        sdisc_step = format!(
            concat!(
                "${{XTEA_MAIN}}  -D --somatic -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} ",
                "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} ",
                "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n"
            ),
            idisc_c, ncores, s_user
        );
    }
    if a.b_mosaic && b_sva {
        sclip_step = format!(
            concat!(
                "${{XTEA_MAIN}} --mosaic --sva -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  ",
                "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} ",
                "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5} {6}\n"
            ),
            iclip_c, iclip_c, iclip_rp, ncores, s_cfolder, s_user, s_clean
        );
        sdisc_step = format!(
            concat!(
                "${{XTEA_MAIN}}  -D --sva --somatic -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} ",
                "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} ",
                "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n"
            ),
            idisc_c, ncores, s_user
        );
    }
    let sbarcode_step = format!(
        concat!(
            "${{XTEA_MAIN}} -B -i ${{PREFIX}}\"candidate_list_from_disc.txt\" --nb 400 ",
            "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM1}} -d ${{BARCODE_BAM}} -p ${{TMP}} ",
            "-o ${{PREFIX}}\"candidate_list_barcode.txt\" -n {0} {1}\n"
        ),
        ncores, s_user
    );
    let mut sfilter_10x = format!(
        concat!(
            "${{XTEA_MAIN}} -N --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} ",
            "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_barcode.txt\" ",
            "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} ",
            "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4} {5} {6} {7}\n"
        ),
        iflt_clip, iflt_disc, iflk_len, ncores, s_purity, s_user, s_tumor, s_resume
    );
    if b_sva {
        sfilter_10x = format!(
            concat!(
                "${{XTEA_MAIN}} -N --sva --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} ",
                "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_barcode.txt\" ",
                "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} ",
                "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4} {5} {6} {7}\n"
            ),
            iflt_clip, iflt_disc, iflk_len, ncores, s_purity, s_user, s_tumor, s_resume
        );
    }
    let mut s_filter = format!(
        concat!(
            "${{XTEA_MAIN}} -N --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} ",
            "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_from_disc.txt\" ",
            "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} ",
            "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4} {5} {6} {7}\n"
        ),
        iflt_clip, iflt_disc, iflk_len, ncores, s_purity, s_user, s_tumor, s_resume
    );
    if b_sva {
        s_filter = format!(
            concat!(
                "${{XTEA_MAIN}} -N --sva --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} ",
                "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_from_disc.txt\" ",
                "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} ",
                "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4} {5} {6} {7}\n"
            ),
            iflt_clip, iflt_disc, iflk_len, ncores, s_purity, s_user, s_tumor, s_resume
        );
    }
    let sf_trsdct = format!(
        concat!(
            "${{XTEA_MAIN}} --transduction --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_TNSD}} ",
            "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" ",
            "-r ${{L1_CNS}} --ref ${{REF}} --input2 ${{PREFIX}}\"candidate_list_from_disc.txt.clip_sites_raw_disc.txt\" ",
            "--rtype {4} -a ${{ANNOTATION1}} {5} {6} {7} ",
            "-o ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\"\n"
        ),
        iflt_clip, iflt_disc, iflk_len, ncores, i_rep_type, s_purity, s_tumor, s_resume
    );
    let mut sf_sibling = format!(
        concat!(
            "${{XTEA_MAIN}} --sibling --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_TNSD}} ",
            "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\" ",
            "-r ${{L1_CNS}} --ref ${{REF}} --input2 ${{PREFIX}}\"candidate_list_from_disc.txt.clip_sites_raw_disc.txt\" ",
            "--rtype {4} -a ${{ANNOTATION1}} --blacklist ${{BLACK_LIST}} {5} {6} {7} ",
            "-o ${{PREFIX}}\"candidate_sibling_transduction2.txt\"\n"
        ),
        iflt_clip, iflt_disc, iflk_len, ncores, i_rep_type, s_purity, s_tumor, s_resume
    );
    if b_sva {
        // here skip orphan transduction calling for SVA (as overlap with L1)
        sf_sibling = format!(
            concat!(
                "${{XTEA_MAIN}} --sibling --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_TNSD}} ",
                "--fflank \"\" --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\" ",
                "-r ${{L1_CNS}} --ref ${{REF}} --input2 ${{PREFIX}}\"candidate_list_from_disc.txt.clip_sites_raw_disc.txt\" ",
                "--rtype {4} -a ${{ANNOTATION1}} --blacklist ${{BLACK_LIST}} {5} {6} {7} ",
                "-o ${{PREFIX}}\"candidate_sibling_transduction2.txt\"\n"
            ),
            iflt_clip, iflt_disc, iflk_len, ncores, i_rep_type, s_purity, s_tumor, s_resume
        );
    }
    let sf_collect = format!(
        concat!(
            "${{XTEA_MAIN}} -E --nb 500 -b ${{BAM1}} -d ${{BARCODE_BAM}} --ref ${{REF}} ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\" -p ${{TMP}} -a ${{ANNOTATION}} -n {0} ",
            "--flklen {1}\n"
        ),
        ncores, iflk_len
    );
    let sf_asm = format!(
        concat!(
            "${{XTEA_MAIN}} -A -L -p ${{TMP}} --ref ${{REF}} -n {0} ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\"\n"
        ),
        ncores
    );
    let sf_alg_ctg = format!(
        concat!(
            "${{XTEA_MAIN}} -M -i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\" ",
            "-r ${{REF}} -n {0} -p ${{TMP}} --ref ${{L1_CNS}} ",
            "-o ${{PREFIX}}\"candidate_list_asm.txt\"\n"
        ),
        ncores
    );
    let mut sf_sr_asm = format!(
        concat!(
            "${{XTEA_MAIN}} --sr_asm -b ${{BAM_LIST}} -p ${{TMP_CNS}}\"asm\" ",
            " -n {0} -i ${{PREFIX}}\"{1}\" ",
            " --ref ${{REF}} -o ${{PREFIX}}\"assembled_ins_sequence.fa\" \n"
        ),
        ncores, sr_asm_input
    );
    if a.b_tumor {
        sf_sr_asm = format!(
            concat!(
                "${{XTEA_MAIN}} --sr_asm -b ${{BAM_LIST}} -p ${{TMP_CNS}}\"asm\" ",
                " -n {0} -i ${{PREFIX}}\"candidate_disc_filtered_cns_high_confident_post_filtering_somatic.txt\" ",
                " --ref ${{REF}} -o ${{PREFIX}}\"assembled_ins_sequence.fa\" \n"
            ),
            ncores
        );
    }
    RnaSteps {
        sclip_step,
        sdisc_step,
        sbarcode_step,
        sfilter_10x,
        s_filter,
        sf_trsdct,
        sf_sibling,
        sf_collect,
        sf_asm,
        sf_alg_ctg,
        sf_sr_asm,
    }
}

/// `--sr_asm` input of the hg38 and hg19 scripts.
pub const SR_ASM_INPUT_DEFAULT: &str = "candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt";

/// Post-filter steps (`sf_post_filter`, `sf_post_filter_hc`) of `gnrt_calling_command` in the
/// hg38, CHM13 and cloud scripts (the hg19 script has its own).
pub fn post_filter_steps_hg38(a: &CallArgs) -> (String, String) {
    let (i_rep_type, ncores, s_tumor, s_hard_cut) = (a.i_rep_type, &a.ncores, a.s_tumor(), a.s_hard_cut());
    let mut sf_post_filter = format!(
        concat!(
            "${{XTEA_MAIN}} --postF --rtype {0} -p ${{TMP_CNS}} ",
            "-n {1} -i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\" ",
            "-a ${{ANNOTATION1}} {2} {3} ",
            "-o ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\"\n"
        ),
        i_rep_type, ncores, s_tumor, s_hard_cut
    );
    if a.b_mosaic {
        sf_post_filter = format!(
            concat!(
                "${{XTEA_MAIN}} --postF --postFmosaic --rtype {0} -p ${{TMP_CNS}} ",
                "-n {1} -i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\" ",
                "-a ${{ANNOTATION1}} {2} ",
                "-o ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\"\n"
            ),
            i_rep_type, ncores, s_tumor
        );
    }
    let mut sf_post_filter_hc = format!(
        concat!(
            "${{XTEA_MAIN}} --postF --rtype {0} -p ${{TMP_CNS}} -n {1} ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt.high_confident\" -a ${{ANNOTATION1}} ",
            "--blacklist ${{BLACK_LIST}} {2} {3} ",
            "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt\"\n"
        ),
        i_rep_type, ncores, s_tumor, s_hard_cut
    );
    if a.b_mosaic {
        sf_post_filter_hc = format!(
            concat!(
                "${{XTEA_MAIN}} --postF --postFmosaic --rtype {0} -p ${{TMP_CNS}} -n {1} ",
                "-i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt.high_confident\" -a ${{ANNOTATION1}} {2} ",
                "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt\"\n"
            ),
            i_rep_type, ncores, s_tumor
        );
    }
    (sf_post_filter, sf_post_filter_hc)
}

/// Step commands of `gnrt_calling_command_non_RNA` / `gnrt_calling_command_MT` that are the
/// same in the hg38, hg19 and CHM13 scripts; `s_dna` is `--dna` (non-RNA) or `--dna --mit` (MT).
pub struct DnaSteps {
    pub sclip_step: String,
    pub sdisc_step: String,
    pub sbarcode_step: String,
    pub sfilter_10x: String,
    pub s_filter: String,
    pub sf_collect: String,
    pub sf_asm: String,
    pub sf_alg_ctg: String,
    pub sf_mutation: String,
}

pub fn dna_steps(a: &CallArgs, s_dna: &str) -> DnaSteps {
    let s_user = a.s_user();
    let (iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc) = (a.iclip_c, a.iclip_rp, a.idisc_c, a.iflt_clip, a.iflt_disc);
    let (ncores, iflk_len, min_tei_len, s_cfolder) = (&a.ncores, a.iflk_len, a.min_tei_len, &a.s_cfolder);
    let sclip_step = format!(
        concat!(
            "${{XTEA_MAIN}} -C {dna} -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  ",
            "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} ",
            "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5}\n"
        ),
        iclip_c,
        iclip_c,
        iclip_rp,
        ncores,
        s_cfolder,
        s_user,
        dna = s_dna
    );
    let mut sdisc_step = format!(
        concat!(
            "${{XTEA_MAIN}}  -D {dna} -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} ",
            "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} ",
            "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n"
        ),
        idisc_c,
        ncores,
        s_user,
        dna = s_dna
    );
    if a.b_mosaic {
        sdisc_step = format!(
            concat!(
                "${{XTEA_MAIN}}  -D --somatic {dna} -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} ",
                "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} ",
                "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n"
            ),
            idisc_c,
            ncores,
            s_user,
            dna = s_dna
        );
    }
    let sbarcode_step = format!(
        concat!(
            "${{XTEA_MAIN}} -B {dna} -i ${{PREFIX}}\"candidate_list_from_disc.txt\" --nb 400 ",
            "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM1}} -d ${{BARCODE_BAM}} -p ${{TMP}} ",
            "-o ${{PREFIX}}\"candidate_list_barcode.txt\" -n {0} {1}\n"
        ),
        ncores,
        s_user,
        dna = s_dna
    );
    let sfilter_10x = format!(
        concat!(
            "${{XTEA_MAIN}} -N {dna} --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} ",
            "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_barcode.txt\" ",
            "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} ",
            "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4}\n"
        ),
        iflt_clip,
        iflt_disc,
        iflk_len,
        ncores,
        s_user,
        dna = s_dna
    );
    let s_filter = format!(
        concat!(
            "${{XTEA_MAIN}} -N {dna} --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} ",
            "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_from_disc.txt\" ",
            "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} ",
            "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4}\n"
        ),
        iflt_clip,
        iflt_disc,
        iflk_len,
        ncores,
        s_user,
        dna = s_dna
    );
    let sf_collect = format!(
        concat!(
            "${{XTEA_MAIN}} -E {dna} --nb 500 -b ${{BAM1}} -d ${{BARCODE_BAM}} --ref ${{REF}} ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" -p ${{TMP}} -a ${{ANNOTATION}} -n {0} ",
            "--flklen {1}\n"
        ),
        ncores,
        iflk_len,
        dna = s_dna
    );
    let sf_asm = format!(
        concat!("${{XTEA_MAIN}} -A -L {dna} -p ${{TMP}} --ref ${{REF}} -n {0} ", "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\"\n"),
        ncores,
        dna = s_dna
    );
    let sf_alg_ctg = format!(
        concat!(
            "${{XTEA_MAIN}} -M {dna} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" ",
            "-r ${{REF}} -n {0} -p ${{TMP}} --ref ${{L1_CNS}} ",
            "-o ${{PREFIX}}\"candidate_list_asm.txt\"\n"
        ),
        ncores,
        dna = s_dna
    );
    let sf_mutation = format!(
        concat!(
            "${{XTEA_MAIN}} -I {dna} -p ${{TMP}} -n {0} ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" -r ${{L1_CNS}} ",
            "--teilen {1} -o ${{PREFIX}}\"internal_snp.vcf.gz\"\n"
        ),
        ncores,
        min_tei_len,
        dna = s_dna
    );
    DnaSteps { sclip_step, sdisc_step, sbarcode_step, sfilter_10x, s_filter, sf_collect, sf_asm, sf_alg_ctg, sf_mutation }
}

/// Steps of `gnrt_calling_command_annotation_vcf` / `..._annotation_vcf_tumor_case` whose text
/// is the same in the hg38, hg19 and CHM13 scripts (hg38/CHM13 append BamSnap steps).
pub struct AnnotationSteps {
    pub sf_gene: String,
    pub sf_gntp_classify: String,
    /// gVCF of the genotyped calls (`annotation_vcf`)
    pub sf_cvt_gvcf: String,
    /// gVCF of the calls with genes, not genotyped (`annotation_vcf_tumor_case`)
    pub sf_cvt_gvcf_tumor: String,
}

pub fn annotation_steps(ncores: &str, i_rep_type: i64) -> AnnotationSteps {
    let sf_gene = format!(
        concat!(
            "${{XTEA_MAIN}} --gene -a ${{GENE}} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt\" ",
            " -n {0} -o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt\"\n"
        ),
        ncores
    );
    let sf_gntp_classify = format!(
        concat!(
            "${{XTEA_MAIN}} --gntp_classify -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt\" ",
            " -n {0} --model ${{XTEA_PATH}}\"genotyping/DF21_model_1_2\" ",
            " -o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt\"\n"
        ),
        1
    );
    let sf_cvt_gvcf = format!(
        concat!(
            "${{XTEA_MAIN}} --gVCF -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt\" ",
            " -o ${{PREFIX}} -b ${{BAM_LIST}} --ref ${{REF}} --rtype {0}\n"
        ),
        i_rep_type
    );
    let sf_cvt_gvcf_tumor = format!(
        concat!(
            "${{XTEA_MAIN}} --gVCF -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt\" ",
            " -o ${{PREFIX}} -b ${{BAM_LIST}} --ref ${{REF}} --rtype {0}\n"
        ),
        i_rep_type
    );
    AnnotationSteps { sf_gene, sf_gntp_classify, sf_cvt_gvcf, sf_cvt_gvcf_tumor }
}

/// `gnrt_calling_command_annotation_vcf_somatic` (same in the hg38, hg19 and CHM13 scripts).
pub fn annotation_vcf_somatic(ncores: &str, i_rep_type: i64) -> String {
    let sf_gene = format!(
        concat!(
            "${{XTEA_MAIN}} --gene -a ${{GENE}} -i ${{PREFIX}}\"candidate_disc_filtered_cns_high_confident_post_filtering_somatic.txt\" ",
            " -n {0} -o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_somatic_with_gene.txt\"\n"
        ),
        ncores
    );
    let sf_cvt_gvcf = format!(
        concat!(
            "${{XTEA_MAIN}} --gVCF -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_somatic_with_gene.txt\" ",
            " -o ${{PREFIX}} -b ${{BAM_LIST}} --ref ${{REF}} --rtype {0}\n"
        ),
        i_rep_type
    );
    sf_gene + &sf_cvt_gvcf
}

/// `gnrt_calling_command_somatic_case_control`; the scripts only differ by the reference name
/// given to the igv steps (`igv_ref`: "hg38", "hg19" or "chm13").
pub fn somatic_case_control(a: &CallArgs, sf_ctr_list: &str, sf_case_control_bam_list: &str, igv_ref: &str) -> Result<String> {
    let (ncores, iflt_clip, iflt_disc, iflk_len) = (&a.ncores, a.iflt_clip, a.iflt_disc, a.iflk_len);
    let s_case_control = format!(
        concat!(
            "${{XTEA_MAIN}} --case_control -p ${{TMP_CNS}} -b {0} ",
            "-n {1} -r ${{L1_CNS}} -i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\" ",
            " --ref ${{REF}} --cr {2} --nd {3} --fflank ${{SF_FLANK}} --flklen {4} ",
            "-o ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering_somatic.txt\"\n"
        ),
        sf_ctr_list, ncores, iflt_clip, iflt_disc, iflk_len
    );
    let s_case_control_hc = format!(
        concat!(
            "${{XTEA_MAIN}} --case_control --somatic_hc -p ${{TMP_CNS}} ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt\" ",
            " --input2 ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering_somatic.txt\" -n {0} ",
            "--fflank ${{SF_FLANK}} --flklen {1} ",
            "-o ${{PREFIX}}\"candidate_disc_filtered_cns_high_confident_post_filtering_somatic.txt\"\n"
        ),
        ncores, iflk_len
    );
    let s_igv = format!(
        concat!(
            "${{XTEA_MAIN}} --igv --single_sample -p ${{PREFIX}}\"tmp/igv\" -b {0} ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\" ",
            " --ref \"{1}\" -e 500 ",
            "-o ${{PREFIX}}\"tmp/igv/igv_screenshot.txt\"\n"
        ),
        sf_case_control_bam_list, igv_ref
    );
    let s_igv_hc = format!(
        concat!(
            "${{XTEA_MAIN}} --igv --single_sample -p ${{PREFIX}}\"tmp/igv\" -b {0} ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns_high_confident_post_filtering_somatic.txt\" ",
            " --ref \"{1}\" -e 500 ",
            "-o ${{PREFIX}}\"tmp/igv/igv_screenshot_somatic.txt\"\n"
        ),
        sf_case_control_bam_list, igv_ref
    );
    // BamSnap
    let s_bamsnap = format!(
        concat!(
            "${{XTEA_MAIN}} --igv --bamsnap --single_sample -p ${{PREFIX}}\"tmp/igv\" -b {0} ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\" ",
            " --ref ${{REF}} -e 1000 -n {1} ",
            "-o ${{PREFIX}}\"tmp/igv/bamsnap_screenshot.txt\"\n"
        ),
        sf_case_control_bam_list, ncores
    );
    let s_bamsnap_hc = format!(
        concat!(
            "${{XTEA_MAIN}} --igv --bamsnap --single_sample -p ${{PREFIX}}\"tmp/igv\" -b {0} ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns_high_confident_post_filtering_somatic.txt\" ",
            " --ref ${{REF}} -e 1000 -n {1} ",
            "-o ${{PREFIX}}\"tmp/igv/bamsnap_screenshot_somatic.txt\"\n"
        ),
        sf_case_control_bam_list, ncores
    );
    let mut s_cmd = [s_case_control, s_case_control_hc, s_igv, s_igv_hc, s_bamsnap, s_bamsnap_hc].concat();
    let iflag = need_flag(a.iflag)?;
    if (iflag & 1024) == 1024 {
        s_cmd += &annotation_vcf_somatic(ncores, a.i_rep_type);
    }
    if iflag == 17408 {
        s_cmd = annotation_vcf_somatic(ncores, a.i_rep_type);
    }
    Ok(s_cmd)
}

/// `gnrt_calling_command_denovo_from_trio`; the hg38/CHM13 scripts skip the `--denovo` step
/// when flag 32768 is set (`skip_on_32768`), the hg19 script never does.
pub fn denovo_from_trio(a: &CallArgs, sf_ctr_list: &str, sf_case_control_bam_list: &str, skip_on_32768: bool) -> Result<String> {
    let (ncores, iflt_clip, iflt_disc, iflk_len) = (&a.ncores, a.iflt_clip, a.iflt_disc, a.iflk_len);
    let mut s_case_control = String::new();
    if !skip_on_32768 || (need_flag(a.iflag)? & 32768) != 32768 {
        s_case_control = format!(
            concat!(
                "${{XTEA_MAIN}} --denovo -p ${{TMP_CNS}} -b {0} ",
                "-n {1} -r ${{L1_CNS}} -i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\" ",
                " --input2 ${{PREFIX}}\"null\" --ref ${{REF}} --cr {2} --nd {3} ",
                "--fflank ${{SF_FLANK}} --flklen {4} ",
                "-o ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering_denovo.txt\"\n"
            ),
            sf_ctr_list, ncores, iflt_clip, iflt_disc, iflk_len
        );
    }
    // BamSnap
    let s_bamsnap = format!(
        concat!(
            "${{XTEA_MAIN}} --igv --bamsnap --single_sample -p ${{PREFIX}}\"tmp/igv\" -b {0} ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering_denovo.txt\" ",
            " --ref ${{REF}} -e 1000 -n {1} ",
            "-o ${{PREFIX}}\"tmp/igv/bamsnap_screenshot.txt\"\n"
        ),
        sf_case_control_bam_list, ncores
    );
    Ok(s_case_control + &s_bamsnap)
}

// ------------------------------------------------------------------------------------------
// Driver of the local generators (hg38, hg19, CHM13)
// ------------------------------------------------------------------------------------------

/// The parsed options of a local generator (`--hard` is always false for hg19, which lacks it).
#[derive(Debug, Clone)]
pub struct LocalOptions {
    pub decompress: bool,
    pub mosaic: bool,
    pub case_control: bool,
    pub denovo: bool,
    pub user: bool,
    pub force: bool,
    pub hard: bool,
    pub tumor: bool,
    pub purity: f64,
    pub lsf: bool,
    pub slurm: bool,
    pub resume: bool,
    pub version: bool,
    pub id: Option<String>,
    pub lib: Option<String>,
    pub bam: Option<String>,
    pub x10: Option<String>,
    pub cores: i64,
    pub memory: i64,
    pub partition: Option<String>,
    pub time: String,
    pub wfolder: Option<String>,
    pub r#ref: Option<String>,
    pub gene: Option<String>,
    pub xtea: Option<String>,
    pub flag: Option<i64>,
    pub rep_type: i64,
    pub flklen: i64,
    pub nclip: i64,
    pub cliprep: i64,
    pub ndisc: i64,
    pub nfilterclip: i64,
    pub nfilterdisc: i64,
    pub teilen: i64,
    pub output: String,
    pub blacklist: String,
}

/// What differs between the three local generators.
pub trait LocalVariant {
    /// Reference name printed by `-V` ("hg38", "hg19", "chm13").
    fn genome(&self) -> &'static str;
    /// Repeat types in `-y` bit order.
    fn rep_types(&self) -> &'static [RepType];
    /// `split_bam_list` skips BAM-list lines with fewer than two fields (hg19 only; the others
    /// fail on an empty line).
    fn split_bam_list_skips_short_lines(&self) -> bool;
    /// In de novo mode the hg19 script passes the `.case_list` file as BAM list, the others the
    /// original list.
    fn denovo_uses_case_list(&self) -> bool;
    /// `gnrt_calling_command` and the per-repeat-type overrides of `gnrt_running_shell`.
    fn calling_command(&self, a: &CallArgs, rep: &RepType) -> Result<String>;
    /// `gnrt_calling_command_somatic_case_control`
    fn somatic_case_control(&self, a: &CallArgs, sf_ctr_list: &str, sf_case_control_bam_list: &str) -> Result<String>;
    /// `gnrt_calling_command_denovo_from_trio`
    fn denovo_from_trio(&self, a: &CallArgs, sf_ctr_list: &str, sf_case_control_bam_list: &str) -> Result<String>;
}

/// `prepare_case_control_bam`: split "id case_bam ctrl_bam..." lines into the case list and the
/// control list.
pub fn prepare_case_control_bam(sf_ori_bam: &str, sf_sprt_bam: &str, sf_control_bam: &str) -> Result<()> {
    let data = std::fs::read(sf_ori_bam).map_err(|e| os_error(&e, sf_ori_bam))?;
    write_py_file(sf_sprt_bam, "")?;
    write_py_file(sf_control_bam, "")?;
    let mut s_sprt = String::new();
    let mut s_ctrl = String::new();
    for line in split_universal_lines(&decode_utf8(data)?) {
        let fields = py_split(&line);
        if fields.len() < 3 {
            py_print(&format!("Not in right format. Less fields: {}", line));
            continue;
        }
        s_sprt += &format!("{}\t{}\n", fields[0], fields[1]);
        s_ctrl += fields[0];
        for sf_tmp in &fields[2..] {
            s_ctrl += "\t";
            s_ctrl += sf_tmp;
        }
        s_ctrl += "\n";
    }
    write_py_file(sf_sprt_bam, &s_sprt)?;
    write_py_file(sf_control_bam, &s_ctrl)
}

/// The `__main__` block of the local generators.
pub fn run_local(v: &dyn LocalVariant, o: &LocalOptions) -> Result<()> {
    if o.version {
        py_print(&format!("xTea {} for short and linked reads on {}\n", S_VERSION, v.genome()));
        return Ok(());
    }
    let xtea_main = xtea_main_path(XTEA_MAIN_EXE);
    let s_wfolder = py_add_slash(o.wfolder.as_deref())?;
    if !exists(&s_wfolder) {
        popen(&format!("mkdir {}", s_wfolder))?;
    }
    // os.path.isfile(None) raises
    let isfile_or_null = |o: &Option<String>| -> Result<String> {
        let p = o.as_deref().ok_or_else(|| type_error("stat: path should be string, bytes, os.PathLike or integer, not NoneType"))?;
        Ok(if is_file(p) { p.to_string() } else { "null".to_string() })
    };
    let sf_bams = isfile_or_null(&o.bam)?;
    let sf_bams_10x = isfile_or_null(&o.x10)?;

    let mut sf_folder_rep = o.lib.clone();
    let mut sf_ref = o.r#ref.clone();
    if o.decompress {
        decompress(py_str_or_none(o.lib.as_deref()), &s_wfolder)?;
        decompress(py_str_or_none(o.r#ref.as_deref()), &s_wfolder)?;
        sf_folder_rep = Some(format!("{}rep_lib_annotation/", s_wfolder));
        sf_ref = Some(format!("{}genome.fa", s_wfolder));
    }
    let l_rep_type = select_rep_types(v.rep_types(), o.rep_type);
    let run = |sf_bams_arg: &str, sf_case_control_bam_list: &str, b_denovo: bool| {
        let ctx = RunCtx {
            v,
            o,
            sf_bams_10x_global: &sf_bams_10x,
            l_rep_type: &l_rep_type,
            sf_folder_rep: sf_folder_rep.as_deref(),
            sf_ref: sf_ref.as_deref(),
            xtea_main: &xtea_main,
        };
        gnrt_running_shell(&ctx, sf_bams_arg, &sf_bams_10x, &s_wfolder, sf_case_control_bam_list, b_denovo)
    };
    if !o.case_control {
        // non-case-ctrl mode (mainly for germline)
        run(&sf_bams, "null", false)
    } else {
        // case-control / de novo mode: prepare the case and the control lists first
        let sf_sprt_bams = format!("{}{}", sf_bams, CASE_LIST_SUFFIX);
        let sf_control_bams = format!("{}{}", sf_bams, CTRL_LIST_SUFFIX);
        prepare_case_control_bam(&sf_bams, &sf_sprt_bams, &sf_control_bams)?;
        if !o.denovo {
            run(&sf_sprt_bams, &sf_bams, false)
        } else {
            let sf_list = if v.denovo_uses_case_list() { &sf_sprt_bams } else { &sf_bams };
            run(sf_list, &sf_bams, true)
        }
    }
}

/// Values `gnrt_running_shell` reads from the `__main__` scope.
struct RunCtx<'a> {
    v: &'a dyn LocalVariant,
    o: &'a LocalOptions,
    /// the module-global `sf_bams_10X` (read by `split_bam_list`)
    sf_bams_10x_global: &'a str,
    l_rep_type: &'a [&'static RepType],
    sf_folder_rep: Option<&'a str>,
    sf_ref: Option<&'a str>,
    xtea_main: &'a str,
}

/// `gnrt_running_shell` of the local generators.
fn gnrt_running_shell(
    ctx: &RunCtx,
    sf_bams: &str,
    sf_10x_bams: &str,
    s_wfolder: &str,
    sf_case_control_bam_list: &str,
    b_denovo: bool,
) -> Result<()> {
    let (v, o) = (ctx.v, ctx.o);
    let s_wfolder = py_add_slash(Some(s_wfolder))?;
    if !exists(&s_wfolder) {
        popen(&format!("mkdir {}", s_wfolder))?;
    }
    let sf_ids = o.id.as_deref().ok_or_else(|| type_error("expected str, bytes or os.PathLike object, not NoneType"))?;
    let mut m_id: IndexSet<String> = IndexSet::new();
    for line in read_py_lines(sf_ids)? {
        let sid = py_rstrip(&line).to_string();
        m_id.insert(sid.clone());
        let sf_folder = format!("{}{}", s_wfolder, sid); // first create folder
        if exists(&sf_folder) {
            continue;
        }
        popen(&format!("mkdir {}", sf_folder))?;
        let sf_pub_clip = format!("{}/{}/", sf_folder, PUB_CLIP);
        popen(&format!("mkdir {}", sf_pub_clip))?;
    }

    // the sample, bam, x10 bam files
    split_bam_list(ctx, &m_id, sf_bams, sf_10x_bams, &s_wfolder, sf_case_control_bam_list)?;

    let mut l_sh: Vec<String> = Vec::new();
    for sid_tmp in &m_id {
        let sf_sample_folder = format!("{}{}/", s_wfolder, sid_tmp);
        let sf_pub_clip = format!("{}{}/", sf_sample_folder, PUB_CLIP);
        gnrt_lib_config(
            ctx.l_rep_type,
            ctx.sf_folder_rep,
            ctx.sf_ref,
            o.gene.as_deref(),
            &o.blacklist,
            o.xtea.as_deref(),
            &sf_sample_folder,
            ctx.xtea_main,
        )?;
        for rep in ctx.l_rep_type {
            let sf_config = format!("{}{}.config", sf_sample_folder, rep.name);
            let s_wfolder_rep = format!("{}{}", sf_sample_folder, rep.name);
            if !exists(&s_wfolder_rep) {
                run_cmd(&format!("mkdir {}", s_wfolder_rep))?;
            }
            let sf_rep_sample_id = format!("{}/sample_id.txt", s_wfolder_rep);
            let mut sf_rep_bam = format!("{}/bam_list1.txt", s_wfolder_rep);
            if !is_file(&sf_rep_bam) {
                sf_rep_bam = "null".to_string();
            }
            let mut sf_rep_x10_bam = format!("{}/x10_bam_list1.txt", s_wfolder_rep);
            if !exists(&sf_rep_x10_bam) {
                sf_rep_x10_bam = "null".to_string();
            }
            let mut s_head = String::new();
            if o.lsf {
                s_head = gnrt_script_head_lsf(o.partition.as_deref(), o.cores, &o.time, o.memory, sid_tmp)?;
            } else if o.slurm {
                s_head = gnrt_script_head(o.partition.as_deref(), o.cores, &o.time, o.memory, sid_tmp);
            }
            let l_libs = load_par_config(&sf_config)?;
            let s_libs = gnrt_parameters(&l_libs);
            let a = CallArgs {
                iclip_c: o.nclip,
                iclip_rp: o.cliprep,
                idisc_c: o.ndisc,
                iflt_clip: o.nfilterclip,
                iflt_disc: o.nfilterdisc,
                ncores: o.cores.to_string(),
                iflk_len: o.flklen,
                min_tei_len: o.teilen,
                iflag: o.flag,
                b_mosaic: o.mosaic,
                b_user_par: o.user,
                b_force: o.force,
                b_tumor: o.tumor,
                b_hard_cut: o.hard,
                b_resume: o.resume,
                f_purity: o.purity,
                i_rep_type: rep.bit,
                s_cfolder: sf_pub_clip.clone(),
            };
            let mut s_calling_cmd = v.calling_command(&a, rep)?;
            if sf_case_control_bam_list != "null" {
                let sf_ctr_list = format!("{}/bam_list2.txt", s_wfolder_rep);
                let sf_case_ctrl_list = format!("{}/bam_list3.txt", s_wfolder_rep);
                s_calling_cmd += &v.somatic_case_control(&a, &sf_ctr_list, &sf_case_ctrl_list)?;
                if b_denovo {
                    s_calling_cmd = v.denovo_from_trio(&a, &sf_ctr_list, &sf_case_ctrl_list)?;
                }
            }
            let l_tmp_sh = gnrt_pipelines(
                &s_head,
                &s_libs,
                &s_calling_cmd,
                &sf_rep_sample_id,
                &sf_rep_bam,
                &sf_rep_x10_bam,
                &sf_sample_folder,
                rep.name,
            )?;
            l_sh.extend(l_tmp_sh);
        }
    }
    let mut s_submit = String::from("#!/bin/bash\n\n");
    for s_sh in &l_sh {
        if !o.lsf {
            s_submit += &format!("sbatch < {}\n", s_sh);
        } else {
            s_submit += &format!("bsub < {}\n", s_sh);
        }
    }
    write_py_file(&o.output, &s_submit)
}

/// `split_bam_list`: write `sample_id.txt`, `bam_list1.txt` (case BAMs), `bam_list2.txt`
/// (controls), `bam_list3.txt` (case + controls) and `x10_bam_list1.txt` in every
/// sample/repeat folder.
fn split_bam_list(
    ctx: &RunCtx,
    m_ids: &IndexSet<String>,
    sf_bams: &str,
    sf_x10_bams: &str,
    s_wfolder: &str,
    sf_case_control_bam_list: &str,
) -> Result<()> {
    // load in sf_bams
    let mut m_bams: IndexMap<String, Vec<String>> = IndexMap::new();
    if sf_bams != "null" {
        for line in read_py_lines(sf_bams)? {
            let fields = py_split(&line);
            if ctx.v.split_bam_list_skips_short_lines() && fields.len() < 2 {
                continue;
            }
            let sid = field(&fields, 0)?;
            if !m_ids.contains(sid) {
                continue;
            }
            let s_bam = field(&fields, 1)?;
            m_bams.entry(sid.to_string()).or_default().push(s_bam.to_string());
        }
    }
    // load in sf_x10_bams; the original tests the module-global `sf_bams_10X` here, not its
    // parameter (both always hold the same value)
    let mut m_bams_10x: IndexMap<String, (String, String)> = IndexMap::new();
    if ctx.sf_bams_10x_global != "null" {
        for line in read_py_lines(sf_x10_bams)? {
            let fields = py_split(&line);
            let sid = field(&fields, 0)?;
            if !m_ids.contains(sid) {
                continue;
            }
            let s_bam = field(&fields, 1)?;
            let s_barcode_bam = field(&fields, 2)?;
            m_bams_10x.insert(sid.to_string(), (s_bam.to_string(), s_barcode_bam.to_string()));
        }
    }
    // load in the control bams
    let mut m_control: IndexMap<String, Vec<String>> = IndexMap::new();
    let mut m_case_control: IndexMap<String, Vec<String>> = IndexMap::new();
    if sf_case_control_bam_list != "null" {
        for line in read_py_lines(sf_case_control_bam_list)? {
            let fields = py_split(&line);
            if fields.len() < 3 {
                continue;
            }
            let s_id = fields[0].to_string();
            m_control.insert(s_id.clone(), fields[2..].iter().map(|s| s.to_string()).collect());
            m_case_control.insert(s_id, fields[1..].iter().map(|s| s.to_string()).collect());
        }
    }

    for sid_tmp in m_ids {
        let sf_sample_folder = format!("{}{}/", s_wfolder, sid_tmp);
        if !exists(&sf_sample_folder) {
            run_cmd(&format!("mkdir {}", sf_sample_folder))?;
        }
        for rep in ctx.l_rep_type {
            let s_wfolder_rep = format!("{}{}", sf_sample_folder, rep.name);
            if !exists(&s_wfolder_rep) {
                run_cmd(&format!("mkdir {}", s_wfolder_rep))?;
            }
            write_py_file(&format!("{}/sample_id.txt", s_wfolder_rep), sid_tmp)?;
            if let Some(bams) = m_bams.get(sid_tmp) {
                let s: String = bams.iter().map(|b| format!("{}\t{}\n", sid_tmp, b)).collect();
                write_py_file(&format!("{}/bam_list1.txt", s_wfolder_rep), &s)?;
            }
            if let Some(ctrl) = m_control.get(sid_tmp) {
                // the control file list
                let s: String = ctrl.iter().map(|b| format!("{}\n", b)).collect();
                write_py_file(&format!("{}/bam_list2.txt", s_wfolder_rep), &s)?;
            }
            if m_control.contains_key(sid_tmp) {
                // the case + control file list
                let mut s = sid_tmp.clone();
                for b in &m_case_control[sid_tmp] {
                    s += "\t";
                    s += b;
                }
                s += "\n";
                write_py_file(&format!("{}/bam_list3.txt", s_wfolder_rep), &s)?;
            }
            if let Some((s_bam, s_barcode_bam)) = m_bams_10x.get(sid_tmp) {
                write_py_file(&format!("{}/x10_bam_list1.txt", s_wfolder_rep), &format!("{}\t{}\t{}\n", sid_tmp, s_bam, s_barcode_bam))?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn true_division_like_python() {
        assert_eq!(float_repr(py_int_truediv(20, 1).unwrap()), "20.0");
        assert_eq!(float_repr(py_int_truediv(25, 8).unwrap()), "3.125");
        assert_eq!(float_repr(py_int_truediv(20, 3).unwrap()), "6.666666666666667");
        assert_eq!(float_repr(py_int_truediv(-7, 2).unwrap()), "-3.5");
        assert_eq!(float_repr(py_int_truediv(0, -5).unwrap()), "-0.0");
        assert_eq!(float_repr(py_int_truediv(1, 3).unwrap()), "0.3333333333333333");
        assert_eq!(float_repr(py_int_truediv(i64::MAX, 3).unwrap()), "3.0744573456182584e+18");
        assert_eq!(float_repr(py_int_truediv(9007199254740993, 1).unwrap()), "9007199254740992.0");
        assert_eq!(float_repr(py_int_truediv(9007199254740995, 1).unwrap()), "9007199254740996.0");
        assert_eq!(float_repr(py_int_truediv(i64::MIN, 7).unwrap()), "-1.3176245766935393e+18");
        assert!(py_int_truediv(1, 0).is_err());
    }

    #[test]
    fn universal_newlines() {
        assert_eq!(split_universal_lines("a\r\nb\rc\nd"), vec!["a\n", "b\n", "c\n", "d"]);
        assert_eq!(split_universal_lines("\r\r\n\n"), vec!["\n", "\n", "\n"]);
        assert!(split_universal_lines("").is_empty());
    }

    #[test]
    fn python_whitespace() {
        assert_eq!(py_split(" a\u{a0}b\x1cc\u{3000}d \n"), vec!["a", "b", "c", "d"]);
        assert_eq!(py_rstrip("ab \x1f\u{2028}\t\n"), "ab");
        assert_eq!(py_split("a\u{feff}b"), vec!["a\u{feff}b"]);
    }

    #[test]
    fn list_repr() {
        assert_eq!(py_list_repr(&[]), "[]");
        assert_eq!(py_list_repr(&["/a b/x.sh".to_string(), "it's".to_string()]), "['/a b/x.sh', \"it's\"]");
    }
}
