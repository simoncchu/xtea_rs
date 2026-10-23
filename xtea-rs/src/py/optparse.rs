//! Value conversions of Python's `optparse`, used as clap value parsers so that the
//! command lines of the Rust binaries accept exactly what the original scripts accepted.

/// Python `int(s, radix)` for a string (whitespace, sign, base prefix of the same radix and
/// single underscores between digits are accepted).
fn py_int_radix(s: &str, radix: u32) -> Option<i64> {
    let t = s.trim();
    let (neg, body) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    let lower = body.to_ascii_lowercase();
    let body = match radix {
        16 if lower.starts_with("0x") => &body[2..],
        8 if lower.starts_with("0o") => &body[2..],
        2 if lower.starts_with("0b") => &body[2..],
        _ => body,
    };
    // underscores are allowed between digits (and right after a base prefix)
    let had_prefix = body.len() != lower.len();
    let mut digits = String::with_capacity(body.len());
    let mut prev_us = !had_prefix; // a leading underscore is only valid after a prefix
    for (i, c) in body.chars().enumerate() {
        if c == '_' {
            if prev_us || i + 1 == body.len() {
                return None;
            }
            prev_us = true;
            continue;
        }
        prev_us = false;
        digits.push(c);
    }
    if digits.is_empty() {
        return None;
    }
    let v = i64::from_str_radix(&digits, radix).ok()?;
    Some(if neg { -v } else { v })
}

/// optparse `type="int"` (`optparse._parse_num`): "0x.." is hexadecimal, "0b.." binary and
/// any other value starting with "0" octal.
pub fn parse_int(s: &str) -> Result<i64, String> {
    let lower2: String = s.chars().take(2).collect::<String>().to_ascii_lowercase();
    let first = s.chars().next();
    let v = if lower2 == "0x" {
        py_int_radix(s, 16)
    } else if lower2 == "0b" {
        let rest = &s[2..];
        py_int_radix(if rest.is_empty() { "0" } else { rest }, 2)
    } else if first == Some('0') {
        py_int_radix(s, 8)
    } else {
        py_int_radix(s, 10)
    };
    v.ok_or_else(|| format!("invalid integer value: {}", crate::py::fmt::py_str_repr(s)))
}

/// optparse `type="float"` (Python `float(s)`).
pub fn parse_float(s: &str) -> Result<f64, String> {
    crate::py::text::py_float(s).map_err(|_| format!("invalid floating-point value: {}", crate::py::fmt::py_str_repr(s)))
}

// ------------------------------------------------------------------------------------------
// optparse command-line semantics on top of clap
// ------------------------------------------------------------------------------------------

use std::ffi::OsString;
use std::io::Write;
use std::path::Path;


/// The long option (name, takes a value) and the short options of a clap command, in definition
/// order.  The generators define their own `-h/--help` flag (clap's is disabled) so that clap's
/// prefix matching sees it, as optparse's does (`--h` is ambiguous with `--hard`).
fn option_specs(cmd: &clap::Command) -> (Vec<(String, bool)>, Vec<(char, usize)>, Vec<(char, bool)>) {
    let mut cmd = cmd.clone();
    cmd.build();
    let mut specs = Vec::new();
    let mut shorts = Vec::new();
    let mut short_only = Vec::new();
    for a in cmd.get_arguments() {
        if a.is_positional() {
            continue;
        }
        if let Some(long) = a.get_long() {
            if let Some(s) = a.get_short() {
                shorts.push((s, specs.len()));
            }
            specs.push((long.to_string(), a.get_action().takes_values()));
        } else if let Some(s) = a.get_short() {
            short_only.push((s, a.get_action().takes_values()));
        }
    }
    (specs, shorts, short_only)
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
    /// the message of the first usage error optparse would report (unknown or ambiguous
    /// option, missing value, value given to a flag)
    pub error: Option<String>,
    /// how each value-taking option was first spelled on the command line ("-n" or "--cores"),
    /// keyed by its long name, for optparse's "option -n: invalid ... value" messages
    pub spelled: std::collections::HashMap<String, String>,
}

/// optparse's `_match_abbrev` failure message for a long option.
fn long_option_error(specs: &[(String, bool)], name: &str) -> String {
    let mut poss: Vec<&str> = specs.iter().filter(|s| s.0.starts_with(name)).map(|s| s.0.as_str()).collect();
    if poss.is_empty() {
        format!("no such option: --{}", name)
    } else {
        poss.sort();
        let l: Vec<String> = poss.iter().map(|p| format!("--{}", p)).collect();
        format!("ambiguous option: --{} ({}?)", name, l.join(", "))
    }
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
    let (specs, shorts, short_only) = option_specs(cmd);
    let mut it = argv.into_iter();
    let mut out: Vec<OsString> = it.next().into_iter().collect();
    let rest: Vec<OsString> = it.collect();
    let args: Vec<String> = match rest.iter().map(|a| a.to_str().map(str::to_string)).collect::<Option<Vec<_>>>() {
        Some(v) => v,
        None => {
            out.extend(rest);
            return OptparseArgv { argv: out, help: false, error: None, spelled: Default::default() };
        }
    };
    let mut spelled: std::collections::HashMap<String, String> = Default::default();
    let done = |mut out: Vec<OsString>, from: &[String], spelled: std::collections::HashMap<String, String>| {
        out.extend(from.iter().map(OsString::from));
        OptparseArgv { argv: out, help: false, error: None, spelled }
    };
    let fail = |mut out: Vec<OsString>, from: &[String], msg: String| {
        out.extend(from.iter().map(OsString::from));
        OptparseArgv { argv: out, help: false, error: Some(msg), spelled: Default::default() }
    };
    let help = |out: Vec<OsString>| OptparseArgv { argv: out, help: true, error: None, spelled: Default::default() };
    let mut i = 0;
    while i < args.len() {
        let arg = &args[i];
        if arg == "--" {
            return done(out, &args[i..], spelled);
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
                        None => return fail(out, &args[i..], format!("--{} option requires 1 argument", long)),
                    };
                    spelled.entry(long.clone()).or_insert_with(|| format!("--{}", long));
                    out.push(format!("--{}={}", long, value).into());
                }
                Some((long, false)) if explicit.is_none() => {
                    if long == "help" {
                        return help(out);
                    }
                    out.push(format!("--{}", long).into());
                }
                Some((long, false)) => {
                    return fail(out, &args[i..], format!("--{} option does not take a value", long));
                }
                None => {
                    let msg = long_option_error(&specs, name);
                    return fail(out, &args[i..], msg);
                }
            }
        } else if arg.len() > 1 && arg.starts_with('-') {
            let chars: Vec<char> = arg[1..].chars().collect();
            let mut j = 0;
            while j < chars.len() {
                let (long, takes_value) = match shorts.iter().find(|(c, _)| *c == chars[j]) {
                    Some(&(_, idx)) => &specs[idx],
                    None => match short_only.iter().find(|(c, _)| *c == chars[j]) {
                        // an option without a long name: passed on as "-c" (and its value as a
                        // separate argument, which clap takes verbatim)
                        Some(&(c, so_takes_value)) => {
                            if so_takes_value {
                                let value = if j + 1 < chars.len() {
                                    chars[j + 1..].iter().collect::<String>()
                                } else if i + 1 < args.len() {
                                    i += 1;
                                    args[i].clone()
                                } else {
                                    out.push(format!("-{}", c).into());
                                    return fail(out, &[], format!("-{} option requires 1 argument", c));
                                };
                                spelled.entry(c.to_string()).or_insert_with(|| format!("-{}", c));
                                out.push(format!("-{}", c).into());
                                out.push(value.into());
                                break;
                            }
                            out.push(format!("-{}", c).into());
                            j += 1;
                            continue;
                        }
                        None => {
                            out.push(format!("-{}", chars[j..].iter().collect::<String>()).into());
                            return fail(out, &args[i + 1..], format!("no such option: -{}", chars[j]));
                        }
                    },
                };
                if *takes_value {
                    let value = if j + 1 < chars.len() {
                        chars[j + 1..].iter().collect::<String>()
                    } else if i + 1 < args.len() {
                        i += 1;
                        args[i].clone()
                    } else {
                        out.push(format!("-{}", chars[j]).into());
                        return fail(out, &[], format!("-{} option requires 1 argument", chars[j]));
                    };
                    spelled.entry(long.clone()).or_insert_with(|| format!("-{}", chars[j]));
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
    OptparseArgv { argv: out, help: false, error: None, spelled }
}

/// optparse's `parser.error(msg)`: usage line and message on stderr, exit status 2.
fn optparse_exit(prog: &str, msg: &str) -> ! {
    let _ = std::io::stdout().flush();
    eprint!("Usage: {} [options]\n\n{}: error: {}\n", prog, prog, msg);
    std::process::exit(2);
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
    if let Some(msg) = &parsed.error {
        optparse_exit(&prog, msg);
    }
    let spelled = parsed.spelled.clone();
    let matches = match cmd.clone().try_get_matches_from(parsed.argv) {
        Ok(m) => m,
        Err(e) => {
            use clap::error::{ContextKind, ContextValue, ErrorKind};
            if matches!(e.kind(), ErrorKind::ValueValidation | ErrorKind::InvalidValue) {
                // "option --cores: invalid integer value: '08'"
                let opt = match e.get(ContextKind::InvalidArg) {
                    Some(ContextValue::String(a)) => a.split_whitespace().next().unwrap_or("").to_string(),
                    _ => String::new(),
                };
                let opt = spelled.get(opt.trim_start_matches('-')).cloned().unwrap_or(opt);
                let detail = std::error::Error::source(&e).map(|s| s.to_string());
                if let Some(d) = detail {
                    optparse_exit(&prog, &format!("option {}: {}", opt, d));
                }
            }
            e.exit()
        }
    };
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
    let t = v.trim_matches(crate::py::text::is_py_space);
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
    let blank = |c: &str| c.chars().all(crate::py::text::is_py_space);
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


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ints_like_optparse() {
        assert_eq!(parse_int("010"), Ok(8));
        assert_eq!(parse_int("0x1f"), Ok(31));
        assert_eq!(parse_int("0b101"), Ok(5));
        assert_eq!(parse_int(" 12"), Ok(12));
        assert_eq!(parse_int("+5"), Ok(5));
        assert_eq!(parse_int("-7"), Ok(-7));
        assert_eq!(parse_int("-010"), Ok(-10));
        assert_eq!(parse_int("1_000"), Ok(1000));
        assert_eq!(parse_int("0"), Ok(0));
        assert_eq!(parse_int("0o17"), Ok(15));
        assert!(parse_int("08").is_err());
        assert!(parse_int("1__0").is_err());
        assert!(parse_int("abc").is_err());
        assert!(parse_int("").is_err());
        assert_eq!(parse_int("5907"), Ok(5907));
    }

    #[test]
    fn floats_like_python() {
        assert_eq!(parse_float("0.45"), Ok(0.45));
        assert_eq!(parse_float(" 1e-3 "), Ok(0.001));
        assert_eq!(parse_float("1_0.5"), Ok(10.5));
        assert!(parse_float("x").is_err());
    }
}
