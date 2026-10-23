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
    let t = s.trim();
    let cleaned: String;
    let t = if t.contains('_') {
        // Python accepts single underscores between digits
        let b = t.as_bytes();
        for (i, &c) in b.iter().enumerate() {
            if c == b'_' && (i == 0 || i + 1 == b.len() || !b[i - 1].is_ascii_digit() || !b[i + 1].is_ascii_digit()) {
                return Err(format!("invalid floating-point value: {}", crate::py::fmt::py_str_repr(s)));
            }
        }
        cleaned = t.replace('_', "");
        cleaned.as_str()
    } else {
        t
    };
    crate::py::fmt::parse_py_float(t).ok_or_else(|| format!("invalid floating-point value: {}", crate::py::fmt::py_str_repr(s)))
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
