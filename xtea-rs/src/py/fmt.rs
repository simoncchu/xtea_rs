//! Python-compatible number formatting and rounding helpers.

/// Split the shortest round-trip representation of a finite, non-zero |x| into
/// (digits, decpt) such that |x| = 0.d1d2d3... * 10^decpt (David Gay's dtoa convention).
fn shortest_digits(x: f64) -> (String, i32) {
    let s = format!("{:e}", x.abs());
    let (mant, exp) = s.split_once('e').expect("exponent format");
    let exp: i32 = exp.parse().expect("exponent value");
    let digits: String = mant.chars().filter(|c| *c != '.').collect();
    (digits, exp + 1)
}

/// `repr(float)` / `str(float)` as produced by CPython 3.
pub fn float_repr(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".to_string() } else { "-inf".to_string() };
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0.0".to_string() } else { "0.0".to_string() };
    }
    let sign = if x < 0.0 { "-" } else { "" };
    let (digits, decpt) = shortest_digits(x);
    let nd = digits.len() as i32;
    let body = if decpt <= -4 || decpt > 16 {
        // exponential notation
        let mut m = String::new();
        m.push_str(&digits[0..1]);
        if nd > 1 {
            m.push('.');
            m.push_str(&digits[1..]);
        }
        let e = decpt - 1;
        let esign = if e < 0 { '-' } else { '+' };
        format!("{}e{}{:02}", m, esign, e.abs())
    } else if decpt <= 0 {
        format!("0.{}{}", "0".repeat((-decpt) as usize), digits)
    } else if decpt >= nd {
        format!("{}{}.0", digits, "0".repeat((decpt - nd) as usize))
    } else {
        let (a, b) = digits.split_at(decpt as usize);
        format!("{}.{}", a, b)
    };
    format!("{}{}", sign, body)
}

/// Python's built-in `round(x, ndigits)` for floats (correctly rounded, ties-to-even on the
/// exact binary value).
pub fn py_round(x: f64, ndigits: i32) -> f64 {
    if !x.is_finite() {
        return x;
    }
    if ndigits >= 0 {
        let s = format!("{:.*}", ndigits as usize, x);
        s.parse::<f64>().unwrap_or(x)
    } else {
        let p = 10f64.powi(-ndigits);
        let y = x / p;
        // round half to even of the exact quotient approximation
        let r = round_half_even(y);
        r * p
    }
}

/// Python's built-in `round(x)` (no ndigits): banker's rounding to an integer.
pub fn py_round_to_int(x: f64) -> i64 {
    round_half_even(x) as i64
}

/// Round half to even (like C `rint` in the default rounding mode).
pub fn round_half_even(x: f64) -> f64 {
    let r = x.round();
    if (x - x.trunc()).abs() == 0.5 {
        // tie: choose even
        let t = x.trunc();
        if (t as i64) % 2 == 0 {
            t
        } else {
            r
        }
    } else {
        r
    }
}

/// numpy's `round(np.float64(x), decimals)` / `np.around`: multiply, rint, divide.
pub fn np_round(x: f64, decimals: i32) -> f64 {
    if decimals >= 0 {
        let f = 10f64.powi(decimals);
        round_half_even(x * f) / f
    } else {
        let f = 10f64.powi(-decimals);
        round_half_even(x / f) * f
    }
}

/// Python `int(x)` for a float: truncation toward zero.
#[inline]
pub fn py_int(x: f64) -> i64 {
    x.trunc() as i64
}

/// Python integer floor division `a // b`.
#[inline]
pub fn floordiv(a: i64, b: i64) -> i64 {
    let q = a / b;
    if (a % b != 0) && ((a < 0) != (b < 0)) {
        q - 1
    } else {
        q
    }
}

/// Python integer modulo `a % b` (result has the sign of `b`).
#[inline]
pub fn pymod(a: i64, b: i64) -> i64 {
    let r = a % b;
    if r != 0 && ((r < 0) != (b < 0)) {
        r + b
    } else {
        r
    }
}

/// `repr(str)` as produced by CPython 3 (quote selection and escapes of `unicode_repr`).
pub fn py_str_repr(s: &str) -> String {
    let has_sq = s.contains('\'');
    let has_dq = s.contains('"');
    let quote = if has_sq && !has_dq { '"' } else { '\'' };
    let mut out = String::with_capacity(s.len() + 2);
    out.push(quote);
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            c if (c as u32) < 0x20 || c as u32 == 0x7f => out.push_str(&format!("\\x{:02x}", c as u32)),
            c if (c as u32) < 0x7f => out.push(c),
            c if c.is_control() => {
                let v = c as u32;
                if v <= 0xff {
                    out.push_str(&format!("\\x{:02x}", v));
                } else if v <= 0xffff {
                    out.push_str(&format!("\\u{:04x}", v));
                } else {
                    out.push_str(&format!("\\U{:08x}", v));
                }
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}

/// Python `float(s)` parser for the subset of inputs xTea produces (accepts surrounding spaces,
/// `inf`, `nan`).
pub fn parse_py_float(s: &str) -> Option<f64> {
    let t = s.trim();
    match t.to_ascii_lowercase().as_str() {
        "inf" | "+inf" | "infinity" | "+infinity" => return Some(f64::INFINITY),
        "-inf" | "-infinity" => return Some(f64::NEG_INFINITY),
        "nan" | "+nan" | "-nan" => return Some(f64::NAN),
        _ => {}
    }
    t.parse::<f64>().ok()
}

/// Python `int(s)` parser (accepts surrounding whitespace, optional sign and underscores are
/// not accepted here since xTea never writes them).
pub fn parse_py_int(s: &str) -> Option<i64> {
    let t = s.trim();
    let t = t.strip_prefix('+').unwrap_or(t);
    t.parse::<i64>().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn str_repr_matches_python() {
        assert_eq!(py_str_repr("abc"), "'abc'");
        assert_eq!(py_str_repr("a'b"), "\"a'b\"");
        assert_eq!(py_str_repr("a'b\""), "'a\\'b\"'");
        assert_eq!(py_str_repr("x\ty\n"), "'x\\ty\\n'");
        assert_eq!(py_str_repr("\x01\\"), "'\\x01\\\\'");
    }

    #[test]
    fn repr_matches_python() {
        let cases: &[(f64, &str)] = &[
            (0.125, "0.125"),
            (1e16, "1e+16"),
            (1e15, "1000000000000000.0"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (123.456, "123.456"),
            (5e-324, "5e-324"),
            (1.7976931348623157e308, "1.7976931348623157e+308"),
            (0.1 + 0.2, "0.30000000000000004"),
            (1.0, "1.0"),
            (-2.5, "-2.5"),
            (100.0, "100.0"),
            (12345678901234567.0, "1.2345678901234568e+16"),
            (0.0, "0.0"),
        ];
        for (x, s) in cases {
            assert_eq!(float_repr(*x), *s, "x={}", x);
        }
    }

    #[test]
    fn rounding() {
        assert_eq!(py_round(0.125, 2), 0.12);
        assert_eq!(py_round(2.675, 2), 2.67);
        assert_eq!(py_round(1.005, 2), 1.0);
        assert_eq!(py_round_to_int(2.5), 2);
        assert_eq!(py_round_to_int(3.5), 4);
        assert_eq!(np_round(0.125, 2), 0.12);
        assert_eq!(floordiv(-7, 2), -4);
        assert_eq!(pymod(-7, 2), 1);
    }
}
