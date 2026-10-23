//! Port of `x_polyA.py`.

use std::collections::BTreeSet;

#[inline]
fn is_a(c: u8) -> bool {
    c == b'A' || c == b'a'
}
#[inline]
fn is_t(c: u8) -> bool {
    c == b'T' || c == b't'
}

/// `is_poly_A_T`
pub fn is_poly_a_t(seq: &str) -> bool {
    let b = seq.as_bytes();
    let cnt_a = b.iter().filter(|&&c| is_a(c)).count() as f64;
    let cnt_t = b.iter().filter(|&&c| is_t(c)).count() as f64;
    let cutoff = b.len() as f64 * 0.75;
    cnt_a > cutoff || cnt_t > cutoff
}

/// `contain_poly_A_T`
pub fn contain_poly_a_t(seq: &str, n_min_cnt: i64) -> bool {
    let (mut max_a, mut max_t, mut cum_a, mut cum_t) = (0i64, 0i64, 0i64, 0i64);
    for &ch in seq.as_bytes() {
        if is_a(ch) {
            cum_a += 1;
        } else {
            if cum_a > max_a {
                max_a = cum_a;
            }
            cum_a = 0;
        }
        if is_t(ch) {
            cum_t += 1;
        } else {
            if cum_t > max_t {
                max_t = cum_t;
            }
            cum_t = 0;
        }
    }
    if cum_a > max_a {
        max_a = cum_a;
    }
    if cum_t > max_t {
        max_t = cum_t;
    }
    max_a >= n_min_cnt || max_t >= n_min_cnt
}

/// `is_consecutive_polyA_T`
pub fn is_consecutive_polya_t(seq: &str) -> bool {
    seq.contains("AAAAA") || seq.contains("TTTTT") || seq.contains("AATAA") || seq.contains("TTATT")
}

/// `is_consecutive_polyA_T2`
pub fn is_consecutive_polya_t2(seq: &str) -> bool {
    seq.contains("AAAAAA") || seq.contains("TTTTTT") || seq.contains("AAATAA") || seq.contains("TTTATT")
}

/// `is_consecutive_polyA`
pub fn is_consecutive_polya(seq: &str) -> bool {
    seq.contains("AAAAA") || seq.contains("AATAA")
}

/// `is_consecutive_polyA_T_with_ori` -> (found, is_polyA)
pub fn is_consecutive_polya_t_with_ori(seq: &str) -> (bool, bool) {
    if seq.contains("AAAAA") || seq.contains("AATAA") {
        (true, true)
    } else if seq.contains("TTTTT") || seq.contains("TTATT") {
        (true, false)
    } else {
        (false, false)
    }
}

/// `is_consecutive_polyA_T_with_oritation`
pub fn is_consecutive_polya_t_with_oritation(seq: &str, b_left: bool) -> bool {
    if b_left {
        seq.contains("AAAAA") || seq.contains("AATAA")
    } else {
        seq.contains("TTTTT") || seq.contains("TTATT")
    }
}

/// `contain_enough_A_T`
pub fn contain_enough_a_t(seq: &str, n_min_a_t: i64) -> bool {
    let b = seq.as_bytes();
    let n_a = b.iter().filter(|&&c| is_a(c)).count() as i64;
    let n_t = b.iter().filter(|&&c| is_t(c)).count() as i64;
    n_a >= n_min_a_t || n_t >= n_min_a_t
}

/// `contain_polyA_T`
pub fn contain_polya_t(s_seq: &str, b_rc: bool) -> bool {
    if !b_rc {
        s_seq.contains("AATAAA") || s_seq.contains("AAAAAA") || s_seq.contains("ATTAAA")
    } else {
        s_seq.contains("TTTATT") || s_seq.contains("TTTTTT") || s_seq.contains("TTTAAT")
    }
}

fn find_all_nonoverlapping(s: &str, pat: &str) -> Vec<usize> {
    // re.finditer semantics: non-overlapping matches scanning left to right
    let mut v = Vec::new();
    let mut i = 0;
    let sb = s.as_bytes();
    let pb = pat.as_bytes();
    while i + pb.len() <= sb.len() {
        if &sb[i..i + pb.len()] == pb {
            v.push(i);
            i += pb.len();
        } else {
            i += 1;
        }
    }
    v
}

/// `search_multi_polyA_locations`
pub fn search_multi_polya_locations(s_seq: &str, b_rc: bool) -> Vec<usize> {
    let pats = if b_rc { ["TTTATT", "TTTTTT", "TTTAAT"] } else { ["AATAAA", "AAAAAA", "ATTAAA"] };
    let mut v = Vec::new();
    for p in pats.iter() {
        v.extend(find_all_nonoverlapping(s_seq, p));
    }
    v
}

/// `is_dominant_polyA`
pub fn is_dominant_polya(seq: &str, f_ratio: f64) -> bool {
    let b = seq.as_bytes();
    let n_len = b.len();
    if n_len == 0 {
        return false;
    }
    let n_a = b.iter().filter(|&&c| is_a(c)).count() as f64;
    let n_t = b.iter().filter(|&&c| is_t(c)).count() as f64;
    n_a / n_len as f64 > f_ratio || n_t / n_len as f64 > f_ratio
}

/// `is_dominant_A`
pub fn is_dominant_a(seq: &str, f_ratio: f64) -> bool {
    let b = seq.as_bytes();
    if b.is_empty() {
        return false;
    }
    let n_a = b.iter().filter(|&&c| is_a(c)).count() as f64;
    n_a / b.len() as f64 > f_ratio
}

fn rotations(s: &str) -> Vec<String> {
    let n = s.len();
    let t = format!("{}{}", s, s);
    (0..n).map(|i| t[i..i + n].to_string()).collect()
}

fn build_set(base: &[&str]) -> BTreeSet<String> {
    let mut m = BTreeSet::new();
    for s in base {
        m.insert(s.to_string());
        // s_tmp[1:-2] strips "(" and ")n"
        let inner = &s[1..s.len() - 2];
        for r in rotations(inner) {
            m.insert(format!("({})n", r));
        }
    }
    m
}

/// `get_pre_defined_polyA_in_rmsk`
pub fn get_pre_defined_polya_in_rmsk() -> BTreeSet<String> {
    build_set(&[
        "(AAAAAC)n", "(AAAAC)n", "(AAAATAA)n", "(AAAAT)n", "(AAATA)n", "(AAAT)n", "(AACAAA)n", "(AAC)n", "(AATAAAA)n", "(AATAAA)n",
        "(AATAA)n", "(AATA)n", "(AAT)n", "(ACA)n", "(A)n",
    ])
}

/// `get_pre_defined_polyT_in_rmsk`
pub fn get_pre_defined_polyt_in_rmsk() -> BTreeSet<String> {
    build_set(&[
        "(GTTTTT)n", "(GTTTT)n", "(TTATTTT)n", "(ATTTT)n", "(TATTT)n", "(ATTT)n", "(TTTGTT)n", "(GTT)n", "(TTTTGTT)n", "(TTTATT)n",
        "(TTATT)n", "(TATT)n", "(ATT)n", "(TGT)n", "(T)n",
    ])
}
