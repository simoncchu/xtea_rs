//! Bit-exact emulation of the numpy reductions xTea relies on.

const PW_BLOCKSIZE: usize = 128;

/// numpy's pairwise summation (`@TYPE@_pairwise_sum` in loops_utils.h).
pub fn pairwise_sum(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut res = 0.0f64;
        for &x in a {
            res += x;
        }
        // numpy starts the inner reduction loop from the identity (0.0); `-0.0 + 0.0` etc.
        // are irrelevant for xTea's data.
        res
    } else if n <= PW_BLOCKSIZE {
        let mut r = [0.0f64; 8];
        r.copy_from_slice(&a[0..8]);
        let mut i = 8;
        while i < n - (n % 8) {
            for j in 0..8 {
                r[j] += a[i + j];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        pairwise_sum(&a[..n2]) + pairwise_sum(&a[n2..])
    }
}

/// `np.add.reduce` over a contiguous float64 array (identity 0.0 + pairwise sum).
pub fn np_sum(a: &[f64]) -> f64 {
    0.0 + pairwise_sum(a)
}

/// `np.mean` of a float64 array.
pub fn np_mean(a: &[f64]) -> f64 {
    np_sum(a) / (a.len() as f64)
}

/// `np.std(a)` (ddof = 0) following numpy/core/_methods.py `_var`.
pub fn np_std(a: &[f64]) -> f64 {
    let n = a.len() as f64;
    let mean = np_sum(a) / n;
    let sq: Vec<f64> = a
        .iter()
        .map(|&x| {
            let d = x - mean;
            d * d
        })
        .collect();
    (np_sum(&sq) / n).sqrt()
}

/// `np.std` of an integer list (converted to float64 first, as numpy does).
pub fn np_std_i64(a: &[i64]) -> f64 {
    let v: Vec<f64> = a.iter().map(|&x| x as f64).collect();
    np_std(&v)
}

/// `np.mean` of an integer list.
pub fn np_mean_i64(a: &[i64]) -> f64 {
    let v: Vec<f64> = a.iter().map(|&x| x as f64).collect();
    np_mean(&v)
}

/// `np.median` of a list of floats (sorted copy; average of the two middle values).
pub fn np_median(a: &[f64]) -> f64 {
    let mut v = a.to_vec();
    v.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 1 {
        v[n / 2]
    } else {
        // numpy computes mean of the two middle elements via np.mean -> (a+b)/2 pairwise
        np_mean(&v[n / 2 - 1..n / 2 + 1])
    }
}
