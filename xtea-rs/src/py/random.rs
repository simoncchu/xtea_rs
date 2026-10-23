//! CPython-compatible `random` module subset (MT19937 + randint/randrange/random/sample).
//!
//! Python's `random` is seeded from OS entropy unless `random.seed()` is called.  The Rust
//! port seeds from entropy as well, unless a seed is supplied (`--seed` / `XTEA_SEED`), in
//! which case the produced sequence is identical to CPython after `random.seed(seed)`.

const N: usize = 624;
const M: usize = 397;
const MATRIX_A: u32 = 0x9908b0df;
const UPPER_MASK: u32 = 0x80000000;
const LOWER_MASK: u32 = 0x7fffffff;

#[derive(Clone)]
pub struct PyRandom {
    mt: [u32; N],
    mti: usize,
}

impl PyRandom {
    fn init_genrand(&mut self, s: u32) {
        self.mt[0] = s;
        for i in 1..N {
            let prev = self.mt[i - 1];
            self.mt[i] = 1812433253u32
                .wrapping_mul(prev ^ (prev >> 30))
                .wrapping_add(i as u32);
        }
        self.mti = N;
    }

    fn init_by_array(&mut self, key: &[u32]) {
        self.init_genrand(19650218);
        let mut i: usize = 1;
        let mut j: usize = 0;
        let key_length = key.len().max(1);
        let mut k = if N > key_length { N } else { key_length };
        while k > 0 {
            let prev = self.mt[i - 1];
            let kj = if key.is_empty() { 0 } else { key[j] };
            self.mt[i] = (self.mt[i] ^ ((prev ^ (prev >> 30)).wrapping_mul(1664525)))
                .wrapping_add(kj)
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= N {
                self.mt[0] = self.mt[N - 1];
                i = 1;
            }
            if j >= key_length {
                j = 0;
            }
            k -= 1;
        }
        k = N - 1;
        while k > 0 {
            let prev = self.mt[i - 1];
            self.mt[i] = (self.mt[i] ^ ((prev ^ (prev >> 30)).wrapping_mul(1566083941)))
                .wrapping_sub(i as u32);
            i += 1;
            if i >= N {
                self.mt[0] = self.mt[N - 1];
                i = 1;
            }
            k -= 1;
        }
        self.mt[0] = 0x80000000;
    }

    /// Equivalent of `random.seed(n)` for a non-negative integer `n`.
    pub fn from_seed(seed: u64) -> Self {
        let mut r = PyRandom { mt: [0; N], mti: N + 1 };
        let mut key: Vec<u32> = Vec::new();
        let mut v = seed;
        if v == 0 {
            key.push(0);
        }
        while v > 0 {
            key.push((v & 0xffff_ffff) as u32);
            v >>= 32;
        }
        r.init_by_array(&key);
        r
    }

    /// Seed from OS entropy (like an unseeded CPython `random`).
    pub fn from_entropy() -> Self {
        use std::io::Read;
        let mut buf = [0u8; 8];
        let ok = std::fs::File::open("/dev/urandom")
            .and_then(|mut f| f.read_exact(&mut buf))
            .is_ok();
        let seed = if ok {
            u64::from_le_bytes(buf)
        } else {
            let t = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_nanos() as u64)
                .unwrap_or(42);
            t ^ (std::process::id() as u64).rotate_left(32)
        };
        Self::from_seed(seed)
    }

    pub fn genrand_uint32(&mut self) -> u32 {
        let mag01 = [0u32, MATRIX_A];
        if self.mti >= N {
            let mut kk = 0;
            while kk < N - M {
                let y = (self.mt[kk] & UPPER_MASK) | (self.mt[kk + 1] & LOWER_MASK);
                self.mt[kk] = self.mt[kk + M] ^ (y >> 1) ^ mag01[(y & 1) as usize];
                kk += 1;
            }
            while kk < N - 1 {
                let y = (self.mt[kk] & UPPER_MASK) | (self.mt[kk + 1] & LOWER_MASK);
                self.mt[kk] = self.mt[kk + M - N] ^ (y >> 1) ^ mag01[(y & 1) as usize];
                kk += 1;
            }
            let y = (self.mt[N - 1] & UPPER_MASK) | (self.mt[0] & LOWER_MASK);
            self.mt[N - 1] = self.mt[M - 1] ^ (y >> 1) ^ mag01[(y & 1) as usize];
            self.mti = 0;
        }
        let mut y = self.mt[self.mti];
        self.mti += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >> 18;
        y
    }

    /// `random.getrandbits(k)` for 0 < k <= 64.
    pub fn getrandbits(&mut self, k: u32) -> u64 {
        if k == 0 {
            return 0;
        }
        if k <= 32 {
            return (self.genrand_uint32() >> (32 - k)) as u64;
        }
        // little-endian 32-bit words; the last word is truncated
        let mut result: u64 = 0;
        let mut remaining = k;
        let mut shift = 0;
        while remaining > 0 {
            let mut r = self.genrand_uint32();
            if remaining < 32 {
                r >>= 32 - remaining;
            }
            result |= (r as u64) << shift;
            shift += 32;
            remaining = remaining.saturating_sub(32);
        }
        result
    }

    /// `random._randbelow(n)` (getrandbits based).
    pub fn randbelow(&mut self, n: u64) -> u64 {
        if n == 0 {
            return 0;
        }
        let k = 64 - n.leading_zeros();
        let mut r = self.getrandbits(k);
        while r >= n {
            r = self.getrandbits(k);
        }
        r
    }

    /// `random.randint(a, b)` (inclusive bounds).
    pub fn randint(&mut self, a: i64, b: i64) -> i64 {
        let width = (b - a + 1) as u64;
        a + self.randbelow(width) as i64
    }

    /// `random.random()`.
    pub fn random(&mut self) -> f64 {
        let a = (self.genrand_uint32() >> 5) as f64;
        let b = (self.genrand_uint32() >> 6) as f64;
        (a * 67108864.0 + b) * (1.0 / 9007199254740992.0)
    }
}

/// Create the process-wide generator honouring `XTEA_SEED` when set.
pub fn new_rng(seed: Option<u64>) -> PyRandom {
    if let Some(s) = seed {
        return PyRandom::from_seed(s);
    }
    if let Ok(v) = std::env::var("XTEA_SEED") {
        if let Ok(s) = v.trim().parse::<u64>() {
            return PyRandom::from_seed(s);
        }
    }
    PyRandom::from_entropy()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn matches_cpython_seed_42() {
        // python3 -c "import random; random.seed(42); print([random.randint(0, 24) for _ in range(5)], random.random())"
        let mut r = PyRandom::from_seed(42);
        let v: Vec<i64> = (0..5).map(|_| r.randint(0, 24)).collect();
        assert_eq!(v, vec![20, 3, 0, 23, 8]);
    }
}
