//! Emulation of CPython's `set` hash table (Objects/setobject.c, CPython 3.8 - 3.13).
//!
//! xTea frequently iterates over Python `set`s and keeps the first/best element, so the
//! iteration order influences results whenever there are ties.  Python sets of ints and of
//! tuples of ints iterate in a deterministic order that depends only on the hash values and
//! on the sequence of insertions/deletions.  This module reproduces that order exactly so the
//! Rust port can produce byte-identical output.

/// Number of linear probes performed before applying the perturbation (CPython: LINEAR_PROBES).
const LINEAR_PROBES: usize = 9;
/// CPython: PERTURB_SHIFT.
const PERTURB_SHIFT: u32 = 5;
/// CPython: PySet_MINSIZE.
const MINSIZE: usize = 8;

/// Something that has a CPython-compatible hash value.
pub trait PyHash {
    fn py_hash(&self) -> i64;
}

/// CPython hash of an integer (`hash(int)`).
#[inline]
pub fn py_hash_int(v: i64) -> i64 {
    const MODULUS: u64 = (1u64 << 61) - 1;
    let mut h = (v.unsigned_abs() % MODULUS) as i64;
    if v < 0 {
        h = -h;
    }
    if h == -1 {
        h = -2;
    }
    h
}

/// CPython (>= 3.8) hash of a tuple whose items have the given hashes (xxHash based `tuplehash`).
#[inline]
pub fn py_hash_tuple(item_hashes: &[i64]) -> i64 {
    const XXPRIME_1: u64 = 11400714785074694791;
    const XXPRIME_2: u64 = 14029467366897019727;
    const XXPRIME_5: u64 = 2870177450012600261;
    let mut acc: u64 = XXPRIME_5;
    for &h in item_hashes {
        let lane = h as u64;
        acc = acc.wrapping_add(lane.wrapping_mul(XXPRIME_2));
        acc = acc.rotate_left(31);
        acc = acc.wrapping_mul(XXPRIME_1);
    }
    acc = acc.wrapping_add((item_hashes.len() as u64) ^ (XXPRIME_5 ^ 3527539u64));
    if acc == u64::MAX {
        return 1546275796;
    }
    acc as i64
}

impl PyHash for i64 {
    #[inline]
    fn py_hash(&self) -> i64 {
        py_hash_int(*self)
    }
}

impl PyHash for (i64, i64) {
    #[inline]
    fn py_hash(&self) -> i64 {
        py_hash_tuple(&[py_hash_int(self.0), py_hash_int(self.1)])
    }
}

#[derive(Clone, Debug)]
enum Slot<K> {
    Empty,
    Dummy,
    Active(K, i64),
}

/// A set with CPython-identical iteration order.
#[derive(Clone, Debug)]
pub struct PySet<K> {
    table: Vec<Slot<K>>,
    mask: usize,
    fill: usize,
    used: usize,
}

impl<K: PyHash + PartialEq + Clone> Default for PySet<K> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: PyHash + PartialEq + Clone> PySet<K> {
    pub fn new() -> Self {
        PySet {
            table: vec![Slot::Empty; MINSIZE],
            mask: MINSIZE - 1,
            fill: 0,
            used: 0,
        }
    }

    /// `set(iterable)` for a non-set iterable (keys are added one by one, in order).
    pub fn from_iter_ordered<I: IntoIterator<Item = K>>(it: I) -> Self {
        let mut s = Self::new();
        for k in it {
            s.add(k);
        }
        s
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.used
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.used == 0
    }

    /// Iterate in CPython order (table slot order).
    pub fn iter(&self) -> impl Iterator<Item = &K> + '_ {
        self.table.iter().filter_map(|s| match s {
            Slot::Active(k, _) => Some(k),
            _ => None,
        })
    }

    /// Items in iteration order (cloned).
    pub fn to_vec(&self) -> Vec<K> {
        self.iter().cloned().collect()
    }

    fn active_entries(&self) -> impl Iterator<Item = (&K, i64)> + '_ {
        self.table.iter().filter_map(|s| match s {
            Slot::Active(k, h) => Some((k, *h)),
            _ => None,
        })
    }

    /// `set.clear()`: back to an empty minimum-size table.
    pub fn clear(&mut self) {
        self.table = vec![Slot::Empty; MINSIZE];
        self.mask = MINSIZE - 1;
        self.fill = 0;
        self.used = 0;
    }

    #[inline]
    fn slot_matches(slot: &Slot<K>, key: &K, hash: i64) -> bool {
        match slot {
            Slot::Active(k, h) => *h == hash && k == key,
            _ => false,
        }
    }

    /// `set.add(key)` (CPython `set_add_entry`).
    pub fn add(&mut self, key: K) {
        let hash = key.py_hash();
        self.add_with_hash(key, hash);
    }

    fn add_with_hash(&mut self, key: K, hash: i64) {
        let mask = self.mask;
        let mut i = (hash as u64 as usize) & mask;
        let mut entry = i;
        if let Slot::Empty = self.table[entry] {
            self.found_unused(entry, key, hash);
            return;
        }
        let mut freeslot: Option<usize> = None;
        let mut perturb: u64 = hash as u64;
        loop {
            match &self.table[entry] {
                Slot::Active(k, h) => {
                    if *h == hash && *k == key {
                        return;
                    }
                }
                Slot::Dummy => freeslot = Some(entry),
                Slot::Empty => {}
            }
            if i + LINEAR_PROBES <= mask {
                for _ in 0..LINEAR_PROBES {
                    entry += 1;
                    match &self.table[entry] {
                        Slot::Empty => {
                            self.found_unused_or_dummy(entry, freeslot, key, hash);
                            return;
                        }
                        Slot::Active(k, h) => {
                            if *h == hash && *k == key {
                                return;
                            }
                        }
                        Slot::Dummy => freeslot = Some(entry),
                    }
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb as usize)) & mask;
            entry = i;
            if let Slot::Empty = self.table[entry] {
                self.found_unused_or_dummy(entry, freeslot, key, hash);
                return;
            }
        }
    }

    fn found_unused_or_dummy(&mut self, entry: usize, freeslot: Option<usize>, key: K, hash: i64) {
        match freeslot {
            None => self.found_unused(entry, key, hash),
            Some(fs) => {
                self.used += 1;
                self.table[fs] = Slot::Active(key, hash);
            }
        }
    }

    fn found_unused(&mut self, entry: usize, key: K, hash: i64) {
        self.fill += 1;
        self.used += 1;
        self.table[entry] = Slot::Active(key, hash);
        if self.fill * 5 < self.mask * 3 {
            return;
        }
        let minused = if self.used > 50000 { self.used * 2 } else { self.used * 4 };
        self.resize(minused);
    }

    fn insert_clean(table: &mut [Slot<K>], mask: usize, key: K, hash: i64) {
        let mut perturb: u64 = hash as u64;
        let mut i = (hash as u64 as usize) & mask;
        loop {
            if let Slot::Empty = table[i] {
                table[i] = Slot::Active(key, hash);
                return;
            }
            if i + LINEAR_PROBES <= mask {
                for j in 1..=LINEAR_PROBES {
                    if let Slot::Empty = table[i + j] {
                        table[i + j] = Slot::Active(key, hash);
                        return;
                    }
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb as usize)) & mask;
        }
    }

    /// CPython `set_table_resize`.
    fn resize(&mut self, minused: usize) {
        let mut newsize = MINSIZE;
        while newsize <= minused {
            newsize <<= 1;
        }
        if newsize == MINSIZE && self.table.len() == MINSIZE && self.fill == self.used {
            // No dummies, small table: nothing to do.
            return;
        }
        let old = std::mem::replace(&mut self.table, vec![Slot::Empty; newsize]);
        self.mask = newsize - 1;
        self.fill = self.used;
        let mask = self.mask;
        for slot in old.into_iter() {
            if let Slot::Active(k, h) = slot {
                Self::insert_clean(&mut self.table, mask, k, h);
            }
        }
    }

    /// CPython `set_lookkey`: returns the index of the slot holding `key`, if present.
    fn lookkey(&self, key: &K, hash: i64) -> Option<usize> {
        let mask = self.mask;
        let mut i = (hash as u64 as usize) & mask;
        let mut entry = i;
        if let Slot::Empty = self.table[entry] {
            return None;
        }
        let mut perturb: u64 = hash as u64;
        loop {
            if Self::slot_matches(&self.table[entry], key, hash) {
                return Some(entry);
            }
            if i + LINEAR_PROBES <= mask {
                for _ in 0..LINEAR_PROBES {
                    entry += 1;
                    if let Slot::Empty = self.table[entry] {
                        return None;
                    }
                    if Self::slot_matches(&self.table[entry], key, hash) {
                        return Some(entry);
                    }
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb as usize)) & mask;
            entry = i;
            if let Slot::Empty = self.table[entry] {
                return None;
            }
        }
    }

    pub fn contains(&self, key: &K) -> bool {
        self.lookkey(key, key.py_hash()).is_some()
    }

    /// `set.discard(key)`; returns true when the key was present.
    pub fn discard(&mut self, key: &K) -> bool {
        match self.lookkey(key, key.py_hash()) {
            Some(idx) => {
                self.table[idx] = Slot::Dummy;
                self.used -= 1;
                true
            }
            None => false,
        }
    }

    /// `set.remove(key)`; panics (like Python's KeyError) when absent.
    pub fn remove(&mut self, key: &K) {
        if !self.discard(key) {
            panic!("KeyError in PySet::remove");
        }
    }

    /// `set(other)` - copy constructor for a set argument (CPython `set_merge` into an empty set).
    pub fn copy_of(other: &PySet<K>) -> Self {
        let mut s = Self::new();
        s.merge(other);
        s
    }

    /// `self.update(other_set)` / `self |= other_set` (CPython `set_merge`).
    pub fn merge(&mut self, other: &PySet<K>) {
        if other.used == 0 {
            return;
        }
        if (self.fill + other.used) * 5 >= self.mask * 3 {
            self.resize((self.used + other.used) * 2);
        }
        if self.fill == 0 && self.mask == other.mask && other.fill == other.used {
            for (dst, src) in self.table.iter_mut().zip(other.table.iter()) {
                if let Slot::Active(k, h) = src {
                    *dst = Slot::Active(k.clone(), *h);
                }
            }
            self.fill = other.fill;
            self.used = other.used;
            return;
        }
        if self.fill == 0 {
            self.fill = other.used;
            self.used = other.used;
            let mask = self.mask;
            for (k, h) in other.active_entries() {
                Self::insert_clean(&mut self.table, mask, k.clone(), h);
            }
            return;
        }
        let entries: Vec<(K, i64)> = other.active_entries().map(|(k, h)| (k.clone(), h)).collect();
        for (k, h) in entries {
            self.add_with_hash(k, h);
        }
    }

    /// `self.update(iterable)` for a non-set iterable.
    pub fn update_iter<I: IntoIterator<Item = K>>(&mut self, it: I) {
        for k in it {
            self.add(k);
        }
    }

    /// `self -= other_set` (CPython `set_difference_update_internal`).
    pub fn difference_update(&mut self, other: &PySet<K>) {
        let keys: Vec<(K, i64)> = other.active_entries().map(|(k, h)| (k.clone(), h)).collect();
        for (k, h) in keys {
            if let Some(idx) = self.lookkey(&k, h) {
                self.table[idx] = Slot::Dummy;
                self.used -= 1;
            }
        }
        if self.fill - self.used <= self.mask / 4 {
            return;
        }
        let minused = if self.used > 50000 { self.used * 2 } else { self.used * 4 };
        self.resize(minused);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn int_hash() {
        assert_eq!(py_hash_int(5), 5);
        assert_eq!(py_hash_int(-1), -2);
        assert_eq!(py_hash_int(-5), -5);
    }

    #[test]
    fn tuple_hash_matches_cpython() {
        // Values obtained from CPython 3.8: hash((1, 2)), hash((0, 0)), hash(())
        assert_eq!(py_hash_tuple(&[py_hash_int(1), py_hash_int(2)]), -3550055125485641917);
        assert_eq!(py_hash_tuple(&[]), 5740354900026072187);
    }

    #[test]
    fn small_int_order() {
        let mut s: PySet<i64> = PySet::new();
        for v in [1007, 1008, 1015] {
            s.add(v);
        }
        // CPython: list({1007, 1008, 1015}) built in that order
        let got: Vec<i64> = s.to_vec();
        assert_eq!(got, vec![1008, 1015, 1007]);
    }
}

/// `frexp`: v = m * 2^e with 0.5 <= |m| < 1 (m = 0, e = 0 for v = 0).
fn frexp(v: f64) -> (f64, i32) {
    if v == 0.0 || !v.is_finite() {
        return (v, 0);
    }
    let bits = v.to_bits();
    let exp = ((bits >> 52) & 0x7ff) as i32;
    if exp == 0 {
        // subnormal: scale up first
        let (m, e) = frexp(v * f64::from_bits(0x4350000000000000)); // 2^54
        return (m, e - 54);
    }
    let m = f64::from_bits((bits & !(0x7ffu64 << 52)) | (1022u64 << 52));
    (m, exp - 1022)
}

/// CPython hash of a float (`_Py_HashDouble`); equals `py_hash_int` for integral values.
pub fn py_hash_float(v: f64) -> i64 {
    const BITS: i32 = 61;
    const MODULUS: u64 = (1u64 << 61) - 1;
    if !v.is_finite() {
        if v.is_infinite() {
            return if v > 0.0 { 314159 } else { -314159 };
        }
        return 0;
    }
    let (mut m, mut e) = frexp(v);
    let mut sign: i64 = 1;
    if m < 0.0 {
        sign = -1;
        m = -m;
    }
    let mut x: u64 = 0;
    while m != 0.0 {
        x = ((x << 28) & MODULUS) | (x >> (BITS - 28));
        m *= 268435456.0;
        e -= 28;
        let y = m as u64;
        m -= y as f64;
        x += y;
        if x >= MODULUS {
            x -= MODULUS;
        }
    }
    let e = if e >= 0 { e % BITS } else { BITS - 1 - ((-1 - e) % BITS) };
    x = ((x << e) & MODULUS) | (x >> (BITS - e));
    let mut h = (x as i64).wrapping_mul(sign);
    if h == -1 {
        h = -2;
    }
    h
}

#[cfg(test)]
mod float_hash_tests {
    use super::*;
    #[test]
    fn matches_cpython() {
        // python3 -c 'print(hash(v), hash((v, v+1)))'
        let cases: [(f64, i64, i64); 9] = [
            (0.0, 0, -1950498447580522560),
            (1.5, 1152921504606846977, -1826667736619738083),
            (-1.5, -1152921504606846977, 6303917944917004362),
            (412.763, 1759358216029995420, -4589340060353359633),
            (123456.5, 1152921504606970432, -2883518688192945867),
            (10000000000.25, 576460762303423488, 5010364700553156067),
            (3.14159, 326484311674566659, -665526534675067629),
            (0.1, 230584300921369408, 7298586803571078865),
            (99999999.987, 2275867058165375231, -576834598107756429),
        ];
        for (v, h, ht) in cases {
            assert_eq!(py_hash_float(v), h, "hash({})", v);
            assert_eq!(py_hash_tuple(&[py_hash_float(v), py_hash_float(v + 1.0)]), ht, "hash(({}, {}))", v, v + 1.0);
        }
        for i in [-5i64, 0, 1, 7, 123456789, 1 << 40] {
            assert_eq!(py_hash_float(i as f64), py_hash_int(i));
        }
    }
}
