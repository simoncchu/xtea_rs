//! Faithful port of the `intervaltree` package vendored by xTea (intervaltree.py / node.py /
//! interval.py, version 3.x by Chaim Leib Halbert).
//!
//! Only the operations used by xTea are provided: `addi(begin, end)` and point queries
//! (`tree[pos]`).  The tree structure (AVL balancing with interval promotion) and the
//! iteration order of the CPython sets holding the intervals are reproduced exactly, so that
//! "take the first hit" logic in xTea returns the same interval as the Python implementation.

use crate::py::set::{py_hash_float, py_hash_tuple};
use std::collections::HashMap;

const NIL: u32 = u32::MAX;
const EMPTY: u32 = u32::MAX;
const DUMMY: u32 = u32::MAX - 1;
const LINEAR_PROBES: usize = 9;
const PERTURB_SHIFT: u32 = 5;
const MINSIZE: usize = 8;

/// A CPython-ordered set of interval ids (hashes looked up in the owning tree).
#[derive(Clone, Debug)]
struct IvSet {
    table: Box<[u32]>,
    fill: u32,
    used: u32,
}

impl IvSet {
    fn new() -> Self {
        IvSet { table: vec![EMPTY; MINSIZE].into_boxed_slice(), fill: 0, used: 0 }
    }
    #[inline]
    fn mask(&self) -> usize {
        self.table.len() - 1
    }
    fn len(&self) -> usize {
        self.used as usize
    }
    fn is_empty(&self) -> bool {
        self.used == 0
    }
    fn iter(&self) -> impl Iterator<Item = u32> + '_ {
        self.table.iter().copied().filter(|&v| v != EMPTY && v != DUMMY)
    }
    fn to_vec(&self) -> Vec<u32> {
        self.iter().collect()
    }

    fn add(&mut self, id: u32, hashes: &[i64]) {
        let hash = hashes[id as usize];
        let mask = self.mask();
        let mut i = (hash as u64 as usize) & mask;
        let mut entry = i;
        if self.table[entry] == EMPTY {
            self.found_unused(entry, id, hashes);
            return;
        }
        let mut freeslot: Option<usize> = None;
        let mut perturb: u64 = hash as u64;
        loop {
            let v = self.table[entry];
            if v == DUMMY {
                freeslot = Some(entry);
            } else if v == id {
                return;
            }
            if i + LINEAR_PROBES <= mask {
                for _ in 0..LINEAR_PROBES {
                    entry += 1;
                    let v = self.table[entry];
                    if v == EMPTY {
                        self.found_unused_or_dummy(entry, freeslot, id, hashes);
                        return;
                    }
                    if v == DUMMY {
                        freeslot = Some(entry);
                    } else if v == id {
                        return;
                    }
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb as usize)) & mask;
            entry = i;
            if self.table[entry] == EMPTY {
                self.found_unused_or_dummy(entry, freeslot, id, hashes);
                return;
            }
        }
    }

    fn found_unused_or_dummy(&mut self, entry: usize, freeslot: Option<usize>, id: u32, hashes: &[i64]) {
        match freeslot {
            None => self.found_unused(entry, id, hashes),
            Some(fs) => {
                self.used += 1;
                self.table[fs] = id;
            }
        }
    }

    fn found_unused(&mut self, entry: usize, id: u32, hashes: &[i64]) {
        self.fill += 1;
        self.used += 1;
        self.table[entry] = id;
        let mask = self.mask();
        if (self.fill as usize) * 5 < mask * 3 {
            return;
        }
        let used = self.used as usize;
        let minused = if used > 50000 { used * 2 } else { used * 4 };
        self.resize(minused, hashes);
    }

    fn insert_clean(table: &mut [u32], mask: usize, id: u32, hash: i64) {
        let mut perturb: u64 = hash as u64;
        let mut i = (hash as u64 as usize) & mask;
        loop {
            if table[i] == EMPTY {
                table[i] = id;
                return;
            }
            if i + LINEAR_PROBES <= mask {
                for j in 1..=LINEAR_PROBES {
                    if table[i + j] == EMPTY {
                        table[i + j] = id;
                        return;
                    }
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb as usize)) & mask;
        }
    }

    fn resize(&mut self, minused: usize, hashes: &[i64]) {
        let mut newsize = MINSIZE;
        while newsize <= minused {
            newsize <<= 1;
        }
        if newsize == MINSIZE && self.table.len() == MINSIZE && self.fill == self.used {
            return;
        }
        let old = std::mem::replace(&mut self.table, vec![EMPTY; newsize].into_boxed_slice());
        self.fill = self.used;
        let mask = newsize - 1;
        for &v in old.iter() {
            if v != EMPTY && v != DUMMY {
                Self::insert_clean(&mut self.table, mask, v, hashes[v as usize]);
            }
        }
    }

    fn lookkey(&self, id: u32, hashes: &[i64]) -> Option<usize> {
        let hash = hashes[id as usize];
        let mask = self.mask();
        let mut i = (hash as u64 as usize) & mask;
        let mut entry = i;
        if self.table[entry] == EMPTY {
            return None;
        }
        let mut perturb: u64 = hash as u64;
        loop {
            if self.table[entry] == id {
                return Some(entry);
            }
            if i + LINEAR_PROBES <= mask {
                for _ in 0..LINEAR_PROBES {
                    entry += 1;
                    let v = self.table[entry];
                    if v == EMPTY {
                        return None;
                    }
                    if v == id {
                        return Some(entry);
                    }
                }
            }
            perturb >>= PERTURB_SHIFT;
            i = (i.wrapping_mul(5).wrapping_add(1).wrapping_add(perturb as usize)) & mask;
            entry = i;
            if self.table[entry] == EMPTY {
                return None;
            }
        }
    }

    fn contains(&self, id: u32, hashes: &[i64]) -> bool {
        self.lookkey(id, hashes).is_some()
    }

    fn discard(&mut self, id: u32, hashes: &[i64]) -> bool {
        match self.lookkey(id, hashes) {
            Some(idx) => {
                self.table[idx] = DUMMY;
                self.used -= 1;
                true
            }
            None => false,
        }
    }

    /// `set(other)` copy (CPython set_merge into a fresh set).
    fn copy_of(other: &IvSet, hashes: &[i64]) -> IvSet {
        let mut s = IvSet::new();
        if other.used == 0 {
            return s;
        }
        if (s.fill + other.used) as usize * 5 >= s.mask() * 3 {
            s.resize(((s.used + other.used) * 2) as usize, hashes);
        }
        if s.fill == 0 && s.table.len() == other.table.len() && other.fill == other.used {
            s.table.copy_from_slice(&other.table);
            s.fill = other.fill;
            s.used = other.used;
            return s;
        }
        // s.fill == 0 here
        s.fill = other.used;
        s.used = other.used;
        let mask = s.mask();
        for v in other.iter() {
            Self::insert_clean(&mut s.table, mask, v, hashes[v as usize]);
        }
        s
    }

    /// `self -= other`
    fn difference_update(&mut self, other: &IvSet, hashes: &[i64]) {
        let ids: Vec<u32> = other.iter().collect();
        for id in ids {
            if let Some(idx) = self.lookkey(id, hashes) {
                self.table[idx] = DUMMY;
                self.used -= 1;
            }
        }
        if (self.fill - self.used) as usize <= self.mask() / 4 {
            return;
        }
        let used = self.used as usize;
        let minused = if used > 50000 { used * 2 } else { used * 4 };
        self.resize(minused, hashes);
    }
}

#[derive(Clone, Debug)]
struct Node {
    x_center: f64,
    s_center: IvSet,
    left: u32,
    right: u32,
    depth: i32,
    balance: i32,
}

/// Interval tree keyed by half-open intervals `[begin, end)`.  Coordinates are stored as f64
/// so that the float intervals the original sometimes builds (e.g. an annotation extended by a
/// mean insert size) are handled; integer intervals behave exactly as Python ints do (CPython
/// hashes integral floats like the equal ints).
#[derive(Clone, Debug, Default)]
pub struct IntervalTree {
    nodes: Vec<Node>,
    begins: Vec<f64>,
    ends: Vec<f64>,
    hashes: Vec<i64>,
    index: HashMap<(u64, u64), u32>,
    top: Option<u32>,
}

#[inline]
fn fkey(v: f64) -> u64 {
    // +0.0 and -0.0 are the same key (equal in Python)
    if v == 0.0 {
        0
    } else {
        v.to_bits()
    }
}

impl IntervalTree {
    pub fn new() -> Self {
        IntervalTree {
            nodes: Vec::new(),
            begins: Vec::new(),
            ends: Vec::new(),
            hashes: Vec::new(),
            index: HashMap::new(),
            top: None,
        }
    }

    pub fn len(&self) -> usize {
        self.begins.len()
    }

    pub fn is_empty(&self) -> bool {
        self.begins.is_empty()
    }

    #[inline]
    fn get(&self, n: u32, right: bool) -> u32 {
        if right {
            self.nodes[n as usize].right
        } else {
            self.nodes[n as usize].left
        }
    }

    #[inline]
    fn set(&mut self, n: u32, right: bool, v: u32) {
        if right {
            self.nodes[n as usize].right = v;
        } else {
            self.nodes[n as usize].left = v;
        }
    }

    #[inline]
    fn contains_point(&self, iv: u32, p: f64) -> bool {
        self.begins[iv as usize] <= p && p < self.ends[iv as usize]
    }

    #[inline]
    fn center_hit(&self, n: u32, iv: u32) -> bool {
        self.contains_point(iv, self.nodes[n as usize].x_center)
    }

    #[inline]
    fn hit_branch(&self, n: u32, iv: u32) -> bool {
        self.begins[iv as usize] > self.nodes[n as usize].x_center
    }

    fn new_node<I: IntoIterator<Item = u32>>(&mut self, x_center: f64, s_center: I) -> u32 {
        let mut set = IvSet::new();
        for id in s_center {
            set.add(id, &self.hashes);
        }
        let node = Node { x_center, s_center: set, left: NIL, right: NIL, depth: 0, balance: 0 };
        self.nodes.push(node);
        let idx = (self.nodes.len() - 1) as u32;
        // Node.__init__ calls self.rotate(); without children this only refreshes the balance.
        self.rotate(idx);
        idx
    }

    fn refresh_balance(&mut self, n: u32) {
        let l = self.nodes[n as usize].left;
        let r = self.nodes[n as usize].right;
        let ld = if l != NIL { self.nodes[l as usize].depth } else { 0 };
        let rd = if r != NIL { self.nodes[r as usize].depth } else { 0 };
        let node = &mut self.nodes[n as usize];
        node.depth = 1 + ld.max(rd);
        node.balance = rd - ld;
    }

    fn rotate(&mut self, n: u32) -> u32 {
        self.refresh_balance(n);
        let bal = self.nodes[n as usize].balance;
        if bal.abs() < 2 {
            return n;
        }
        let my_heavy = bal > 0;
        let child = self.get(n, my_heavy);
        let child_bal = self.nodes[child as usize].balance;
        let child_heavy = child_bal > 0;
        if my_heavy == child_heavy || child_bal == 0 {
            self.srotate(n)
        } else {
            self.drotate(n)
        }
    }

    fn srotate(&mut self, n: u32) -> u32 {
        let heavy = self.nodes[n as usize].balance > 0;
        let light = !heavy;
        let save = self.get(n, heavy);
        let save_light = self.get(save, light);
        self.set(n, heavy, save_light);
        let rotated = self.rotate(n);
        self.set(save, light, rotated);

        // Some intervals may overlap both self.x_center and save.x_center: promote them.
        let sl = self.get(save, light);
        let save_x = self.nodes[save as usize].x_center;
        let promotees: Vec<u32> = self.nodes[sl as usize]
            .s_center
            .iter()
            .filter(|&iv| self.contains_point(iv, save_x))
            .collect();
        if !promotees.is_empty() {
            for &iv in &promotees {
                let cur = self.get(save, light);
                let newn = self.remove(cur, iv);
                self.set(save, light, newn);
            }
            let hashes = std::mem::take(&mut self.hashes);
            for &iv in &promotees {
                self.nodes[save as usize].s_center.add(iv, &hashes);
            }
            self.hashes = hashes;
        }
        self.refresh_balance(save);
        save
    }

    fn drotate(&mut self, n: u32) -> u32 {
        let my_heavy = self.nodes[n as usize].balance > 0;
        let child = self.get(n, my_heavy);
        let new_child = self.srotate(child);
        self.set(n, my_heavy, new_child);
        self.refresh_balance(n);
        self.srotate(n)
    }

    fn node_add(&mut self, n: u32, iv: u32) -> u32 {
        if self.center_hit(n, iv) {
            let hashes = std::mem::take(&mut self.hashes);
            self.nodes[n as usize].s_center.add(iv, &hashes);
            self.hashes = hashes;
            n
        } else {
            let direction = self.hit_branch(n, iv);
            let child = self.get(n, direction);
            if child == NIL {
                let b = self.begins[iv as usize];
                let newn = self.new_node(b, [iv]);
                self.set(n, direction, newn);
                self.refresh_balance(n);
                n
            } else {
                let newc = self.node_add(child, iv);
                self.set(n, direction, newc);
                self.rotate(n)
            }
        }
    }

    /// Node.remove (raises if absent).
    fn remove(&mut self, n: u32, iv: u32) -> u32 {
        let mut done = false;
        self.remove_helper(n, iv, &mut done, true)
    }

    fn remove_helper(&mut self, n: u32, iv: u32, done: &mut bool, should_raise: bool) -> u32 {
        if self.center_hit(n, iv) {
            if !should_raise && !self.nodes[n as usize].s_center.contains(iv, &self.hashes) {
                *done = true;
                return n;
            }
            let hashes = std::mem::take(&mut self.hashes);
            let removed = self.nodes[n as usize].s_center.discard(iv, &hashes);
            self.hashes = hashes;
            if !removed {
                panic!("intervaltree: KeyError while removing interval");
            }
            if !self.nodes[n as usize].s_center.is_empty() {
                *done = true;
                return n;
            }
            self.prune(n)
        } else {
            let direction = self.hit_branch(n, iv);
            let child = self.get(n, direction);
            if child == NIL {
                if should_raise {
                    panic!("intervaltree: ValueError while removing interval");
                }
                *done = true;
                return n;
            }
            let newc = self.remove_helper(child, iv, done, should_raise);
            self.set(n, direction, newc);
            if !*done {
                return self.rotate(n);
            }
            n
        }
    }

    /// Node.prune: returns the new subtree root (may be NIL).
    fn prune(&mut self, n: u32) -> u32 {
        let l = self.nodes[n as usize].left;
        let r = self.nodes[n as usize].right;
        if l == NIL || r == NIL {
            // direction = not self[0]
            let direction = l == NIL;
            return self.get(n, direction);
        }
        let (heir, new_left) = self.pop_greatest_child(l);
        self.set(n, false, new_left);
        let nl = self.nodes[n as usize].left;
        let nr = self.nodes[n as usize].right;
        self.set(heir, false, nl);
        self.set(heir, true, nr);
        self.refresh_balance(heir);
        self.rotate(heir)
    }

    /// Node.pop_greatest_child: returns (greatest_child, new_subtree_root).
    fn pop_greatest_child(&mut self, n: u32) -> (u32, u32) {
        if self.nodes[n as usize].right == NIL {
            // sorted(self.s_center, key=attrgetter('end', 'begin'))
            let mut ivs: Vec<u32> = self.nodes[n as usize].s_center.to_vec();
            ivs.sort_by(|&a, &b| {
                (self.ends[a as usize], self.begins[a as usize])
                    .partial_cmp(&(self.ends[b as usize], self.begins[b as usize]))
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let max_iv = ivs.pop().expect("pop_greatest_child on empty node");
            let mut new_x_center = self.nodes[n as usize].x_center;
            while let Some(next_max_iv) = ivs.pop() {
                if self.ends[next_max_iv as usize] == self.ends[max_iv as usize] {
                    continue;
                }
                // Python max(): keeps the first argument on ties
                let cand = self.ends[next_max_iv as usize];
                if cand > new_x_center {
                    new_x_center = cand;
                }
            }
            let members: Vec<u32> = self.nodes[n as usize]
                .s_center
                .iter()
                .filter(|&iv| self.contains_point(iv, new_x_center))
                .collect();
            let child = self.new_node(new_x_center, members);
            let hashes = std::mem::take(&mut self.hashes);
            let child_set = self.nodes[child as usize].s_center.clone();
            self.nodes[n as usize].s_center.difference_update(&child_set, &hashes);
            self.hashes = hashes;
            if !self.nodes[n as usize].s_center.is_empty() {
                (child, n)
            } else {
                (child, self.nodes[n as usize].left)
            }
        } else {
            let r = self.nodes[n as usize].right;
            let (greatest_child, new_r) = self.pop_greatest_child(r);
            self.set(n, true, new_r);
            // Move any overlaps into greatest_child
            let copy = IvSet::copy_of(&self.nodes[n as usize].s_center, &self.hashes);
            let gx = self.nodes[greatest_child as usize].x_center;
            for iv in copy.iter() {
                if self.contains_point(iv, gx) {
                    let hashes = std::mem::take(&mut self.hashes);
                    let ok = self.nodes[n as usize].s_center.discard(iv, &hashes);
                    self.hashes = hashes;
                    if !ok {
                        panic!("intervaltree: KeyError in pop_greatest_child");
                    }
                    // the return value of greatest_child.add(iv) is ignored in the original
                    let _ = self.node_add(greatest_child, iv);
                }
            }
            if !self.nodes[n as usize].s_center.is_empty() {
                self.refresh_balance(n);
                let new_self = self.rotate(n);
                (greatest_child, new_self)
            } else {
                let new_self = self.prune(n);
                (greatest_child, new_self)
            }
        }
    }

    /// `IntervalTree.addi(begin, end)` (data = None) for integer coordinates.
    pub fn addi(&mut self, begin: i64, end: i64) {
        self.addf(begin as f64, end as f64)
    }

    /// `IntervalTree.addi(begin, end)` with float coordinates.
    pub fn addf(&mut self, begin: f64, end: f64) {
        if self.index.contains_key(&(fkey(begin), fkey(end))) {
            return;
        }
        if begin >= end {
            panic!("IntervalTree: Null Interval objects not allowed in IntervalTree: Interval({}, {})", begin, end);
        }
        let id = self.begins.len() as u32;
        self.begins.push(begin);
        self.ends.push(end);
        self.hashes.push(py_hash_tuple(&[py_hash_float(begin), py_hash_float(end)]));
        self.index.insert((fkey(begin), fkey(end)), id);
        match self.top {
            None => {
                let n = self.new_node(begin, [id]);
                self.top = Some(n);
            }
            Some(t) => {
                let n = self.node_add(t, id);
                self.top = Some(n);
            }
        }
    }

    /// `tree[p]` / `tree.at(p)`: all intervals containing `p`, in CPython set iteration order
    /// (integer coordinates).
    pub fn at(&self, p: i64) -> Vec<(i64, i64)> {
        self.at_f(p as f64).into_iter().map(|(b, e)| (b as i64, e as i64)).collect()
    }

    /// `tree[p]` with float coordinates.
    pub fn at_f(&self, p: f64) -> Vec<(f64, f64)> {
        let mut result = IvSet::new();
        let mut cur = match self.top {
            None => return Vec::new(),
            Some(t) => t,
        };
        // the result set is grown with `add`, which needs mutable hashes access only for reads
        loop {
            let node = &self.nodes[cur as usize];
            for iv in node.s_center.iter() {
                if self.contains_point(iv, p) {
                    result.add(iv, &self.hashes);
                }
            }
            if p < node.x_center && node.left != NIL {
                cur = node.left;
            } else if p > node.x_center && node.right != NIL {
                cur = node.right;
            } else {
                break;
            }
        }
        result.iter().map(|iv| (self.begins[iv as usize], self.ends[iv as usize])).collect()
    }

    /// Number of intervals containing `p` (order independent).
    pub fn count_at(&self, p: i64) -> usize {
        self.at(p).len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_point_queries() {
        let mut t = IntervalTree::new();
        t.addi(10, 20);
        t.addi(15, 25);
        t.addi(30, 40);
        let mut hits = t.at(16);
        hits.sort();
        assert_eq!(hits, vec![(10, 20), (15, 25)]);
        assert!(t.at(25).is_empty());
        assert_eq!(t.at(30), vec![(30, 40)]);
    }
}
