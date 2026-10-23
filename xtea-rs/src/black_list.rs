//! Port of `x_black_list.py`.

use crate::annotation::process_chrm_name;
use crate::intervaltree::IntervalTree;
use crate::utils::{self, split_ws};
use anyhow::Result;
use std::collections::HashMap;

pub struct XBlackList {
    pub b_with_chr: bool,
    pub m_interval_tree: HashMap<String, IntervalTree>,
}

impl Default for XBlackList {
    fn default() -> Self {
        Self::new()
    }
}

impl XBlackList {
    pub fn new() -> Self {
        XBlackList { b_with_chr: true, m_interval_tree: HashMap::new() }
    }

    pub fn set_with_chr(&mut self, b: bool) {
        self.b_with_chr = b;
    }

    /// `load_index_regions`
    pub fn load_index_regions(&mut self, sf_bed: &str) -> Result<()> {
        if !utils::is_file(sf_bed) {
            return Ok(());
        }
        for line in utils::read_lines(sf_bed)? {
            let fields = split_ws(&line);
            let chrm = fields[0].to_string();
            let istart = utils::pint(fields[1])?;
            let iend = utils::pint(fields[2])?;
            self.m_interval_tree.entry(chrm).or_insert_with(IntervalTree::new).addi(istart, iend);
        }
        Ok(())
    }

    /// `fall_in_region` -> (hit, start of the first interval in set order)
    pub fn fall_in_region(&self, ins_chrm1: &str, ins_pos: i64) -> (bool, i64) {
        let chrm = process_chrm_name(ins_chrm1, self.b_with_chr);
        match self.m_interval_tree.get(&chrm) {
            None => (false, -1),
            Some(t) => match t.at(ins_pos).first() {
                Some((b, _)) => (true, *b),
                None => (false, -1),
            },
        }
    }
}
