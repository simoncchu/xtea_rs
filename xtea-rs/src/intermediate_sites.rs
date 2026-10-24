//! Port of `x_intermediate_sites.py` (XIntemediateSites).

use crate::chrom::is_decoy_contig_chrms_sites;
use crate::globals as gv;
use crate::num::PyVal;
use crate::py::fmt::np_round;
use crate::py::numpy::np_std_i64;
use crate::py::set::PySet;
use crate::utils::{self, split_ws};
use anyhow::Result;
use indexmap::IndexMap;
use std::io::Write;

/// chrm -> pos -> values
pub type SiteMap<T> = IndexMap<String, IndexMap<i64, T>>;

/// `parse_sites_with_clip_cutoff_for_chrm`
pub fn parse_sites_with_clip_cutoff_for_chrm(
    m_clip_pos_freq: &IndexMap<i64, Vec<i64>>,
    cutoff_left_clip: i64,
    cutoff_right_clip: i64,
    cutoff_clip_mate_in_rep: i64,
) -> IndexMap<i64, Vec<PyVal>> {
    let mut out = IndexMap::new();
    for (&pos, v) in m_clip_pos_freq.iter() {
        let (mut nl, mut nr, mut nm) = (0i64, 0i64, 0i64);
        for i in -gv::NEARBY_REGION..gv::NEARBY_REGION {
            if let Some(t) = m_clip_pos_freq.get(&(pos + i)) {
                nl += t[0];
                nr += t[1];
                nm += t[2] + t[3] + t[4];
            }
        }
        if (nl >= cutoff_left_clip || nr >= cutoff_right_clip) && nm >= cutoff_clip_mate_in_rep {
            out.insert(pos, vec![PyVal::Int(v[0]), PyVal::Int(v[1]), PyVal::Int(v[2])]);
        }
    }
    out
}

/// `parse_sites_with_clip_cutoff_for_chrm_with_polyA`
pub fn parse_sites_with_clip_cutoff_for_chrm_with_polya(
    m_clip_pos_freq: &IndexMap<i64, Vec<i64>>,
    cutoff_left_clip: i64,
    cutoff_right_clip: i64,
    cutoff_clip_mate_in_rep: i64,
    cutoff_polya: i64,
) -> IndexMap<i64, Vec<PyVal>> {
    let mut out = IndexMap::new();
    for (&pos, v) in m_clip_pos_freq.iter() {
        let (mut nl, mut nr, mut nm) = (0i64, 0i64, 0i64);
        for i in -gv::NEARBY_REGION..gv::NEARBY_REGION {
            if let Some(t) = m_clip_pos_freq.get(&(pos + i)) {
                nl += t[0];
                nr += t[1];
                nm += t[2] + t[3] + t[4];
            }
        }
        if (nl >= cutoff_left_clip || nr >= cutoff_right_clip) && nm >= cutoff_clip_mate_in_rep {
            let lp = v[5];
            let rp = v[6];
            if lp + rp < cutoff_polya {
                continue;
            }
            out.insert(pos, vec![PyVal::Int(v[0]), PyVal::Int(v[1]), PyVal::Int(v[2]), PyVal::Int(lp), PyVal::Int(rp)]);
        }
    }
    out
}

/// `output_candidate_sites` (skips decoy chromosomes, `XIntemediateSites` rule).
pub fn output_candidate_sites(m: &SiteMap<Vec<PyVal>>, sf_out: &str) -> Result<()> {
    let mut f = utils::create(sf_out)?;
    for (chrm, mp) in m.iter() {
        if is_decoy_contig_chrms_sites(chrm) {
            continue;
        }
        for (pos, vals) in mp.iter() {
            write!(f, "{}\t{}\t", chrm, pos)?;
            for v in vals {
                write!(f, "{}\t", v.py_str())?;
            }
            writeln!(f)?;
        }
    }
    f.flush()?;
    Ok(())
}

/// `is_decoy_contig_chrms` of XIntemediateSites
pub fn is_decoy_contig_chrms(chrm: &str) -> bool {
    is_decoy_contig_chrms_sites(chrm)
}

/// `load_in_candidate_list` (integer values).
pub fn load_in_candidate_list(sf: &str) -> Result<SiteMap<Vec<i64>>> {
    let mut m: SiteMap<Vec<i64>> = IndexMap::new();
    for line in utils::read_lines(sf)? {
        let fields = split_ws(&line);
        if fields.len() < 3 {
            println!("{:?} does not have enough fields", fields);
            continue;
        }
        let chrm = fields[0].to_string();
        let pos = utils::pint(fields[1])?;
        let e = m.entry(chrm).or_default().entry(pos).or_default();
        for v in &fields[2..] {
            e.push(utils::pint(v)?);
        }
    }
    Ok(m)
}

/// `load_in_candidate_list_one_line`
pub fn load_in_candidate_list_one_line(sf: &str) -> Result<SiteMap<String>> {
    let mut m: SiteMap<String> = IndexMap::new();
    for line in utils::read_lines(sf)? {
        let fields = split_ws(&line);
        if fields.len() < 3 {
            println!("{:?} does not have enough fields", fields);
            continue;
        }
        let pos = utils::pint(fields[1])?;
        m.entry(fields[0].to_string()).or_default().insert(pos, utils::rstrip(&line).to_string());
    }
    Ok(m)
}

/// `load_in_candidate_list_str_version`
pub fn load_in_candidate_list_str_version(sf: &str) -> Result<SiteMap<Vec<String>>> {
    let mut m: SiteMap<Vec<String>> = IndexMap::new();
    for line in utils::read_lines(sf)? {
        let fields = split_ws(&line);
        if fields.len() < 3 {
            println!("{:?} does not have enough fields", fields);
            continue;
        }
        let pos = utils::pint(fields[1])?;
        let e = m.entry(fields[0].to_string()).or_default().entry(pos).or_default();
        for v in &fields[2..] {
            e.push(v.to_string());
        }
    }
    Ok(m)
}

/// `load_in_candidate_list2` (string values, no field-count check)
pub fn load_in_candidate_list2(sf: &str) -> Result<SiteMap<Vec<String>>> {
    let mut m: SiteMap<Vec<String>> = IndexMap::new();
    for line in utils::read_lines(sf)? {
        let fields = split_ws(&line);
        let pos = utils::pint(fields[1])?;
        let e = m.entry(fields[0].to_string()).or_default().entry(pos).or_default();
        for v in &fields[2..] {
            e.push(v.to_string());
        }
    }
    Ok(m)
}

/// `is_in_existing_list` (checks `i + ins_pos`, as the original does)
pub fn is_in_existing_list<T>(ins_chrm: &str, ins_pos: i64, m_sites: &SiteMap<T>, i_win: i64) -> bool {
    let m = match m_sites.get(ins_chrm) {
        Some(m) => m,
        None => return false,
    };
    for i in (ins_pos - i_win)..(ins_pos + i_win) {
        if m.contains_key(&(i + ins_pos)) {
            return true;
        }
    }
    false
}

/// `calc_std_derivation`: round(np.std(positions), 2) or -1 for an empty list.
pub fn calc_std_derivation(l_pos: &[i64]) -> PyVal {
    if l_pos.is_empty() {
        return PyVal::Int(-1);
    }
    PyVal::Float(np_round(np_std_i64(l_pos), 2))
}

/// Generic peak calling over sorted positions; clusters are split when the gap is larger than
/// `peak_window`.  `pick` receives the cluster as a CPython-ordered set.
fn cluster_positions<F: FnMut(&PySet<i64>)>(l_pos_sorted: &[i64], peak_window: i64, mut on_cluster: F) {
    let mut pre_pos: i64 = -1;
    let mut set_cluster: PySet<i64> = PySet::new();
    for &pos in l_pos_sorted {
        if pre_pos == -1 {
            pre_pos = pos;
            set_cluster.add(pre_pos);
            continue;
        }
        if pos - pre_pos > peak_window {
            on_cluster(&set_cluster);
            set_cluster.clear();
        }
        pre_pos = pos;
        set_cluster.add(pre_pos);
    }
    on_cluster(&set_cluster);
}

/// `call_peak_candidate_sites` (values: [left_clip, right_clip, ...] as ints)
pub fn call_peak_candidate_sites(m: &SiteMap<Vec<i64>>, peak_window: i64) -> SiteMap<Vec<PyVal>> {
    let mut out: SiteMap<Vec<PyVal>> = IndexMap::new();
    for (chrm, mp) in m.iter() {
        let mut l_pos: Vec<i64> = mp.keys().cloned().collect();
        l_pos.sort();
        let mut res: Vec<(i64, i64)> = Vec::new();
        cluster_positions(&l_pos, peak_window, |set| {
            let mut max_clip = 0;
            let mut cand = 0;
            for &p in set.iter() {
                let v = &mp[&p];
                let all = v[0] + v[1];
                if max_clip < all {
                    cand = p;
                    max_clip = all;
                }
            }
            res.push((cand, max_clip));
        });
        let e = out.entry(chrm.clone()).or_default();
        for (cand, max_clip) in res {
            e.entry(cand).or_insert_with(|| vec![PyVal::Int(max_clip)]);
        }
    }
    out
}

/// `call_peak_candidate_sites_with_std_derivation`
pub fn call_peak_candidate_sites_with_std_derivation(m: &SiteMap<Vec<i64>>, peak_window: i64) -> SiteMap<Vec<PyVal>> {
    let mut out: SiteMap<Vec<PyVal>> = IndexMap::new();
    for (chrm, mp) in m.iter() {
        let mut l_pos: Vec<i64> = mp.keys().cloned().collect();
        l_pos.sort();
        let mut res: Vec<(i64, Vec<PyVal>)> = Vec::new();
        cluster_positions(&l_pos, peak_window, |set| {
            let mut max_clip = 0;
            let mut cand = 0;
            let mut l_lclip: Vec<i64> = Vec::new();
            let mut l_rclip: Vec<i64> = Vec::new();
            for &p in set.iter() {
                let v = &mp[&p];
                let tl = v[0];
                let tr = v[1];
                let all = tl + tr;
                if max_clip < all {
                    cand = p;
                    max_clip = all;
                }
                for _ in 0..tl.max(0) {
                    l_lclip.push(p);
                }
                for _ in 0..tr.max(0) {
                    l_rclip.push(p);
                }
            }
            let f_l = calc_std_derivation(&l_lclip);
            let f_r = calc_std_derivation(&l_rclip);
            l_rclip.extend(l_lclip.iter().cloned());
            let f_all = calc_std_derivation(&l_rclip);
            res.push((cand, vec![PyVal::Int(max_clip), f_l, f_r, f_all]));
        });
        let e = out.entry(chrm.clone()).or_default();
        for (cand, v) in res {
            e.entry(cand).or_insert(v);
        }
    }
    out
}

/// `call_peak_candidate_sites_all_features` (values are (sum_clip, sum_all))
pub fn call_peak_candidate_sites_all_features(m: &SiteMap<(i64, i64)>, peak_window: i64) -> SiteMap<Vec<PyVal>> {
    let mut out: SiteMap<Vec<PyVal>> = IndexMap::new();
    for (chrm, mp) in m.iter() {
        let mut l_pos: Vec<i64> = mp.keys().cloned().collect();
        l_pos.sort();
        let mut res = Vec::new();
        cluster_positions(&l_pos, peak_window, |set| {
            let mut max_clip = 0;
            let mut max_all = 0;
            let mut cand = 0;
            for &p in set.iter() {
                let (c, a) = mp[&p];
                if max_clip < c || (max_clip == c && max_all < a) {
                    cand = p;
                    max_clip = c;
                    max_all = a;
                }
            }
            res.push((cand, max_clip));
        });
        let e = out.entry(chrm.clone()).or_default();
        for (cand, mc) in res {
            e.entry(cand).or_insert_with(|| vec![PyVal::Int(mc)]);
        }
    }
    out
}

/// `merge_clip_disc`
pub fn merge_clip_disc(sf_disc_tmp: &str, sf_clip: &str, sf_out: &str) -> Result<()> {
    let mut m_disc: IndexMap<String, std::collections::HashMap<String, (String, String)>> = IndexMap::new();
    for line in utils::read_lines(sf_disc_tmp)? {
        let f = split_ws(&line);
        m_disc.entry(f[0].to_string()).or_default().insert(f[1].to_string(), (f[2].to_string(), f[3].to_string()));
    }
    let mut fout = utils::create(sf_out)?;
    for line in utils::read_lines(sf_clip)? {
        let mut f: Vec<String> = split_ws(&line).iter().map(|s| s.to_string()).collect();
        let chrm = f[0].clone();
        let m = match m_disc.get(&chrm) {
            Some(m) => m,
            None => {
                println!("Error happen at merge clip and disc feature step: {} not exist", chrm);
                continue;
            }
        };
        let (l, r) = match m.get(&f[1]) {
            Some(x) => x.clone(),
            None => continue,
        };
        f.push(l);
        f.push(r);
        writeln!(fout, "{}", f.join("\t"))?;
    }
    fout.flush()?;
    Ok(())
}

/// `combine_closing_sites`
pub fn combine_closing_sites(sf_input: &str, iwindow: i64, sf_out: &str) -> Result<()> {
    let mut m: SiteMap<(i64, i64)> = IndexMap::new();
    let lines = utils::read_raw_lines(sf_input)?;
    for line in lines.iter() {
        let f = split_ws(line);
        let pos = utils::pint(f[1])?;
        let clip = utils::pint(f[2])? + utils::pint(f[3])? + utils::pint(f[4])?;
        let all = clip + utils::pint(f[5])? + utils::pint(f[6])?;
        m.entry(f[0].to_string()).or_default().insert(pos, (clip, all));
    }
    let peaks = call_peak_candidate_sites_all_features(&m, iwindow);
    let mut fout = utils::create(sf_out)?;
    for line in lines.iter() {
        let f = split_ws(line);
        let pos = utils::pint(f[1])?;
        if peaks.get(f[0]).map(|x| x.contains_key(&pos)).unwrap_or(false) {
            fout.write_all(line.as_bytes())?;
        }
    }
    fout.flush()?;
    Ok(())
}

/// `are_sites_close`
pub fn are_sites_close(pos1: i64, pos2: i64, iwin: i64) -> bool {
    pos2 >= pos1 - iwin && pos2 < pos1 + iwin
}
