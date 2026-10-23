//! Ports of `x_cluster_consistency.py` (ClusterChecker) and `disc_cluster.py` (DiscCluster).

use crate::py::set::PySet;
use indexmap::IndexMap;

/// `ClusterChecker.find_max_cover` (positions must be sorted).
pub fn find_max_cover(l_pos: &[i64], win_size: i64) -> (usize, usize) {
    let mut i_right: usize = 0;
    let i_total = l_pos.len();
    let mut i_cnt: i64 = 0;
    let mut i_max: i64 = 0;
    let mut i_start: usize = 0;
    let mut i_end: usize = 0;
    for i_left in 0..i_total {
        let cur_val = l_pos[i_left];
        if i_right < i_left {
            i_right = i_left;
        }
        let mut broke = false;
        for i_temp in i_right..i_total {
            let right_val = l_pos[i_temp];
            i_right = i_temp;
            if right_val - cur_val > win_size {
                if i_max < i_cnt {
                    i_max = i_cnt;
                    i_start = i_left;
                    i_end = i_right - 1;
                }
                broke = true;
                break;
            }
            i_cnt += 1;
        }
        let _ = broke;
        if i_max < i_cnt {
            i_max = i_cnt;
            i_start = i_left;
            i_end = i_right;
        }
        i_cnt -= 1;
        if i_cnt < 0 {
            i_cnt = 0;
        }
        if i_right == i_total - 1 {
            break;
        }
    }
    (i_start, i_end)
}

/// `ClusterChecker._is_disc_cluster` -> (is_cluster, left_boundary, right_boundary).
/// The input list is sorted in place like the original.
pub fn is_disc_cluster(l_pos: &mut Vec<i64>, i_dist: i64, ratio: f64) -> (bool, i64, i64) {
    l_pos.sort();
    let n_all = l_pos.len();
    if n_all == 0 {
        return (false, -1, -1);
    }
    let (i_max_left, i_max_right) = find_max_cover(l_pos, i_dist);
    let peak_win_cnt = i_max_right as i64 - i_max_left as i64;
    let l_boundary = l_pos[i_max_left];
    let r_boundary = l_pos[i_max_right];
    if (peak_win_cnt as f64) / (n_all as f64) < ratio {
        return (false, -1, -1);
    }
    (true, l_boundary, r_boundary)
}

/// `ClusterChecker._is_disc_orientation_consistency`
pub fn is_disc_orientation_consistency(n_l_rc: i64, n_l_non_rc: i64, n_r_rc: i64, n_r_non_rc: i64, ratio: f64) -> bool {
    let mut b_l = true;
    if n_l_rc + n_l_non_rc > 0 && (n_l_non_rc as f64) / ((n_l_rc + n_l_non_rc) as f64) < ratio {
        b_l = false;
    }
    let mut b_r = true;
    if n_r_rc + n_r_non_rc > 0 && (n_r_rc as f64) / ((n_r_rc + n_r_non_rc) as f64) < ratio {
        b_r = false;
    }
    b_l && b_r
}

/// `ClusterChecker.find_first_second_peak`.  `mpos` maps ref position -> number of reads.
/// Returns (peak_pos, peak_nbrs, peak_acm, second_pos, second_nbrs, second_acm).
pub fn find_first_second_peak(mpos: &IndexMap<i64, usize>, bin_size: i64) -> (i64, Vec<i64>, i64, i64, Vec<i64>, i64) {
    let mut peak_pos: i64 = -1;
    let mut peak_acm: i64 = 0;
    let mut scnd_peak_pos: i64 = -1;
    let mut scnd_peak_acm: i64 = 0;
    let istart = -(bin_size / 2);
    let iend = bin_size / 2;
    for (&pos, _) in mpos.iter() {
        if pos <= 0 {
            continue;
        }
        let mut cur_acm: i64 = 0;
        for offset in istart..iend {
            if let Some(n) = mpos.get(&(pos + offset)) {
                cur_acm += *n as i64;
            }
        }
        if cur_acm > peak_acm {
            scnd_peak_pos = peak_pos;
            scnd_peak_acm = peak_acm;
            peak_acm = cur_acm;
            peak_pos = pos;
        } else if cur_acm > scnd_peak_acm {
            scnd_peak_acm = cur_acm;
            scnd_peak_pos = pos;
        }
    }
    let mut l_peak_nbr = Vec::new();
    let mut l_scnd_peak_nbr = Vec::new();
    for offset in istart..iend {
        if mpos.contains_key(&(peak_pos + offset)) {
            l_peak_nbr.push(peak_pos + offset);
        }
        if mpos.contains_key(&(scnd_peak_pos + offset)) {
            l_scnd_peak_nbr.push(scnd_peak_pos + offset);
        }
    }
    (peak_pos, l_peak_nbr, peak_acm, scnd_peak_pos, l_scnd_peak_nbr, scnd_peak_acm)
}

/// `DiscCluster.form_one_side_cluster`: m_pos maps mate chromosome -> positions (sorted in place).
pub fn form_one_side_cluster(m_pos: &mut IndexMap<String, Vec<i64>>, i_is: i64, f_ratio: f64) -> (bool, Option<String>, Option<i64>) {
    let mut max_cnt: usize = 0;
    let mut max_pos: i64 = -1;
    let mut max_chrm = String::new();
    let mut all_cnt: usize = 0;
    if m_pos.is_empty() {
        return (false, None, None);
    }
    for (chrm, l_pos) in m_pos.iter_mut() {
        l_pos.sort();
        let mut pre_pos: i64 = -1;
        let mut set_cluster: Vec<i64> = Vec::new();
        for &pos in l_pos.iter() {
            if pre_pos == -1 {
                pre_pos = pos;
                set_cluster.push(pre_pos);
                continue;
            }
            if pos - pre_pos > i_is {
                let set_size = set_cluster.len();
                let mid_pos = set_cluster[set_size / 2];
                if max_cnt < set_size {
                    max_cnt = set_size;
                    max_pos = mid_pos;
                    max_chrm = chrm.clone();
                }
                all_cnt += set_size;
                set_cluster.clear();
            }
            pre_pos = pos;
            set_cluster.push(pre_pos);
        }
        let set_size = set_cluster.len();
        if set_size == 0 {
            // the original would raise IndexError here
            continue;
        }
        let mid_pos = set_cluster[set_size / 2];
        if max_cnt < set_size {
            max_cnt = set_size;
            max_pos = mid_pos;
            max_chrm = chrm.clone();
        }
        all_cnt += set_size;
    }
    if all_cnt > 0 && (max_cnt as f64) / (all_cnt as f64) > f_ratio {
        return (true, Some(max_chrm), Some(max_pos));
    }
    (false, None, None)
}

/// Helper: iterate a Python-ordered int set built from the given values.
pub fn pyset_from(values: &[i64]) -> PySet<i64> {
    PySet::from_iter_ordered(values.iter().cloned())
}
