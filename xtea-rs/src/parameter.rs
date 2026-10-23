//! Port of `x_parameter.py`: coverage dependent cutoffs.

use std::collections::BTreeMap;

/// A parameter record keyed by coverage.
#[derive(Clone, Debug)]
pub struct CovPars<T: Clone> {
    pub m_cov_par: BTreeMap<i64, T>,
}

impl<T: Clone> CovPars<T> {
    /// `PrototypeParameters.get_par_by_cov`
    pub fn get_par_by_cov(&self, icov: f64) -> T {
        let i_cur_cov = icov.trunc() as i64;
        for i in 0..10 {
            let i_minus = i_cur_cov - i;
            if let Some(v) = self.m_cov_par.get(&i_minus) {
                return v.clone();
            }
            let i_plus = i_cur_cov + i;
            if let Some(v) = self.m_cov_par.get(&i_plus) {
                return v.clone();
            }
        }
        // i_tmp = int(icov)/10 * 10 (float in Python 3)
        let i_tmp = (i_cur_cov as f64) / 10.0 * 10.0;
        if i_tmp < 10.0 {
            self.m_cov_par[&5].clone()
        } else if i_tmp > 200.0 {
            self.m_cov_par[&200].clone()
        } else {
            self.m_cov_par[&100].clone()
        }
    }
}

/// `Parameters` (germline): (nclip, ndisc, nclip_disc)
pub fn parameters() -> CovPars<(i64, i64, i64)> {
    let mut m = BTreeMap::new();
    for (k, v) in [
        (5, (1, 3, 0)),
        (10, (2, 3, 0)),
        (15, (2, 3, 0)),
        (20, (2, 4, 0)),
        (25, (2, 4, 1)),
        (30, (3, 4, 1)),
        (35, (3, 5, 1)),
        (40, (3, 6, 1)),
        (45, (4, 6, 1)),
        (50, (4, 7, 1)),
        (60, (5, 8, 2)),
        (70, (5, 9, 2)),
        (80, (6, 10, 2)),
        (90, (7, 11, 3)),
        (100, (8, 12, 3)),
        (110, (9, 12, 3)),
        (120, (9, 13, 3)),
        (130, (10, 14, 3)),
        (140, (11, 14, 4)),
        (150, (12, 15, 4)),
        (175, (13, 17, 5)),
        (200, (14, 20, 5)),
        (250, (20, 25, 7)),
        (300, (25, 30, 8)),
    ] {
        m.insert(k, v);
    }
    CovPars { m_cov_par: m }
}

/// `XParameters` (10X)
pub fn x_parameters() -> CovPars<(i64, i64, i64)> {
    let mut m = BTreeMap::new();
    for (k, v) in [
        (25, (4, 5, 1)),
        (30, (4, 6, 1)),
        (35, (4, 6, 1)),
        (40, (4, 7, 1)),
        (45, (5, 7, 1)),
        (50, (5, 7, 1)),
        (60, (5, 7, 2)),
        (70, (5, 7, 2)),
    ] {
        m.insert(k, v);
    }
    CovPars { m_cov_par: m }
}

/// `XParameters.get_barcode_cov_cutoff`
pub fn get_barcode_cov_cutoff(fcov: f64) -> i64 {
    let icov = fcov.trunc() as i64;
    if (icov - 30).abs() > (icov - 60).abs() {
        500
    } else {
        400
    }
}

/// `LRDParameters`
pub fn lrd_parameters() -> CovPars<i64> {
    let mut m = BTreeMap::new();
    for (k, v) in [
        (5, 1),
        (10, 3),
        (15, 4),
        (20, 4),
        (25, 5),
        (30, 6),
        (35, 6),
        (40, 7),
        (45, 8),
        (50, 9),
        (60, 10),
        (70, 11),
        (80, 12),
        (90, 13),
        (100, 14),
        (110, 15),
        (120, 16),
        (130, 17),
        (140, 18),
        (150, 19),
        (200, 20),
    ] {
        m.insert(k, v);
    }
    CovPars { m_cov_par: m }
}

/// `CaseControlFilterPars`
pub fn case_control_filter_pars() -> CovPars<(i64, i64, i64)> {
    let mut m = BTreeMap::new();
    for (k, v) in [
        (5, (1, 1, 0)),
        (10, (1, 2, 0)),
        (15, (2, 2, 0)),
        (20, (2, 3, 1)),
        (25, (2, 4, 1)),
        (30, (2, 4, 1)),
        (35, (2, 5, 1)),
        (40, (3, 5, 1)),
        (45, (3, 6, 1)),
        (50, (4, 6, 1)),
        (60, (5, 7, 2)),
        (70, (5, 8, 2)),
        (80, (6, 8, 2)),
        (90, (7, 9, 2)),
        (100, (8, 10, 3)),
        (200, (14, 20, 5)),
    ] {
        m.insert(k, v);
    }
    CovPars { m_cov_par: m }
}
