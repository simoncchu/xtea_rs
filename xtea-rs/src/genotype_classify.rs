//! Genotype classification of the called insertions (port of `x_genotype_classify.py`).
//!
//! The Python version predicts genotypes with a deep-forest (DF21, `deep-forest` 0.1.x)
//! `CascadeForestClassifier`.  The trained model (`genotyping/DF21_model_1_2`) is a
//! directory of joblib pickles; this module reads it directly and re-implements the
//! inference path of deep-forest (binning, forest prediction, cascade layers, probability
//! merging) so that predictions are identical to the Python package.

use crate::py::fmt::{float_repr, parse_py_float};
use crate::py::pickle::{self, Value};
use crate::utils;
use anyhow::{anyhow, bail, Context, Result};
use std::io::Write;
use std::path::Path;

/// Feature binning (`deepforest._binner.Binner`).
#[derive(Debug, Clone)]
pub struct Binner {
    pub thresholds: Vec<Vec<f64>>,
    pub missing_bin: u8,
    pub n_features: usize,
}

impl Binner {
    fn from_pickle(v: &Value) -> Result<Self> {
        let thresholds = v
            .req_attr("bin_thresholds_")?
            .as_list()?
            .iter()
            .map(|a| a.as_array()?.to_f64_vec())
            .collect::<Result<Vec<_>>>()?;
        let missing_bin = v.req_attr("missing_values_bin_idx_")?.as_i64()? as u8;
        let n_features = v.req_attr("n_bins_non_missing_")?.as_array()?.len();
        Ok(Binner { thresholds, missing_bin, n_features })
    }

    /// `Binner.transform`: index of the first threshold >= x (NaN -> missing bin).
    pub fn transform(&self, x: &[Vec<f64>]) -> Result<Vec<Vec<u8>>> {
        let mut out = Vec::with_capacity(x.len());
        for row in x {
            if row.len() != self.n_features {
                bail!(
                    "The binner was fitted with {} features but {} features got passed to `transform`.",
                    self.n_features,
                    row.len()
                );
            }
            let mut b = Vec::with_capacity(row.len());
            for (j, &v) in row.iter().enumerate() {
                if v.is_nan() {
                    b.push(self.missing_bin);
                    continue;
                }
                let t = &self.thresholds[j];
                let (mut left, mut right) = (0usize, t.len());
                while left < right {
                    let middle = (right + left - 1) / 2;
                    if v <= t[middle] {
                        right = middle;
                    } else {
                        left = middle + 1;
                    }
                }
                b.push(left as u8);
            }
            out.push(b);
        }
        Ok(out)
    }
}

/// One tree of a deep-forest forest: internal nodes plus leaf values.
#[derive(Debug, Clone)]
struct Tree {
    feature: Vec<usize>,
    threshold: Vec<u8>,
    /// (left, right); values <= 0 are leaves (leaf index = -value)
    children: Vec<(i64, i64)>,
    /// leaf values, `n_classes` per leaf
    value: Vec<f64>,
    n_leaf: usize,
}

impl Tree {
    /// `_forest._apply_region` for one sample.
    fn leaf(&self, x: &[u8]) -> usize {
        if self.feature.is_empty() {
            return 0;
        }
        let mut node = 0usize;
        loop {
            let c = if x[self.feature[node]] <= self.threshold[node] {
                self.children[node].0
            } else {
                self.children[node].1
            };
            if c <= 0 {
                return (-c) as usize;
            }
            node = c as usize;
        }
    }
}

/// `deepforest.forest.RandomForestClassifier` / `ExtraTreesClassifier` (inference only).
#[derive(Debug, Clone)]
struct Forest {
    trees: Vec<Tree>,
    n_classes: usize,
}

impl Forest {
    fn from_pickle(est: &Value, sf: &str) -> Result<Self> {
        if let Some(b) = est.getattr("backend") {
            if b.as_str().map(|s| s != "custom").unwrap_or(false) {
                bail!("{}: only deep-forest estimators with the 'custom' backend are supported", sf);
            }
        }
        if let Some(c) = est.getattr("is_classifier") {
            if !c.truthy() {
                bail!("{}: only classification models are supported", sf);
            }
        }
        let forest = est.req_attr("estimator_")?;
        let features = forest.req_attr("features")?.as_list()?;
        let thresholds = forest.req_attr("thresholds")?.as_list()?;
        let childrens = forest.req_attr("childrens")?.as_list()?;
        let values = forest.req_attr("values")?.as_list()?;
        let n_classes = forest.req_attr("n_classes_")?.as_i64()? as usize;
        if features.len() != thresholds.len() || features.len() != childrens.len() || features.len() != values.len() {
            bail!("{}: inconsistent forest arrays", sf);
        }
        let mut trees = Vec::with_capacity(features.len());
        for i in 0..features.len() {
            let f = features[i].as_array()?;
            let t = thresholds[i].as_array()?;
            let c = childrens[i].as_array()?;
            let v = values[i].as_array()?;
            let n_int = f.len();
            if t.len() != n_int || c.len() != 2 * n_int {
                bail!("{}: inconsistent tree arrays (tree {})", sf, i);
            }
            let feature = f.to_i64_vec()?.into_iter().map(|x| x as usize).collect();
            let threshold = t.to_i64_vec()?.into_iter().map(|x| x as u8).collect();
            let cv = c.to_i64_vec()?;
            let children = (0..n_int).map(|k| (cv[2 * k], cv[2 * k + 1])).collect();
            let value = v.to_f64_vec()?;
            let ncol = if v.shape.len() >= 2 { v.shape[1] } else { 1 };
            if ncol != n_classes {
                bail!("{}: leaf values have {} columns, expected {}", sf, ncol, n_classes);
            }
            let n_leaf = if v.shape.is_empty() { 0 } else { v.shape[0] };
            trees.push(Tree { feature, threshold, children, value, n_leaf });
        }
        Ok(Forest { trees, n_classes })
    }

    /// `ForestClassifier.predict_proba`: mean of the leaf distributions over the trees
    /// (accumulated in tree order, then divided by the number of trees).
    fn predict_proba(&self, xb: &[Vec<u8>]) -> Vec<Vec<f64>> {
        let mut out = vec![vec![0.0f64; self.n_classes]; xb.len()];
        for t in &self.trees {
            for (i, row) in xb.iter().enumerate() {
                let mut leaf = t.leaf(row);
                if t.n_leaf == 0 {
                    continue;
                }
                if leaf >= t.n_leaf {
                    leaf = t.n_leaf - 1; // numpy take(mode='clip')
                }
                let base = leaf * self.n_classes;
                for k in 0..self.n_classes {
                    out[i][k] += t.value[base + k];
                }
            }
        }
        let n = self.trees.len() as f64;
        for r in out.iter_mut() {
            for v in r.iter_mut() {
                *v /= n;
            }
        }
        out
    }
}

/// `deepforest.CascadeForestClassifier` loaded from a model directory.
pub struct Df21Model {
    pub n_layers: usize,
    pub n_features: usize,
    pub n_outputs: usize,
    pub n_estimators: usize,
    pub verbose: i64,
    classes: Option<Vec<i64>>,
    binners: Vec<Binner>,
    /// per layer: estimators in the order deep-forest loads them (`rf` first, then `erf`)
    layers: Vec<Vec<Forest>>,
}

/// `deepforest._utils.ctime()`: "[YYYY-mm-dd HH:MM:SS.mmm]" in UTC.
pub fn df_ctime() -> String {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default();
    let secs = now.as_secs() as i64;
    let ms = now.subsec_millis();
    let days = secs.div_euclid(86400);
    let rem = secs.rem_euclid(86400);
    // civil-from-days (H. Hinnant)
    let z = days + 719468;
    let era = z.div_euclid(146097);
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!(
        "[{:04}-{:02}-{:02} {:02}:{:02}:{:02}.{:03}]",
        y,
        m,
        d,
        rem / 3600,
        (rem % 3600) / 60,
        rem % 60,
        ms
    )
}

impl Df21Model {
    /// `CascadeForestClassifier.load(dirname)`.
    pub fn load(dirname: &str) -> Result<Self> {
        if !Path::new(dirname).is_dir() {
            bail!("Cannot find the target directory: {}.", dirname);
        }
        let join = |f: &str| Path::new(dirname).join(f).to_string_lossy().to_string();
        let d = pickle::load_file(&join("param.pkl"))?;
        let geti = |k: &str| -> Result<i64> {
            d.get_item(k).ok_or_else(|| anyhow!("param.pkl: missing '{}'", k))?.as_i64()
        };
        let n_estimators = geti("n_estimators")? as usize;
        let n_layers = geti("n_layers")? as usize;
        let n_features = geti("n_features")? as usize;
        let n_outputs = geti("n_outputs")? as usize;
        let verbose = d.get_item("verbose").map(|v| v.as_i64()).transpose()?.unwrap_or(1);
        let use_predictor = d.get_item("use_predictor").map(|v| v.truthy()).unwrap_or(false);
        let use_custom = d.get_item("use_custom_estimator").map(|v| v.truthy()).unwrap_or(false);
        let is_classifier = d.get_item("is_classifier").map(|v| v.truthy()).unwrap_or(true);
        if use_custom {
            bail!("{}: deep-forest models with custom estimators are not supported", dirname);
        }
        if use_predictor {
            bail!("{}: deep-forest models with a final predictor are not supported", dirname);
        }
        if !is_classifier {
            bail!("{}: deep-forest regression models are not supported", dirname);
        }
        let mut classes = None;
        if let Some(enc) = d.get_item("labels_are_encoded") {
            if enc.truthy() {
                let le = d.get_item("label_encoder").ok_or_else(|| anyhow!("param.pkl: missing 'label_encoder'"))?;
                let c = le.req_attr("classes_")?.as_array()?;
                if c.kind == 'f' {
                    classes = Some(c.to_f64_vec()?.into_iter().map(|x| x as i64).collect());
                } else {
                    classes = Some(c.to_i64_vec()?);
                }
            }
        }
        let b = pickle::load_file(&join("binner.pkl"))?;
        let mut binners = Vec::new();
        for i in 0..n_layers {
            let key = format!("binner_{}", i);
            let bv = b.get_item(&key).ok_or_else(|| anyhow!("binner.pkl: missing '{}'", key))?;
            binners.push(Binner::from_pickle(&bv).with_context(|| format!("binner.pkl: '{}'", key))?);
        }
        let mut layers = Vec::new();
        for layer_idx in 0..n_layers {
            let mut ests = Vec::new();
            for est_type in ["rf", "erf"] {
                for est_idx in 0..n_estimators {
                    let key = format!("{}-{}-{}", layer_idx, est_idx, est_type);
                    let dest = Path::new(dirname).join("estimator").join(format!("{}.est", key));
                    let dest_s = dest.to_string_lossy().to_string();
                    if !dest.is_file() {
                        bail!("Missing estimator in the path: {}.", dest_s);
                    }
                    let est = pickle::load_file(&dest_s)?;
                    ests.push(Forest::from_pickle(&est, &dest_s)?);
                }
            }
            for f in &ests {
                if f.n_classes != n_outputs {
                    bail!("{}: estimator outputs ({}) do not match n_outputs ({})", dirname, f.n_classes, n_outputs);
                }
            }
            layers.push(ests);
        }
        Ok(Df21Model { n_layers, n_features, n_outputs, n_estimators, verbose, classes, binners, layers })
    }

    fn n_aug_features(&self) -> usize {
        2 * self.n_estimators * self.n_outputs
    }

    /// `ClassificationCascadeLayer.predict_full`: concatenated class distributions.
    fn layer_predict_full(&self, layer: usize, xb: &[Vec<u8>]) -> Vec<Vec<f64>> {
        let ests = &self.layers[layer];
        let width = self.n_outputs * ests.len();
        let mut pred = vec![vec![0.0f64; width]; xb.len()];
        for (idx, est) in ests.iter().enumerate() {
            let p = est.predict_proba(xb);
            let left = self.n_outputs * idx;
            for (i, r) in p.iter().enumerate() {
                for k in 0..self.n_outputs {
                    pred[i][left + k] += r[k];
                }
            }
        }
        pred
    }

    /// `_utils.merge_proba`: average of the per-estimator distributions.
    fn merge_proba(&self, probas: &[Vec<f64>]) -> Vec<Vec<f64>> {
        probas
            .iter()
            .map(|row| {
                let mut out = vec![0.0f64; self.n_outputs];
                let mut count = 0usize;
                let mut start = 0usize;
                while start < row.len() {
                    for k in 0..self.n_outputs {
                        out[k] += row[start + k];
                    }
                    start += self.n_outputs;
                    count += 1;
                }
                for v in out.iter_mut() {
                    *v /= count as f64;
                }
                out
            })
            .collect()
    }

    /// `CascadeForestClassifier.predict_proba`.
    pub fn predict_proba(&self, x: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
        if self.verbose > 0 {
            println!("{} Start to evalute the model:", df_ctime());
        }
        let x_test = self.binners[0].transform(x)?;
        let n_aug = self.n_aug_features();
        // X_middle: binned input followed by the binned augmented features
        let mut x_middle: Vec<Vec<u8>> = x_test
            .iter()
            .map(|r| {
                let mut v = r.clone();
                v.resize(r.len() + n_aug, 0);
                v
            })
            .collect();
        let mut x_aug: Vec<Vec<f64>> = Vec::new();
        for layer_idx in 0..self.n_layers {
            if self.verbose > 0 {
                println!("{} Evaluating cascade layer = {:<2}", df_ctime(), layer_idx);
            }
            if layer_idx == 0 {
                x_aug = self.layer_predict_full(0, &x_test);
            } else {
                let xb = self.binners[layer_idx].transform(&x_aug)?;
                for (row, b) in x_middle.iter_mut().zip(xb.iter()) {
                    if self.n_features + b.len() != row.len() {
                        bail!("unexpected layout of the cascade layer input");
                    }
                    row[self.n_features..].copy_from_slice(b);
                }
                if layer_idx < self.n_layers - 1 {
                    x_aug = self.layer_predict_full(layer_idx, &x_middle);
                }
            }
        }
        let proba = if self.n_layers > 1 {
            let p = self.layer_predict_full(self.n_layers - 1, &x_middle);
            self.merge_proba(&p)
        } else {
            self.merge_proba(&x_aug)
        };
        Ok(proba)
    }

    /// `CascadeForestClassifier.predict`: argmax (first maximum) decoded to the labels.
    pub fn predict(&self, x: &[Vec<f64>]) -> Result<Vec<i64>> {
        let proba = self.predict_proba(x)?;
        Ok(proba
            .iter()
            .map(|r| {
                let mut best = 0usize;
                for (k, &v) in r.iter().enumerate() {
                    if v > r[best] || (v.is_nan() && !r[best].is_nan()) {
                        best = k;
                        if v.is_nan() {
                            break;
                        }
                    }
                }
                match &self.classes {
                    Some(c) => c.get(best).copied().unwrap_or(best as i64),
                    None => best as i64,
                }
            })
            .collect())
    }
}

/// `GntpClassifier_DF21`.
pub struct GntpClassifierDf21 {
    pub n_feature: usize,
}

impl Default for GntpClassifierDf21 {
    fn default() -> Self {
        Self::new()
    }
}

/// A feature value kept both as the text written to the ARFF file and as its float value.
fn feat_f(v: f64) -> String {
    float_repr(v)
}

fn pf(s: &str) -> Result<f64> {
    parse_py_float(s).ok_or_else(|| anyhow!("could not convert string to float: '{}'", s))
}

impl GntpClassifierDf21 {
    pub fn new() -> Self {
        GntpClassifierDf21 { n_feature: 15 }
    }

    const ARFF_HEADER: &'static str = "@RELATION\tinsgntp\n\n\
@ATTRIBUTE\tlclipcns\tNUMERIC\n\
@ATTRIBUTE\trclipcns\tNUMERIC\n\
@ATTRIBUTE\tldisccns\tNUMERIC\n\
@ATTRIBUTE\trdisccns\tNUMERIC\n\
@ATTRIBUTE\tpolyA\tNUMERIC\n\
@ATTRIBUTE\tlcov\tNUMERIC\n\
@ATTRIBUTE\trcov\tNUMERIC\n\
@ATTRIBUTE\tclip\tNUMERIC\n\
@ATTRIBUTE\tfullmap\tNUMERIC\n\
@ATTRIBUTE\tclipratio\tNUMERIC\n\
@ATTRIBUTE\tdiscratio\tNUMERIC\n\
@ATTRIBUTE\trawlclip\tNUMERIC\n\
@ATTRIBUTE\trawrclip\tNUMERIC\n\
@ATTRIBUTE\tdiscordant\tNUMERIC\n\
@ATTRIBUTE\tconcordant\tNUMERIC\n\
@ATTRIBUTE\tclass\t{1,2}\n\n\
@DATA\n";

    /// `_parser_features`: the 15 features (as the strings Python writes) of one record.
    pub fn parser_features(&mut self, l_fields: &[&str]) -> Result<Vec<String>> {
        let g = |i: usize| -> Result<&str> {
            l_fields.get(i).copied().ok_or_else(|| anyhow!("list index out of range"))
        };
        let mut l = Vec::with_capacity(15);
        let mut f_lcov = pf(g(11)?)?;
        if f_lcov == 0.0 {
            f_lcov = 0.0000000001;
        }
        let mut f_rcov = pf(g(12)?)?;
        if f_rcov == 0.0 {
            f_rcov = 0.0000000001;
        }
        l.push(feat_f(pf(g(5)?)? / f_lcov));
        l.push(feat_f(pf(g(6)?)? / f_rcov));
        l.push(feat_f(pf(g(7)?)? / f_lcov));
        l.push(feat_f(pf(g(8)?)? / f_rcov));
        l.push(feat_f(pf(g(9)?)? / f_lcov + pf(g(10)?)? / f_rcov));
        l.push(g(11)?.to_string());
        l.push(g(12)?.to_string());
        l.push(feat_f(pf(g(35)?)? / (f_lcov + f_rcov)));
        l.push(feat_f(pf(g(36)?)? / (f_lcov + f_rcov)));
        let (f35, f36) = (pf(g(35)?)?, pf(g(36)?)?);
        if f36 + f35 > 0.0 {
            l.push(feat_f(f35 / (f36 + f35)));
        } else {
            l.push("0".to_string());
        }
        let (f39, f40) = (pf(g(39)?)?, pf(g(40)?)?);
        if f39 + f40 > 0.0 {
            l.push(feat_f(f39 / (f39 + f40)));
        } else {
            l.push("0".to_string());
        }
        l.push(feat_f(pf(g(37)?)? / f_lcov));
        l.push(feat_f(pf(g(38)?)? / f_rcov));
        l.push(feat_f(f39 / (f_lcov + f_rcov)));
        l.push(feat_f(f40 / (f_lcov + f_rcov)));
        self.n_feature = l.len();
        Ok(l)
    }

    /// `load_in_feature_from_xTEA_output` (prediction mode: class "1" appended).
    pub fn load_in_feature_from_xtea_output(&mut self, sf_xtea: &str) -> Result<Vec<Vec<String>>> {
        let mut all = Vec::new();
        if !utils::is_file(sf_xtea) {
            return Ok(all);
        }
        for line in utils::read_raw_lines(sf_xtea)? {
            let fields = utils::split_ws(&line);
            let mut f = self.parser_features(&fields)?;
            f.push("1".to_string());
            all.push(f);
        }
        Ok(all)
    }

    /// `prepare_arff_from_xTEA_output`: write the ARFF file and return the feature matrix as
    /// scipy's ARFF reader would parse it.
    pub fn prepare_arff_from_xtea_output(&mut self, sf_xtea: &str, sf_arff: &str) -> Result<Vec<Vec<f64>>> {
        let feats = self.load_in_feature_from_xtea_output(sf_xtea)?;
        {
            let mut w = utils::create(sf_arff)?;
            w.write_all(Self::ARFF_HEADER.as_bytes())?;
            for rcd in &feats {
                writeln!(w, "{}", rcd.join(","))?;
            }
            w.flush()?;
        }
        let n = self.n_feature;
        let mut x = Vec::with_capacity(feats.len());
        for rcd in &feats {
            let mut row = Vec::with_capacity(n);
            for s in rcd.iter().take(n) {
                // scipy.io.arff NumericAttribute.parse_data
                if s.contains('?') {
                    row.push(f64::NAN);
                } else {
                    row.push(pf(s)?);
                }
            }
            x.push(row);
        }
        Ok(x)
    }

    /// `predict_for_site`: append the predicted genotype to every record.
    pub fn predict_for_site(&mut self, sf_model: &str, sf_xtea: &str, sf_new: &str) -> Result<()> {
        let model = Df21Model::load(sf_model)?;
        let sf_arff = format!("{}.arff", sf_xtea);
        let site_features = self.prepare_arff_from_xtea_output(sf_xtea, &sf_arff)?;
        let preds = if !site_features.is_empty() { Some(model.predict(&site_features)?) } else { None };
        let lines = utils::read_raw_lines(sf_xtea)?;
        let mut w = utils::create(sf_new)?;
        let preds = match preds {
            None => return Ok(()),
            Some(p) => p,
        };
        for (i_idx, line) in lines.iter().enumerate() {
            let sinfo = utils::rstrip(line);
            let p = *preds.get(i_idx).ok_or_else(|| anyhow!("index out of bounds"))?;
            let s_gntp = if p == 1 {
                "0/1"
            } else if p == 2 {
                "1/1"
            } else {
                "0/0"
            };
            writeln!(w, "{}\t{}", sinfo, s_gntp)?;
        }
        w.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ctime_format() {
        let s = df_ctime();
        assert_eq!(s.len(), 25);
        assert!(s.starts_with('[') && s.ends_with(']'));
    }

    #[test]
    fn binner_searchsorted_left() {
        let b = Binner { thresholds: vec![vec![1.0, 2.0, 3.0]], missing_bin: 255, n_features: 1 };
        let x = vec![vec![0.5], vec![1.0], vec![1.5], vec![3.0], vec![9.0], vec![f64::NAN]];
        let r: Vec<u8> = b.transform(&x).unwrap().into_iter().map(|v| v[0]).collect();
        assert_eq!(r, vec![0, 0, 1, 2, 3, 255]);
    }

    #[test]
    fn tree_traversal_matches_deepforest() {
        // same cases as checked against deepforest._forest.predict
        let t = Tree {
            feature: vec![0, 1],
            threshold: vec![9, 100],
            children: vec![(1, -2), (0, -1)],
            value: vec![1.0, 0.0, 0.0, 1.0, 0.5, 0.5],
            n_leaf: 3,
        };
        assert_eq!(t.leaf(&[0, 5]), 0);
        assert_eq!(t.leaf(&[10, 200]), 2);
        assert_eq!(t.leaf(&[255, 3]), 2);
        let empty = Tree { feature: vec![], threshold: vec![], children: vec![], value: vec![0.25, 0.75], n_leaf: 1 };
        assert_eq!(empty.leaf(&[1, 2]), 0);
    }

    /// Compare with predictions made by the Python package (set XTEA_DF21_MODEL to the model
    /// directory and XTEA_DF21_CASES to a TSV of 15 features followed by the 2 probabilities).
    #[test]
    #[ignore]
    fn df21_matches_python() {
        let (Ok(m), Ok(c)) = (std::env::var("XTEA_DF21_MODEL"), std::env::var("XTEA_DF21_CASES")) else {
            return;
        };
        let mut model = Df21Model::load(&m).unwrap();
        model.verbose = 0;
        let mut x = Vec::new();
        let mut expect = Vec::new();
        for line in std::fs::read_to_string(&c).unwrap().lines() {
            let v: Vec<f64> = line.split('\t').map(|s| parse_py_float(s).unwrap()).collect();
            x.push(v[..15].to_vec());
            expect.push(v[15..].to_vec());
        }
        let p = model.predict_proba(&x).unwrap();
        let mut n_diff = 0;
        for (a, b) in p.iter().zip(expect.iter()) {
            if a[0].to_bits() != b[0].to_bits() || a[1].to_bits() != b[1].to_bits() {
                n_diff += 1;
            }
        }
        assert_eq!(n_diff, 0, "{} of {} predictions differ", n_diff, p.len());
    }
}
