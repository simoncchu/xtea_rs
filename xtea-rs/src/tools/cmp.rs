//! Port of `x_cmp.py` (stand-alone tool with its own main): compare TE insertion calls with a
//! benchmark set.
//!
//! * single sample (default): `-i calls --bcmk benchmark -o out [--bed regions] [--exclude
//!   sites] [--extnd N]` writes `out` (benchmark sites that were missed), `out.true_positive`,
//!   `out.false_positive` and the two `._with_info` files;
//! * group (`-L --list list`, each line `result_file repeat_type sample`): per-result hit
//!   counts, `out_ori.csv`, `out_fall_in_hc.csv`, `out_hit.csv`, `out.f1_score`,
//!   `out.tp_fp_barchart(.2)`.
//!
//! Quirks of the original that are kept on purpose:
//! * matches are searched in `[pos - extnd, pos + extnd)` and the smallest matching position is
//!   used; chromosome names of the benchmark are given/stripped the "chr" prefix depending on
//!   whether the benchmark (second dict) contains "chr1" (`process_chrm_name`);
//! * `--bed` runs `bedtools intersect -a <bed> -b <calls>`, i.e. the *regions* overlapping the
//!   calls are compared, with the regions' extra columns;
//! * group mode overwrites `out` / `out.true_positive` for every result and writes
//!   `out<sample>..true_positive`; TP/FP/FN use line counts of the raw files (header lines
//!   included); `s_sample.split("_")[1]` must exist;
//! * `get_joint_info` compares positions as strings and `get_joint_info2` parses every line of
//!   the benchmark (header lines raise, like in the script).

use super::{at, create_file, need, open_text, py_fdiv, py_int, py_repr_str, py_split_ws, read_text_lines};
use crate::py::fmt::float_repr;
use crate::utils::{is_file, run_cmd_print};
use anyhow::{anyhow, Result};
use clap::Parser;
use indexmap::{IndexMap, IndexSet};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::io::Write;

/// `_fp_suffix`
pub const FP_SUFFIX: &str = ".false_positive";
/// `_tp_suffix`
pub const TP_SUFFIX: &str = ".true_positive";

/// Command line of `x_cmp.py` (optparse options, same names and defaults).
#[derive(Parser, Debug, Clone, Default)]
#[command(
    name = "xtea_cmp",
    about = "Compare TE insertion calls against a benchmark (Rust port of x_cmp.py)",
    infer_long_args = true,
    args_override_self = true,
    disable_help_flag = true
)]
pub struct Opts {
    #[arg(short = 'h', long = "help", help = "show this help message and exit")]
    pub help: bool,
    #[arg(short = 'L', help = "Input is a list of results")]
    pub bgroup: bool,
    #[arg(allow_hyphen_values = true, short = 'i', long = "input", value_name = "FILE", help = "Input file")]
    pub input: Option<String>,
    #[arg(allow_hyphen_values = true, long = "list", value_name = "FILE", help = "List of results")]
    pub list: Option<String>,
    #[arg(allow_hyphen_values = true, long = "bcmk", value_name = "FILE", help = "Benchmark file")]
    pub bcmk: Option<String>,
    #[arg(allow_hyphen_values = true, long = "bed", value_name = "FILE", help = "Compare within the given region")]
    pub focal_region: Option<String>,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, long = "extnd", default_value_t = 200, help = "Slack value")]
    pub extnd: i64,
    #[arg(allow_hyphen_values = true, long = "exclude", value_name = "FILE", help = "Sites need to be excluded")]
    pub exclude: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'o', long = "output", value_name = "FILE", help = "The output file")]
    pub output: Option<String>,
    /// Positional arguments are accepted and ignored (optparse collects them separately).
    #[arg(hide = true)]
    pub extra_args: Vec<String>,
}

/// The script's `__main__`.
pub fn run(opts: Opts) -> Result<()> {
    let sf_bcmk = opts.bcmk.as_deref();
    let sf_hc_bed = opts.focal_region.as_deref();
    let i_extd = opts.extnd;
    let sf_unhit = opts.output.as_deref();
    let sf_exclude = opts.exclude.as_deref();
    if opts.bgroup {
        let (m_cnt, l_sample_info, m_tp_fp_fn) = cmp_for_list_of_rslts(opts.list.as_deref(), sf_bcmk, sf_exclude, sf_hc_bed, i_extd, sf_unhit)?;
        let sf_csv = need(sf_unhit, "-o")?;
        export_cnt_to_csv_for_point_chart(&l_sample_info, &m_cnt, sf_csv)?;
        calc_export_f1_score_sensitivity_fdr_csv(&m_tp_fp_fn, &format!("{}.f1_score", sf_csv))?;
        calc_export_f1_score_sensitivity_fdr_csv2(&m_tp_fp_fn, &format!("{}.tp_fp_barchart", sf_csv))?;
    } else {
        let s_type = "None";
        let sf_input = need(opts.input.as_deref(), "-i")?;
        let sf_new_input = format!("{}.tmp", sf_input);
        gnrt_bed_from_input(sf_input, &sf_new_input)?;
        check_one_sample_agnst_bcmk2(&sf_new_input, sf_bcmk, sf_exclude, sf_hc_bed, i_extd, s_type, sf_unhit)?;
        let sf_unhit = need(sf_unhit, "-o")?;

        // the false positive list (with info)
        let sf_slct_sites = format!("{}{}", sf_unhit, FP_SUFFIX);
        let sf_joint = format!("{}._with_info", sf_slct_sites);
        get_joint_info(sf_input, &sf_slct_sites, &sf_joint)?;

        // the true positive list (with info), from the benchmark
        let sf_tp_sites = format!("{}{}", sf_unhit, TP_SUFFIX);
        let sf_joint_tp = format!("{}._with_info", sf_tp_sites);
        get_joint_info2(need(sf_bcmk, "--bcmk")?, &sf_tp_sites, i_extd, &sf_joint_tp)?;
    }
    Ok(())
}

/// Sites of one chromosome: `{pos: [values...]}` in insertion order plus the sorted key set
/// used by the window look-ups.
#[derive(Debug, Default, Clone)]
pub struct ChrSites {
    pub map: IndexMap<i64, Vec<String>>,
    keys: BTreeSet<i64>,
}

impl ChrSites {
    fn entry(&mut self, pos: i64) -> &mut Vec<String> {
        self.keys.insert(pos);
        self.map.entry(pos).or_default()
    }

    fn set(&mut self, pos: i64, values: Vec<String>) {
        self.keys.insert(pos);
        self.map.insert(pos, values);
    }
}

/// `{chrm: {pos: [values...]}}` as returned by `load_sites`.
pub type Sites = IndexMap<String, ChrSites>;
/// `{chrm: {pos: 1}}`
pub type PosDict = IndexMap<String, IndexSet<i64>>;

/// `RsltCMP.load_sites`: "#" lines are skipped (unless the line is just "#"), the remaining
/// columns are appended to the values of the position.
pub fn load_sites(sf_input: &str) -> Result<Sites> {
    let mut m_sites = Sites::new();
    for line in read_text_lines(sf_input)? {
        if line.len() > 1 && line.starts_with('#') {
            continue;
        }
        let fields = py_split_ws(&line);
        let chrm = at(&fields, 0)?;
        let pos = py_int(at(&fields, 1)?)?;
        let values = m_sites.entry(chrm.to_string()).or_default().entry(pos);
        values.extend(fields[2..].iter().map(|s| s.to_string()));
    }
    Ok(m_sites)
}

/// `RsltCMP.process_chrm_name`
pub fn process_chrm_name(chrm: &str, b_with_chr: bool) -> String {
    let b_chrm_with_chr = chrm.len() > 3 && chrm.starts_with("chr");
    match (b_with_chr, b_chrm_with_chr) {
        (true, true) => chrm.to_string(),
        (true, false) => format!("chr{}", chrm),
        (false, true) => chrm[3..].to_string(),
        (false, false) => chrm.to_string(),
    }
}

/// Result of `cmp_dict`: (n_overlap, m_unhit, m_hit, m_hit_coordinate_bcmk).
pub struct CmpResult {
    pub n_overlap: i64,
    pub m_unhit: PosDict,
    pub m_hit: PosDict,
    pub m_hit_coordinate_bcmk: PosDict,
}

/// `RsltCMP.cmp_dict`: for every benchmark site look for a site of `m_dict1` within
/// `[pos - extnd, pos + extnd)` (the global EXTND of the script).
pub fn cmp_dict(m_dict1: &Sites, m_bcmk: &Sites, extnd: i64) -> CmpResult {
    let mut n_overlap = 0;
    let mut m_unhit = PosDict::new();
    let mut m_hit = PosDict::new();
    let mut m_hit_coordinate_bcmk = PosDict::new();
    let b_contain_chr = m_bcmk.contains_key("chr1");
    for (chrm_1, bcmk_sites) in m_bcmk {
        let chrm = process_chrm_name(chrm_1, b_contain_chr);
        let dict1_sites = match m_dict1.get(&chrm) {
            Some(s) => s,
            None => {
                let unhit = m_unhit.entry(chrm).or_default();
                for &tmp_pos in bcmk_sites.map.keys() {
                    unhit.insert(tmp_pos);
                }
                continue;
            }
        };
        for &pos in bcmk_sites.map.keys() {
            match super::first_key_in_window(&dict1_sites.keys, pos, extnd) {
                Some(hit_pos) => {
                    n_overlap += 1;
                    if !m_hit.contains_key(&chrm) {
                        m_hit.insert(chrm.clone(), IndexSet::new());
                        m_hit_coordinate_bcmk.insert(chrm.clone(), IndexSet::new());
                    }
                    m_hit[&chrm].insert(hit_pos);
                    m_hit_coordinate_bcmk[&chrm].insert(pos);
                }
                None => {
                    m_unhit.entry(chrm.clone()).or_default().insert(pos);
                }
            }
        }
    }
    CmpResult { n_overlap, m_unhit, m_hit, m_hit_coordinate_bcmk }
}

/// `RsltCMP.exclude_sites`: drop the sites with an excluded site of the same chromosome name
/// within `[pos - i_extnd, pos + i_extnd)`.
pub fn exclude_sites(m_input: &Sites, m_exclude: &Sites, i_extnd: i64) -> Sites {
    let mut m_new_input = Sites::new();
    for (chrm, chr_sites) in m_input {
        let excluded = m_exclude.get(chrm);
        for (&pos, values) in &chr_sites.map {
            if let Some(ex) = excluded {
                if super::first_key_in_window(&ex.keys, pos, i_extnd).is_some() {
                    continue;
                }
            }
            m_new_input.entry(chrm.clone()).or_default().set(pos, values.clone());
        }
    }
    m_new_input
}

/// Look up `m[chrm][pos]` like the script (KeyError -> Err).
fn site_values<'a>(m: &'a Sites, chrm: &str, pos: i64) -> Result<&'a Vec<String>> {
    m.get(chrm)
        .ok_or_else(|| anyhow!("KeyError: {}", py_repr_str(chrm)))?
        .map
        .get(&pos)
        .ok_or_else(|| anyhow!("KeyError: {}", pos))
}

fn write_site_line<W: Write>(w: &mut W, chrm: &str, pos: i64, values: &[String]) -> Result<()> {
    w.write_all(format!("{}\t{}\t{}\n", chrm, pos, values.join("\t")).as_bytes())?;
    Ok(())
}

fn write_unhit(sf_unhit: &str, m_unhit: &PosDict) -> Result<()> {
    let mut fout_unhit = create_file(sf_unhit)?;
    for (chrm, positions) in m_unhit {
        for pos in positions {
            fout_unhit.write_all(format!("{}\t{}\n", chrm, pos).as_bytes())?;
        }
    }
    fout_unhit.flush()?;
    Ok(())
}

fn write_hits(sf_out: &str, m_hit: &PosDict, m_sites: &Sites) -> Result<()> {
    let mut fout = create_file(sf_out)?;
    for (chrm, positions) in m_hit {
        for &pos in positions {
            let values = site_values(m_sites, chrm, pos)?;
            write_site_line(&mut fout, chrm, pos, values)?;
        }
    }
    fout.flush()?;
    Ok(())
}

/// `RsltCMP.cmp_rslts` (group mode)
pub fn cmp_rslts(sf_input: &str, sf_bcmk: &str, i_extd: i64, sf_unhit: Option<&str>, s_sample: &str, m_exclude: &Sites) -> Result<i64> {
    let m_bcmk1 = load_sites(sf_bcmk)?;
    let m_input1 = load_sites(sf_input)?;
    let m_bcmk = exclude_sites(&m_bcmk1, m_exclude, i_extd);
    let m_input = exclude_sites(&m_input1, m_exclude, i_extd);

    let n_bcmk = cmp_dict(&m_bcmk, &m_bcmk, i_extd).n_overlap;
    let r2 = cmp_dict(&m_input, &m_bcmk, i_extd);
    println!("{}/{} hits/all", r2.n_overlap, n_bcmk);
    let sf_unhit = need(sf_unhit, "-o")?;
    write_unhit(sf_unhit, &r2.m_unhit)?;
    write_hits(&format!("{}{}", sf_unhit, TP_SUFFIX), &r2.m_hit, &m_input)?;
    write_hits(&format!("{}{}.{}", sf_unhit, s_sample, TP_SUFFIX), &r2.m_hit_coordinate_bcmk, &m_bcmk)?;
    Ok(r2.n_overlap)
}

/// `RsltCMP.cmp_rslts2` (single-sample mode, also writes the false positives)
pub fn cmp_rslts2(sf_input: &str, sf_bcmk: Option<&str>, i_extd: i64, sf_unhit: Option<&str>, m_exclude: &Sites) -> Result<i64> {
    let m_bcmk1 = load_sites(need(sf_bcmk, "--bcmk")?)?;
    let m_input1 = load_sites(sf_input)?;
    let m_bcmk = exclude_sites(&m_bcmk1, m_exclude, i_extd);
    let m_input = exclude_sites(&m_input1, m_exclude, i_extd);

    let n_bcmk = cmp_dict(&m_bcmk, &m_bcmk, i_extd).n_overlap;
    let r2 = cmp_dict(&m_input, &m_bcmk, i_extd);
    let r3 = cmp_dict(&m_bcmk, &m_input, i_extd);
    println!("{}/{} hits/all", r2.n_overlap, n_bcmk);
    let sf_unhit = need(sf_unhit, "-o")?;
    write_unhit(sf_unhit, &r2.m_unhit)?;
    write_hits(&format!("{}{}", sf_unhit, TP_SUFFIX), &r2.m_hit, &m_input)?;
    {
        let mut fout_fp = create_file(format!("{}{}", sf_unhit, FP_SUFFIX))?;
        for (chrm, positions) in &r3.m_unhit {
            for &pos in positions {
                let values = match m_input.get(chrm).and_then(|s| s.map.get(&pos)) {
                    Some(v) => v,
                    None => continue,
                };
                write_site_line(&mut fout_fp, chrm, pos, values)?;
            }
        }
        fout_fp.flush()?;
    }
    Ok(r2.n_overlap)
}

/// `RsltCMP.gnrt_ins_within_focal_region`
pub fn gnrt_ins_within_focal_region(sf_var: &str, sf_hc_bed: Option<&str>, sf_new_var: &str) {
    if let Some(sf_hc_bed) = sf_hc_bed {
        let cmd = format!("bedtools intersect  -a {} -b {} > {}", sf_hc_bed, sf_var, sf_new_var);
        run_cmd_print(&cmd);
    }
}

/// `RsltCMP.cnt_number_of_ins`: number of lines of the file.
pub fn cnt_number_of_ins(sf_input: &str) -> Result<i64> {
    Ok(read_text_lines(sf_input)?.len() as i64)
}

/// `{sample: {rep_type: result file}}`
pub type RsltList = IndexMap<String, IndexMap<String, String>>;

/// `RsltCMP.load_in_rslt_list` (lines: `sf_rslt rep_type sample`; first occurrence wins)
pub fn load_in_rslt_list(sf_list: &str) -> Result<RsltList> {
    let mut m_rslts = RsltList::new();
    for line in read_text_lines(sf_list)? {
        let fields = py_split_ws(&line);
        let sf_rslt = at(&fields, 0)?;
        let s_type = at(&fields, 1)?;
        let s_sample = at(&fields, 2)?;
        m_rslts.entry(s_sample.to_string()).or_default().entry(s_type.to_string()).or_insert_with(|| sf_rslt.to_string());
    }
    Ok(m_rslts)
}

/// Python `repr()` of the `{sample: {rep_type: file}}` dict (printed by the group mode).
pub fn repr_rslt_list(m: &RsltList) -> String {
    let items: Vec<String> = m
        .iter()
        .map(|(k, v)| {
            let inner: Vec<String> = v.iter().map(|(a, b)| format!("{}: {}", py_repr_str(a), py_repr_str(b))).collect();
            format!("{}: {{{}}}", py_repr_str(k), inner.join(", "))
        })
        .collect();
    format!("{{{}}}", items.join(", "))
}

fn print_counts(s_type: &str, n_ori_in: i64, n_hc_in: i64) {
    println!("##########{}############", s_type);
    println!("Original call set size: {}\n", n_ori_in);
    println!("Fall in high confident regions: {}\n", n_hc_in);
    println!("#########################");
}

/// `check_one_sample_agnst_bcmk` (group mode) -> (n_ori_in, n_hc_in, n_hit)
#[allow(clippy::too_many_arguments)]
pub fn check_one_sample_agnst_bcmk(
    sf_input: &str,
    sf_bcmk: &str,
    m_exclude_sites: &Sites,
    sf_hc_bed: Option<&str>,
    i_extd: i64,
    s_type: &str,
    sf_unhit: Option<&str>,
    s_sample: &str,
) -> Result<(i64, i64, i64)> {
    let sf_hc_input = if sf_hc_bed.is_none() { sf_input.to_string() } else { format!("{}.focal_region.txt", sf_input) };
    gnrt_ins_within_focal_region(sf_input, sf_hc_bed, &sf_hc_input);
    let n_ori_in = cnt_number_of_ins(sf_input)?;
    let n_hc_in = cnt_number_of_ins(&sf_hc_input)?;
    print_counts(s_type, n_ori_in, n_hc_in);
    let n_hit = cmp_rslts(&sf_hc_input, sf_bcmk, i_extd, sf_unhit, s_sample, m_exclude_sites)?;
    Ok((n_ori_in, n_hc_in, n_hit))
}

/// Python `sf_exclude is not None and os.path.isfile(sf_exclude)` -> `load_sites`, else `{}`.
fn load_exclude_sites(sf_exclude: Option<&str>) -> Result<Sites> {
    match sf_exclude {
        Some(f) if is_file(f) => load_sites(f),
        _ => Ok(Sites::new()),
    }
}

/// `check_one_sample_agnst_bcmk2` (single-sample mode) -> (n_ori_in, n_hc_in, n_hit)
pub fn check_one_sample_agnst_bcmk2(
    sf_input: &str,
    sf_bcmk: Option<&str>,
    sf_exclude: Option<&str>,
    sf_hc_bed: Option<&str>,
    i_extd: i64,
    s_type: &str,
    sf_unhit: Option<&str>,
) -> Result<(i64, i64, i64)> {
    let sf_hc_input = if sf_hc_bed.is_none() { sf_input.to_string() } else { format!("{}.focal_region.txt", sf_input) };
    gnrt_ins_within_focal_region(sf_input, sf_hc_bed, &sf_hc_input);
    let m_exclude_sites = load_exclude_sites(sf_exclude)?;
    let n_ori_in = cnt_number_of_ins(sf_input)?;
    let n_hc_in = cnt_number_of_ins(&sf_hc_input)?;
    print_counts(s_type, n_ori_in, n_hc_in);
    let n_hit = cmp_rslts2(&sf_hc_input, sf_bcmk, i_extd, sf_unhit, &m_exclude_sites)?;
    Ok((n_ori_in, n_hc_in, n_hit))
}

/// `gnrt_bed_from_input`: `chrm pos pos+1` for every non-header line.
pub fn gnrt_bed_from_input(sf_input: &str, sf_new: &str) -> Result<()> {
    let text = open_text(sf_input)?;
    let mut fout_new = create_file(sf_new)?;
    for line in &text.lines {
        if line.starts_with('#') {
            continue;
        }
        let fields = py_split_ws(line);
        let ins_chrm = at(&fields, 0)?;
        let ins_pos = py_int(at(&fields, 1)?)?;
        fout_new.write_all(format!("{}\t{}\t{}\n", ins_chrm, ins_pos, ins_pos as i128 + 1).as_bytes())?;
    }
    fout_new.flush()?;
    text.check()
}

/// `{sample: {rep_type: (a, b, c)}}`
pub type CntMap = IndexMap<String, IndexMap<String, (i64, i64, i64)>>;

/// `cmp_for_list_of_rslts` -> (m_cnt, l_sample_id, m_tp_fp_fn)
pub fn cmp_for_list_of_rslts(
    sf_list: Option<&str>,
    sf_bcmk: Option<&str>,
    sf_exclude: Option<&str>,
    sf_hc_bed: Option<&str>,
    i_extd: i64,
    sf_unhit: Option<&str>,
) -> Result<(CntMap, Vec<String>, CntMap)> {
    let mut m_cnt = CntMap::new();
    let mut l_sample_id: Vec<String> = Vec::new();
    let mut m_tmp: HashSet<String> = HashSet::new();
    let m_rslt_list = load_in_rslt_list(need(sf_list, "--list")?)?;
    println!("{}", repr_rslt_list(&m_rslt_list));
    let sf_bcmk = need(sf_bcmk, "--bcmk")?;
    let n_bcmk = cnt_number_of_ins(sf_bcmk)?;
    let mut m_tp_fp_fn = CntMap::new();
    let m_exclude_sites = load_exclude_sites(sf_exclude)?;

    for (s_sample, m_types) in &m_rslt_list {
        for (s_rep_type, sf_input) in m_types {
            if !m_tmp.contains(s_sample) {
                l_sample_id.push(s_sample.clone());
                m_tmp.insert(s_sample.clone());
            }
            let sf_new_input = format!("{}.tmp", sf_input);
            gnrt_bed_from_input(sf_input, &sf_new_input)?;
            let (n_ori_in, n_hc_in, n_hit) = check_one_sample_agnst_bcmk(&sf_new_input, sf_bcmk, &m_exclude_sites, sf_hc_bed, i_extd, s_rep_type, sf_unhit, s_sample)?;
            m_cnt.entry(s_sample.clone()).or_default().insert(s_rep_type.clone(), (n_ori_in, n_hc_in, n_hit));

            // used for calc F1 score
            let n_sample_cnt = cnt_number_of_ins(sf_input)?;
            let n_tp = n_hit;
            let n_fp = n_sample_cnt - n_tp;
            let n_fn = n_bcmk - n_tp;
            m_tp_fp_fn.entry(s_sample.clone()).or_default().insert(s_rep_type.clone(), (n_tp, n_fp, n_fn));
        }
    }
    Ok((m_cnt, l_sample_id, m_tp_fp_fn))
}

/// `export_cnt_to_csv_for_point_chart`
pub fn export_cnt_to_csv_for_point_chart(l_sample_info: &[String], m_cnt: &CntMap, sf_cnt: &str) -> Result<()> {
    let mut fout_csv1 = create_file(format!("{}_ori.csv", sf_cnt))?;
    let mut fout_csv2 = create_file(format!("{}_fall_in_hc.csv", sf_cnt))?;
    let mut fout_csv3 = create_file(format!("{}_hit.csv", sf_cnt))?;
    let s_head = "Repeat\tMethod\tTP\tFP\n";
    fout_csv1.write_all(s_head.as_bytes())?;
    fout_csv2.write_all(s_head.as_bytes())?;
    fout_csv3.write_all(s_head.as_bytes())?;
    for s_sample in l_sample_info {
        let m_types = m_cnt.get(s_sample).ok_or_else(|| anyhow!("KeyError: {}", py_repr_str(s_sample)))?;
        for (s_type, &(cnt1, cnt2, cnt3)) in m_types {
            fout_csv1.write_all(format!("{}\t{}\t{}\n", s_type, s_sample, cnt1).as_bytes())?;
            fout_csv2.write_all(format!("{}\t{}\t{}\n", s_type, s_sample, cnt2).as_bytes())?;
            fout_csv3.write_all(format!("{}\t{}\t{}\n", s_type, s_sample, cnt3).as_bytes())?;
        }
    }
    fout_csv1.flush()?;
    fout_csv2.flush()?;
    fout_csv3.flush()?;
    Ok(())
}

/// F1, sensitivity, FDR and the (method, coverage) fields of one record, computed in the
/// order of the script (so that the first failing expression is the same).
struct Scores<'a> {
    f1: f64,
    method: &'a str,
    cov: &'a str,
    sensitivity: f64,
    fdr: f64,
}

fn calc_scores(s_sample: &str, tp: i64, fp: i64, fn_: i64) -> Result<Scores<'_>> {
    let (tp, fp, fn_) = (tp as f64, fp as f64, fn_ as f64);
    let f1 = py_fdiv(2.0 * tp, 2.0 * tp + fp + fn_)?;
    let tmp_fields: Vec<&str> = s_sample.split('_').collect();
    let method = tmp_fields[0];
    let cov = at(&tmp_fields, 1)?;
    let sensitivity = py_fdiv(tp, tp + fn_)?;
    let fdr = py_fdiv(fp, fp + tp)?;
    Ok(Scores { f1, method, cov, sensitivity, fdr })
}

/// `calc_export_F1_score_sensitivity_FDR_csv`
pub fn calc_export_f1_score_sensitivity_fdr_csv(m_tp_fp_fn: &CntMap, sf_out: &str) -> Result<()> {
    let mut fout_f1 = create_file(sf_out)?;
    fout_f1.write_all(b"Method\tCoverage\tF1\tSensitivity\tFDR\tTP\tFP\tFN\n")?;
    for (s_sample, m_types) in m_tp_fp_fn {
        for &(tp, fp, fn_) in m_types.values() {
            let s = calc_scores(s_sample, tp, fp, fn_)?;
            fout_f1.write_all(
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    s.method,
                    s.cov,
                    float_repr(s.f1),
                    float_repr(s.sensitivity),
                    float_repr(s.fdr),
                    tp,
                    fp,
                    fn_
                )
                .as_bytes(),
            )?;
        }
    }
    fout_f1.flush()?;
    Ok(())
}

/// `calc_export_F1_score_sensitivity_FDR_csv2` (writes `sf_out` and `sf_out.2`)
pub fn calc_export_f1_score_sensitivity_fdr_csv2(m_tp_fp_fn: &CntMap, sf_out: &str) -> Result<()> {
    let mut fout_f1 = create_file(sf_out)?;
    let mut fout_f2 = create_file(format!("{}.2", sf_out))?;
    fout_f1.write_all(b"Method\tCoverage\tF1\tSensitivity\tFDR\tCategory\tNumber of TE-Insertion\n")?;
    fout_f2.write_all(b"Platform\tRatio\tCategory\n")?;
    for (s_sample, m_types) in m_tp_fp_fn {
        for &(tp, fp, fn_) in m_types.values() {
            let s = calc_scores(s_sample, tp, fp, fn_)?;
            let prefix = format!("{}\t{}\t{}\t{}\t{}", s.method, s.cov, float_repr(s.f1), float_repr(s.sensitivity), float_repr(s.fdr));
            fout_f1.write_all(format!("{}\tTP\t{}\n", prefix, tp).as_bytes())?;
            fout_f1.write_all(format!("{}\tFP\t{}\n", prefix, fp).as_bytes())?;
            fout_f2.write_all(format!("{}\t{}\tRecall\n", s.method, float_repr(s.sensitivity)).as_bytes())?;
            fout_f2.write_all(format!("{}\t{}\t1-FDR\n", s.method, float_repr(1.0 - s.fdr)).as_bytes())?;
        }
    }
    fout_f1.flush()?;
    fout_f2.flush()?;
    Ok(())
}

/// `get_joint_info`: lines of `sf_ori` whose (chrm, pos) *strings* are listed in
/// `sf_slct_sites`.
pub fn get_joint_info(sf_ori: &str, sf_slct_sites: &str, sf_joint: &str) -> Result<()> {
    let mut m_slct: HashMap<String, HashSet<String>> = HashMap::new();
    for line in read_text_lines(sf_slct_sites)? {
        let fields = py_split_ws(&line);
        let chrm = at(&fields, 0)?;
        let pos = at(&fields, 1)?;
        m_slct.entry(chrm.to_string()).or_default().insert(pos.to_string());
    }
    let text = open_text(sf_ori)?;
    let mut fout = create_file(sf_joint)?;
    for line in &text.lines {
        let fields = py_split_ws(line);
        let chrm = at(&fields, 0)?;
        let pos = at(&fields, 1)?;
        if m_slct.get(chrm).is_some_and(|s| s.contains(pos)) {
            fout.write_all(line.as_bytes())?;
        }
    }
    fout.flush()?;
    text.check()
}

/// `get_joint_info2`: lines of `sf_ori` with a listed site of the same chromosome within
/// `[pos - i_slack, pos + i_slack)`.
pub fn get_joint_info2(sf_ori: &str, sf_slct_sites: &str, i_slack: i64, sf_joint: &str) -> Result<()> {
    let mut m_slct: HashMap<String, BTreeSet<i64>> = HashMap::new();
    for line in read_text_lines(sf_slct_sites)? {
        let fields = py_split_ws(&line);
        let chrm = at(&fields, 0)?;
        let pos = py_int(at(&fields, 1)?)?;
        m_slct.entry(chrm.to_string()).or_default().insert(pos);
    }
    let text = open_text(sf_ori)?;
    let mut fout = create_file(sf_joint)?;
    for line in &text.lines {
        let fields = py_split_ws(line);
        let chrm = at(&fields, 0)?;
        let pos = py_int(at(&fields, 1)?)?;
        if let Some(keys) = m_slct.get(chrm) {
            if super::first_key_in_window(keys, pos, i_slack).is_some() {
                fout.write_all(line.as_bytes())?;
            }
        }
    }
    fout.flush()?;
    text.check()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sites(v: &[(&str, i64)]) -> Sites {
        let mut s = Sites::new();
        for (c, p) in v {
            s.entry(c.to_string()).or_default().entry(*p);
        }
        s
    }

    #[test]
    fn chrm_name_normalisation() {
        assert_eq!(process_chrm_name("chr1", true), "chr1");
        assert_eq!(process_chrm_name("1", true), "chr1");
        assert_eq!(process_chrm_name("chr1", false), "1");
        assert_eq!(process_chrm_name("X", false), "X");
        // "chr" alone is not considered prefixed (len must be > 3)
        assert_eq!(process_chrm_name("chr", true), "chrchr");
        assert_eq!(process_chrm_name("chr", false), "chr");
    }

    #[test]
    fn cmp_dict_window_and_naming() {
        // benchmark without "chr1": its "chr2" is looked up as "2"
        let calls = sites(&[("2", 1000), ("2", 1100), ("1", 50)]);
        let bcmk = sites(&[("chr2", 1150), ("1", 250), ("1", 249), ("chr2", 800), ("chr2", 1300)]);
        let r = cmp_dict(&calls, &bcmk, 200);
        // 1150 -> window [950, 1350): smallest hit 1000; 800 -> [600, 1000): no hit;
        // 1300 -> [1100, 1500): 1100; 250 -> [50, 450): 50; 249 -> [49, 449): 50
        assert_eq!(r.n_overlap, 4);
        let hit: Vec<(String, Vec<i64>)> = r.m_hit.iter().map(|(c, s)| (c.clone(), s.iter().copied().collect())).collect();
        assert_eq!(hit, vec![("2".to_string(), vec![1000, 1100]), ("1".to_string(), vec![50])]);
        let unhit: Vec<(String, Vec<i64>)> = r.m_unhit.iter().map(|(c, s)| (c.clone(), s.iter().copied().collect())).collect();
        assert_eq!(unhit, vec![("2".to_string(), vec![800])]);
        let coord: Vec<(String, Vec<i64>)> = r.m_hit_coordinate_bcmk.iter().map(|(c, s)| (c.clone(), s.iter().copied().collect())).collect();
        assert_eq!(coord, vec![("2".to_string(), vec![1150, 1300]), ("1".to_string(), vec![250, 249])]);
        // no slack: nothing matches, not even identical positions
        assert_eq!(cmp_dict(&bcmk, &bcmk, 0).n_overlap, 0);
    }

    #[test]
    fn dict_repr_like_python() {
        let mut m = RsltList::new();
        m.entry("S1_30x".into()).or_default().insert("L1".into(), "a/b.txt".into());
        m.entry("S1_30x".into()).or_default().insert("Alu".into(), "it's.txt".into());
        m.entry("S2".into()).or_default().insert("SVA".into(), "x\"y".into());
        // Python: print({'S1_30x': {'L1': 'a/b.txt', 'Alu': "it's.txt"}, 'S2': {'SVA': 'x"y'}})
        assert_eq!(repr_rslt_list(&m), "{'S1_30x': {'L1': 'a/b.txt', 'Alu': \"it's.txt\"}, 'S2': {'SVA': 'x\"y'}}");
        assert_eq!(repr_rslt_list(&RsltList::new()), "{}");
    }

    #[test]
    fn scores_like_python() {
        // Python: tp=7, fp=3, fn=5 -> f1 = 14/22, sens = 7/12, fdr = 3/10
        let s = calc_scores("xTea_30x", 7, 3, 5).unwrap();
        assert_eq!(float_repr(s.f1), "0.6363636363636364");
        assert_eq!(float_repr(s.sensitivity), "0.5833333333333334");
        assert_eq!(float_repr(s.fdr), "0.3");
        assert_eq!(float_repr(1.0 - s.fdr), "0.7");
        assert_eq!((s.method, s.cov), ("xTea", "30x"));
        assert!(calc_scores("xTea_30x", 0, 0, 0).is_err());
        assert!(calc_scores("noUnderscore", 1, 0, 0).is_err());
        // fp + tp == 0 with a negative fp
        assert!(calc_scores("a_b", 2, -2, 0).is_err());
    }
}
