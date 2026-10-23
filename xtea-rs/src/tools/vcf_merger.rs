//! Port of `x_vcf_merger.py` (stand-alone tool with its own main):
//!
//! * `-S`: merge the L1 / SVA / Alu / HERV VCFs of every sample into `<wfolder>/<sample>.vcf`
//!   (`SampleVCFMerger`);
//! * `-P`: merge the VCFs of many samples into one population VCF with per-population allele
//!   frequencies (`PopVCFMerger`), plus the `.pca` table and the transduction-source side files.
//!
//! The unused `XVCF` class and `VCFMerger.load_rslts_list` of the script are not reachable from
//! its main and are not ported.
//!
//! Quirks of the original that are kept on purpose:
//! * `append_sites_from_vcf`: when a nearby site already exists and the first record wins the
//!   quality comparison, the new record is still added at its own position (the "keep first"
//!   branch falls through), so both records end up in the merged VCF.  The window is
//!   `[pos - w, pos + w)` and the hit is the smallest position in it, including positions whose
//!   record was blanked by an earlier replacement.
//! * When the L1 VCF of a sample is missing the merged VCF gets no header at all (the header is
//!   taken from the first enabled repeat type, whether or not its file exists).
//! * Peak calling iterates a Python `set` of positions; ties and the order of the merged
//!   records follow CPython's set order (reproduced with [`PySet`]).
//! * The `.td_src.tmp_unsorted` / `.td_src.tmp_sorted` files of `-P` are left on disk.

use super::{at, at_back, create_file, need, py_fdiv, py_int, py_rstrip, py_slice_from, py_slice_inner, py_split_ws, read_text_lines, remove_file};
use crate::globals as gv;
use crate::py::fmt::float_repr;
use crate::py::set::PySet;
use crate::rep_type;
use crate::utils::{is_file, run_cmd_small_output};
use anyhow::Result;
use clap::Parser;
use indexmap::IndexMap;
use std::collections::{BTreeSet, HashMap};
use std::io::Write;

/// Command line of `x_vcf_merger.py` (optparse options, same names and defaults).
#[derive(Parser, Debug, Clone, Default)]
#[command(
    name = "xtea_vcf_merger",
    about = "Merge the VCFs (different repeat types) of each sample, or the samples of a population (Rust port of x_vcf_merger.py)",
    infer_long_args = true,
    args_override_self = true
)]
pub struct Opts {
    #[arg(short = 'S', long = "sample", help = "Merge VCFs (different repeats) of each sample")]
    pub sample_merge: bool,
    #[arg(short = 'P', long = "population", help = "Merge VCFs (different repeats) of each sample")]
    pub population: bool,
    #[arg(allow_hyphen_values = true, short = 'i', value_name = "FILE", help = "sample id list file ")]
    pub sample_info: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'x', value_name = "FILE", help = "xTEA result folder list")]
    pub xrslts: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'p', long = "path", help = "Working folder")]
    pub wfolder: Option<String>,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, short = 'n', long = "cores", help = "number of cores")]
    pub cores: Option<i64>,
    #[arg(
        allow_hyphen_values = true,
        value_parser = crate::py::optparse::parse_int,
        short = 'w',
        long = "window",
        default_value_t = 50,
        help = "Maximum window size to merge two sites close to each other"
    )]
    pub window: i64,
    #[arg(allow_hyphen_values = true, short = 'r', long = "ref", help = "reference genome")]
    pub r#ref: Option<String>,
    #[arg(
        allow_hyphen_values = true,
        value_parser = crate::py::optparse::parse_int,
        short = 'y',
        long = "reptype",
        help = "Type of repeats working on: 1-L1, 2-Alu, 4-SVA, 8-HERV, 16-Mitochondrial"
    )]
    pub rep_type: Option<i64>,
    #[arg(allow_hyphen_values = true, short = 'o', long = "output", value_name = "FILE", help = "The output file")]
    pub output: Option<String>,
    /// Positional arguments are accepted and ignored (optparse collects them separately).
    #[arg(hide = true)]
    pub extra_args: Vec<String>,
}

/// The script's `__main__`.
pub fn run(opts: Opts) -> Result<()> {
    let sf_sample_info = opts.sample_info.as_deref();
    let sf_rslt_list = opts.xrslts.as_deref();
    let sf_wfolder = opts.wfolder.as_deref();
    if opts.sample_merge {
        merge_vcf_for_all_rep_types(sf_sample_info, sf_rslt_list, opts.rep_type, opts.window, sf_wfolder)
    } else if opts.population {
        gnrt_population_vcf(sf_sample_info, sf_rslt_list, opts.window, sf_wfolder, opts.output.as_deref())
    } else {
        Ok(())
    }
}

/// sample name -> (sex, biosample id, population code)
pub type SampleInfo = IndexMap<String, (String, String, String)>;

// ------------------------------------------------------------------------------------------
// VCFMerger (shared base class)
// ------------------------------------------------------------------------------------------

/// `load_in_sample_info`: (sample info, number of samples of each population).
pub fn load_in_sample_info(sf_sample_info: Option<&str>) -> Result<(SampleInfo, IndexMap<String, i64>)> {
    let mut m_sample_info = SampleInfo::new();
    let mut m_pop_cnt: IndexMap<String, i64> = IndexMap::new();
    for line in read_text_lines(need(sf_sample_info, "-i")?)? {
        let fields = py_split_ws(&line);
        let sample_name = at(&fields, 0)?;
        let sex = at(&fields, 1)?;
        let biosample_id = at(&fields, 2)?;
        let pop_code = at(&fields, 3)?;
        if pop_code.is_empty() {
            println!("population sample file parse error: {}", py_rstrip(&line));
            continue;
        }
        *m_pop_cnt.entry(pop_code.to_string()).or_insert(0) += 1;
        m_sample_info.insert(sample_name.to_string(), (sex.to_string(), biosample_id.to_string(), pop_code.to_string()));
    }
    Ok((m_sample_info, m_pop_cnt))
}

/// `load_rslts_folders`: sample -> second column (first occurrence wins).
pub fn load_rslts_folders(sf_rslt_folders: Option<&str>) -> Result<IndexMap<String, String>> {
    let mut m_folders: IndexMap<String, String> = IndexMap::new();
    for line in read_text_lines(need(sf_rslt_folders, "-x")?)? {
        let fields = py_split_ws(py_rstrip(&line));
        let s_sample = at(&fields, 0)?;
        let sf_folder = at(&fields, 1)?;
        if !m_folders.contains_key(s_sample) {
            m_folders.insert(s_sample.to_string(), sf_folder.to_string());
        }
    }
    Ok(m_folders)
}

/// `sort_raw_rslt` (through the shell, like `CMD_RUNNER.run_cmd_small_output`).
pub fn sort_raw_rslt(sf_raw_rslt: &str, sf_sorted_rslt: &str) {
    let cmd = format!("sort -k1,1V -k2,2n -o {} {}", sf_sorted_rslt, sf_raw_rslt);
    run_cmd_small_output(&cmd);
}

/// `load_vcf_head`
pub fn load_vcf_head(sf_vcf: &str) -> Result<Vec<String>> {
    let mut l_head = Vec::new();
    if !is_file(sf_vcf) {
        println!("[Get VCF head]: File {} doesn't exist!", sf_vcf);
        return Ok(l_head);
    }
    for line in read_text_lines(sf_vcf)? {
        if line.starts_with('#') {
            if line.contains("contig") && (line.contains("HLA") || line.contains("Un") || line.contains('_')) {
                continue;
            }
            l_head.push(py_rstrip(&line).to_string());
        }
    }
    Ok(l_head)
}

// ------------------------------------------------------------------------------------------
// SampleVCFMerger
// ------------------------------------------------------------------------------------------

/// Sites of one chromosome: `{pos: record}` in insertion order plus the sorted key set used for
/// the window look-ups (keys are never removed; replaced records are blanked with "").
#[derive(Debug, Default, Clone)]
pub struct ChrSites {
    pub map: IndexMap<i64, String>,
    keys: BTreeSet<i64>,
}

impl ChrSites {
    fn set(&mut self, pos: i64, rcd: String) {
        self.keys.insert(pos);
        self.map.insert(pos, rcd);
    }
}

/// `{chrm: {pos: record}}`
pub type SampleSites = IndexMap<String, ChrSites>;

/// `merge_vcf_for_all_rep_types`
pub fn merge_vcf_for_all_rep_types(
    sf_sample_info: Option<&str>,
    sf_rslt_list: Option<&str>,
    i_rep_type: Option<i64>,
    i_single_win: i64,
    sf_wfolder: Option<&str>,
) -> Result<()> {
    let (m_sample_info, _m_pop_cnt) = load_in_sample_info(sf_sample_info)?;
    let m_rslt_folders = load_rslts_folders(sf_rslt_list)?;
    let mut sf_wfolder = need(sf_wfolder, "-p")?.to_string();
    if sf_wfolder.is_empty() {
        return Ok(());
    }
    if !sf_wfolder.ends_with('/') {
        sf_wfolder.push('/');
    }
    for s_sample in m_sample_info.keys() {
        if let Some(folder) = m_rslt_folders.get(s_sample) {
            let mut sf_rslt_folder = folder.clone();
            if !sf_rslt_folder.is_empty() && !sf_rslt_folder.ends_with('/') {
                sf_rslt_folder.push('/');
            }
            let sf_merged = format!("{}{}.vcf", sf_wfolder, s_sample);
            merge_single_sample_ins(i_rep_type, &sf_rslt_folder, s_sample, i_single_win, &sf_merged)?;
        }
    }
    Ok(())
}

/// `merge_single_sample_ins`: L1 first, then SVA, Alu and HERV (appended with conflict check).
pub fn merge_single_sample_ins(i_rep_type: Option<i64>, sf_rslt_folder: &str, s_sample_id: &str, i_win: i64, sf_merged: &str) -> Result<()> {
    println!("Merging sub-type vcf of sample {}", s_sample_id);
    // `None & 1` raises a TypeError in the script when -y is not given
    let i_rep_type = i_rep_type.ok_or_else(|| anyhow::anyhow!("TypeError: unsupported operand type(s) for &: 'NoneType' and 'int' (-y not given)"))?;
    let m_rep_types = rep_type::parse_out_all_rep_type(i_rep_type);
    let has = |s: &str| m_rep_types.contains(&s);
    let mut m_all_sites = SampleSites::new();
    let mut b_head = false;
    let mut l_head: Vec<String> = Vec::new();
    // L1
    if has(rep_type::REP_TYPE_L1) {
        let sf_vcf = format!("{}{}/{}_{}.vcf", sf_rslt_folder, "L1", s_sample_id, rep_type::REP_TYPE_L1);
        l_head = load_vcf_head(&sf_vcf)?;
        b_head = true;
        load_in_sites_from_vcf(&sf_vcf, &mut m_all_sites)?;
    }
    // SVA, Alu, HERV
    for (s_type, s_folder) in [(rep_type::REP_TYPE_SVA, rep_type::REP_TYPE_SVA), (rep_type::REP_TYPE_ALU, "Alu"), (rep_type::REP_TYPE_HERV, rep_type::REP_TYPE_HERV)] {
        if has(s_type) {
            let sf_vcf = format!("{}{}/{}_{}.vcf", sf_rslt_folder, s_folder, s_sample_id, s_type);
            if !b_head {
                l_head = load_vcf_head(&sf_vcf)?;
                b_head = true;
            }
            append_sites_from_vcf(&sf_vcf, i_win, &mut m_all_sites)?;
        }
    }
    let sf_unsorted = format!("{}.unsorted", sf_merged);
    dump_to_unsorted_file(&m_all_sites, &sf_unsorted)?;
    let sf_sorted = format!("{}.sorted", sf_merged);
    sort_raw_rslt(&sf_unsorted, &sf_sorted);
    {
        let mut fout_merged = create_file(sf_merged)?;
        for sline in &l_head {
            fout_merged.write_all(py_rstrip(sline).as_bytes())?;
            fout_merged.write_all(b"\n")?;
        }
        for line in read_text_lines(&sf_sorted)? {
            fout_merged.write_all(py_rstrip(&line).as_bytes())?;
            fout_merged.write_all(b"\n")?;
        }
        fout_merged.flush()?;
    }
    remove_file(&sf_unsorted)?;
    remove_file(&sf_sorted)?;
    Ok(())
}

/// `dump_to_unsorted_file` (blanked records are skipped)
pub fn dump_to_unsorted_file(m_sites: &SampleSites, sf_tmp: &str) -> Result<()> {
    let mut fout_tmp = create_file(sf_tmp)?;
    for chr_sites in m_sites.values() {
        for sinfo in chr_sites.map.values() {
            if !sinfo.is_empty() {
                fout_tmp.write_all(sinfo.as_bytes())?;
                fout_tmp.write_all(b"\n")?;
            }
        }
    }
    fout_tmp.flush()?;
    Ok(())
}

/// `load_in_sites_from_vcf`
pub fn load_in_sites_from_vcf(sf_vcf: &str, m_sites: &mut SampleSites) -> Result<()> {
    if !is_file(sf_vcf) {
        println!("File {} doesn't exist!", sf_vcf);
        return Ok(());
    }
    for line in read_text_lines(sf_vcf)? {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 || line.starts_with('#') {
            continue;
        }
        let ins_pos = py_int(fields[1])?;
        let chr_sites = m_sites.entry(fields[0].to_string()).or_default();
        chr_sites.set(ins_pos, py_rstrip(&line).to_string());
    }
    Ok(())
}

/// `append_sites_from_vcf` (see the module documentation for the kept quirks)
pub fn append_sites_from_vcf(sf_vcf: &str, i_win: i64, m_sites: &mut SampleSites) -> Result<()> {
    if !is_file(sf_vcf) {
        println!("File {} doesn't exist!!", sf_vcf);
        return Ok(());
    }
    for line in read_text_lines(sf_vcf)? {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 || line.starts_with('#') {
            continue;
        }
        let ins_pos = py_int(fields[1])?;
        let chr_sites = m_sites.entry(fields[0].to_string()).or_default();
        let s_rcd = py_rstrip(&line);
        if let Some(i_hit_pos) = super::first_key_in_window(&chr_sites.keys, ins_pos, i_win) {
            let b_keep_first = cmp_quality_of_two_ins(&chr_sites.map[&i_hit_pos], s_rcd);
            if !b_keep_first {
                chr_sites.set(i_hit_pos, String::new());
                chr_sites.set(ins_pos, s_rcd.to_string());
                continue;
            }
        }
        chr_sites.set(ins_pos, s_rcd.to_string());
    }
    Ok(())
}

/// `cmp_quality_of_two_ins`: false only when the second record is a two-side-TPRT call and the
/// first is not.
pub fn cmp_quality_of_two_ins(s_rcd1: &str, s_rcd2: &str) -> bool {
    let s_tprt_both = gv::TWO_SIDE_TPRT_BOTH;
    !(!s_rcd1.contains(s_tprt_both) && s_rcd2.contains(s_tprt_both))
}

// ------------------------------------------------------------------------------------------
// PopVCFMerger
// ------------------------------------------------------------------------------------------

/// (sample id, VCF record)
pub type SampleRcd = (String, String);
/// `{chrm: {pos: [(sample, record), ...]}}`
pub type PopSites = IndexMap<String, IndexMap<i64, Vec<SampleRcd>>>;
/// `{population: {transduction source: count}}`
pub type TdSrc = IndexMap<String, IndexMap<String, i64>>;

/// `cnt_current_pop_cnt`: number of listed samples of each population (in result-list order).
pub fn cnt_current_pop_cnt(m_sample_info: &SampleInfo, m_rslt_list: &IndexMap<String, String>) -> IndexMap<String, i64> {
    let mut m_pop_cnt: IndexMap<String, i64> = IndexMap::new();
    for sample_id in m_rslt_list.keys() {
        if let Some(info) = m_sample_info.get(sample_id) {
            *m_pop_cnt.entry(info.2.clone()).or_insert(0) += 1;
        }
    }
    m_pop_cnt
}

/// `gnrt_population_vcf`
pub fn gnrt_population_vcf(
    sf_sample_info: Option<&str>,
    sf_rslt_list: Option<&str>,
    peak_window: i64,
    sf_wfolder: Option<&str>,
    sf_out: Option<&str>,
) -> Result<()> {
    let (m_sample_info, _m_pop_cnt_total) = load_in_sample_info(sf_sample_info)?;
    let m_rslt_list = load_rslts_folders(sf_rslt_list)?;
    let m_pop_cnt = cnt_current_pop_cnt(&m_sample_info, &m_rslt_list);
    if need(sf_wfolder, "-p")?.is_empty() {
        return Ok(());
    }

    let mut m_all_sites = PopSites::new();
    let mut sf_template_vcf = String::new();
    let mut b_set_template = false;
    for s_sample in m_sample_info.keys() {
        if let Some(sf_sample_vcf) = m_rslt_list.get(s_sample) {
            if !b_set_template {
                sf_template_vcf = sf_sample_vcf.clone();
                b_set_template = true;
            }
            if !is_file(sf_sample_vcf) {
                println!("{} doesn't exist!", sf_sample_vcf);
                continue;
            }
            load_in_sample_vcf(s_sample, sf_sample_vcf, &mut m_all_sites)?;
        }
    }

    let m_peak_sites = call_peak_candidate_sites_vcf(&m_all_sites, peak_window);
    let mut l_samples: Vec<String> = m_rslt_list.keys().cloned().collect();
    l_samples.sort();

    let mut s_column_name = String::from("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for sample_name in &l_samples {
        s_column_name.push('\t');
        s_column_name.push_str(sample_name);
    }
    s_column_name.push('\n');

    let sf_out = need(sf_out, "-o")?;
    let sf_unsorted = format!("{}.tmp_unsorted", sf_out);
    let sf_td_src_unsort = format!("{}.td_src.tmp_unsorted", sf_out);
    {
        let mut fout = create_file(&sf_unsorted)?;
        let mut fout_td_src = create_file(&sf_td_src_unsort)?;
        for (chrm, m_pos) in &m_peak_sites {
            for (pos, l_sample_rcds) in m_pos {
                let (s_info, m_td_src_rcd, s_pos_info) = cvt_one_rcd_to_pop_info(chrm, *pos, &m_sample_info, &l_samples, &m_pop_cnt, l_sample_rcds)?;
                fout.write_all(s_info.as_bytes())?;
                let s_td_src = cvt_td_src_rcd_to_str(&m_td_src_rcd);
                if !s_td_src.is_empty() {
                    fout_td_src.write_all(s_pos_info.as_bytes())?;
                    fout_td_src.write_all(s_td_src.as_bytes())?;
                    fout_td_src.write_all(b"\n")?;
                }
            }
        }
        fout.flush()?;
        fout_td_src.flush()?;
    }

    let sf_sorted = format!("{}.tmp_sorted", sf_out);
    sort_raw_rslt(&sf_unsorted, &sf_sorted);
    let sf_td_src_sort = format!("{}.td_src.tmp_sorted", sf_out);
    sort_raw_rslt(&sf_td_src_unsort, &sf_td_src_sort);

    dump_info_for_pca(&m_sample_info, &l_samples, &m_peak_sites, &format!("{}.pca", sf_out))?;
    {
        let mut fout = create_file(sf_out)?;
        for tmp_head_line in take_chrm_head_from_template_vcf(&sf_template_vcf)? {
            fout.write_all(tmp_head_line.as_bytes())?;
            fout.write_all(b"\n")?;
        }
        fout.write_all(prepare_head_other(&m_pop_cnt).as_bytes())?;
        fout.write_all(s_column_name.as_bytes())?;
        for line in read_text_lines(&sf_sorted)? {
            fout.write_all(line.as_bytes())?;
        }
        fout.flush()?;
    }
    remove_file(&sf_unsorted)?;
    remove_file(&sf_sorted)?;
    Ok(())
}

/// `dump_info_for_pca`: one row per listed sample with 0/1/2 per merged site.
pub fn dump_info_for_pca(m_sample_info: &SampleInfo, l_samples: &[String], m_sites: &PopSites, sf_out: &str) -> Result<()> {
    let mut l_sites: Vec<String> = Vec::new();
    let mut m_sample_sites: HashMap<&str, HashMap<String, String>> = HashMap::new();
    let mut s_head = String::from("sample,population");
    for (chrm, m_pos) in m_sites {
        for (pos, l_rcds) in m_pos {
            let s_site = format!("{}_{}", chrm, pos);
            s_head.push(',');
            s_head.push_str(&s_site);
            for (s_sample, site_info) in l_rcds {
                let site_fields = py_split_ws(py_rstrip(site_info));
                let sgntp = at_back(&site_fields, 1)?;
                m_sample_sites.entry(s_sample.as_str()).or_default().insert(s_site.clone(), sgntp.to_string());
            }
            l_sites.push(s_site);
        }
    }
    s_head.push('\n');

    let mut fout = create_file(sf_out)?;
    fout.write_all(s_head.as_bytes())?;
    for sample in l_samples {
        let sites = match m_sample_sites.get(sample.as_str()) {
            Some(s) => s,
            None => {
                println!("{} not in sample list", sample);
                continue;
            }
        };
        let s_pop = &m_sample_info
            .get(sample)
            .ok_or_else(|| anyhow::anyhow!("KeyError: {}", super::py_repr_str(sample)))?
            .2;
        let mut sinfo = format!("{},{}", sample, s_pop);
        for s_site in &l_sites {
            match sites.get(s_site) {
                Some(g) if g == "1/1" || g == "1|1" => sinfo.push_str(",2"),
                Some(_) => sinfo.push_str(",1"),
                None => sinfo.push_str(",0"),
            }
        }
        sinfo.push('\n');
        fout.write_all(sinfo.as_bytes())?;
    }
    fout.flush()?;
    Ok(())
}

/// `_take_chrm_head_from_template_vcf`: the ##fileformat / ##reference / ##contig lines.
pub fn take_chrm_head_from_template_vcf(sf_vcf: &str) -> Result<Vec<String>> {
    let mut l_heads = Vec::new();
    for line in read_text_lines(sf_vcf)? {
        if line.starts_with('#') {
            if line.contains("contig") && (line.contains("HLA") || line.contains("Un") || line.contains('_')) {
                continue;
            }
            let field0 = line.split('=').next().unwrap_or("");
            if field0 == "##fileformat" || field0 == "##reference" || field0 == "##contig" {
                l_heads.push(py_rstrip(&line).to_string());
            }
        }
    }
    Ok(l_heads)
}

/// `_cvt_td_src_rcd_to_str`: the most frequent transduction source of each population.
pub fn cvt_td_src_rcd_to_str(m_td_src_rcd: &TdSrc) -> String {
    let mut s_td_src_info = String::new();
    for (s_pop_code, m_src) in m_td_src_rcd {
        let mut i_max_cnt = 0;
        let mut s_max_td_src = "";
        for (s_td_src, &cnt) in m_src {
            if cnt > i_max_cnt {
                s_max_td_src = s_td_src;
                i_max_cnt = cnt;
            }
        }
        if !s_max_td_src.is_empty() {
            s_td_src_info.push('\t');
            s_td_src_info.push_str(s_pop_code);
            s_td_src_info.push(':');
            s_td_src_info.push_str(s_max_td_src);
        }
    }
    s_td_src_info
}

/// The key with the highest count (the first one on ties); "" when no count is positive.
fn most_frequent(m: &IndexMap<String, i64>) -> &str {
    let mut n_max = 0;
    let mut s_best = "";
    for (k, &v) in m {
        if n_max < v {
            n_max = v;
            s_best = k;
        }
    }
    s_best
}

fn count(m: &mut IndexMap<String, i64>, key: &str) {
    *m.entry(key.to_string()).or_insert(0) += 1;
}

/// `_cvt_one_rcd_to_pop_info`: one merged site -> (VCF line, transduction sources by
/// population, the leading "CHROM..INFO(SVTYPE)" part).
pub fn cvt_one_rcd_to_pop_info(
    chrm: &str,
    pos: i64,
    m_sample_info: &SampleInfo,
    l_samples: &[String],
    m_pop_cnt: &IndexMap<String, i64>,
    l_sample_rcd: &[SampleRcd],
) -> Result<(String, TdSrc, String)> {
    let mut m_pop_cnt_per_site: IndexMap<String, i64> = IndexMap::new();
    let mut m_sample_gntp: HashMap<&str, &str> = HashMap::new();
    let mut m_het: IndexMap<String, i64> = IndexMap::new();
    let mut m_homo: IndexMap<String, i64> = IndexMap::new();
    let mut m_site_type: IndexMap<String, i64> = IndexMap::new();
    let mut m_end_pos: IndexMap<String, i64> = IndexMap::new();
    let mut m_genes: IndexMap<String, i64> = IndexMap::new();
    let mut m_ref_rep: IndexMap<String, i64> = IndexMap::new();
    let mut m_sv_len: IndexMap<String, i64> = IndexMap::new();
    let mut m_td_src = TdSrc::new();
    for (sample_id, s_gvcf_rcd) in l_sample_rcd {
        let info = match m_sample_info.get(sample_id) {
            Some(i) => i,
            None => {
                println!("[Convert Error]: sample {} not in dictionary", sample_id);
                continue;
            }
        };
        let l_tmp_fields = py_split_ws(s_gvcf_rcd);
        let s_rep_type = at(&l_tmp_fields, 4)?;
        count(&mut m_site_type, s_rep_type);
        let s_gntp = at_back(&l_tmp_fields, 1)?;
        m_sample_gntp.insert(sample_id.as_str(), s_gntp);
        let s_tmp_fields2: Vec<&str> = at(&l_tmp_fields, 7)?.split(';').collect();
        let s_sv_len = at(&s_tmp_fields2, 1)?; // note fixed index
        count(&mut m_sv_len, s_sv_len);
        let s_end = at(&s_tmp_fields2, 2)?;
        count(&mut m_end_pos, s_end);

        let s_pop_code = info.2.as_str();
        let s_td_src = py_slice_from(at(&s_tmp_fields2, 6)?, 7); // note fixed index
        let m_pop_td = m_td_src.entry(s_pop_code.to_string()).or_default();
        if !s_td_src.contains(gv::NOT_TRANSDUCTION) {
            *m_pop_td.entry(s_td_src.to_string()).or_insert(0) += 1;
        }
        let s_ref_rep = at_back(&s_tmp_fields2, 2)?;
        count(&mut m_ref_rep, s_ref_rep);
        let s_gene = at_back(&s_tmp_fields2, 1)?;
        count(&mut m_genes, s_gene);

        count(&mut m_pop_cnt_per_site, s_pop_code);
        match s_gntp {
            "0/1" | "1/0" | "0|1" | "1|0" => count(&mut m_het, s_pop_code),
            "1|1" | "1/1" => count(&mut m_homo, s_pop_code),
            _ => {}
        }
    }

    let mut sinfo = String::new();
    let s_type = most_frequent(&m_site_type);
    let s_position_info = format!("{}\t{}\t.\tN\t{}\t.\tPASS\tSVTYPE={};", chrm, pos, s_type, py_slice_inner(s_type));
    sinfo.push_str(&s_position_info);

    // NOTE (original): the most frequent SVLEN is used, not the longest one
    let s_max_end = most_frequent(&m_end_pos);
    let s_max_sv_len = most_frequent(&m_sv_len);
    sinfo.push_str(s_max_end);
    sinfo.push(';');
    sinfo.push_str(s_max_sv_len);
    sinfo.push(';');

    let mut s_samples_gntp = String::new();
    for sample_id in l_samples {
        s_samples_gntp.push('\t');
        s_samples_gntp.push_str(m_sample_gntp.get(sample_id.as_str()).copied().unwrap_or("0/0"));
    }

    let mut n_total_het: i64 = 0;
    let mut n_total_homo: i64 = 0;
    let mut s_tmp_allele_info = String::new();
    for (s_pop_code, &n_pop_samples) in m_pop_cnt {
        if !m_pop_cnt_per_site.contains_key(s_pop_code) {
            continue;
        }
        let n_pop_an = 2 * n_pop_samples;
        let n_het = m_het.get(s_pop_code).copied().unwrap_or(0);
        let n_homo = m_homo.get(s_pop_code).copied().unwrap_or(0);
        let n_pop_ac = n_het + 2 * n_homo;
        let n_pop_af = py_fdiv(n_pop_ac as f64, n_pop_an as f64)?;
        n_total_het += n_het;
        n_total_homo += n_homo;
        let n_homo_ref = n_pop_samples - n_het - n_homo;
        let n_het_freq = py_fdiv(n_het as f64, n_pop_samples as f64)?;
        let n_homo_freq = py_fdiv(n_homo as f64, n_pop_samples as f64)?;
        let n_homoref_freq = py_fdiv(n_homo_ref as f64, n_pop_samples as f64)?;
        let rcd = AlleleRcd {
            an: n_pop_an,
            ac: n_pop_ac,
            af: n_pop_af,
            n_het,
            n_homo,
            n_homo_ref,
            het_freq: n_het_freq,
            homo_freq: n_homo_freq,
            homo_ref_freq: n_homoref_freq,
            n_samples: n_pop_samples,
        };
        s_tmp_allele_info.push_str(&prepare_one_record(&format!("{}_", s_pop_code), &rcd));
    }

    let n_total_sample = l_samples.len() as i64; // total number of counted samples
    let n_total_an = 2 * n_total_sample;
    let n_total_ac = n_total_het + 2 * n_total_homo;
    let n_total_af = py_fdiv(n_total_ac as f64, n_total_an as f64)?;
    let n_total_homoref = n_total_sample - n_total_het - n_total_homo;
    let n_total_het_freq = py_fdiv(n_total_het as f64, n_total_sample as f64)?;
    let n_total_homo_freq = py_fdiv(n_total_homo as f64, n_total_sample as f64)?;
    let n_total_homo_ref_freq = py_fdiv(n_total_homoref as f64, n_total_sample as f64)?;
    let total_rcd = AlleleRcd {
        an: n_total_an,
        ac: n_total_ac,
        af: n_total_af,
        n_het: n_total_het,
        n_homo: n_total_homo,
        n_homo_ref: n_total_homoref,
        het_freq: n_total_het_freq,
        homo_freq: n_total_homo_freq,
        homo_ref_freq: n_total_homo_ref_freq,
        n_samples: n_total_sample,
    };
    sinfo.push_str(&prepare_one_record("", &total_rcd)); // total AF
    sinfo.push_str(&s_tmp_allele_info); // sub-population AF

    sinfo.push_str(most_frequent(&m_ref_rep));
    sinfo.push(';');
    sinfo.push_str(most_frequent(&m_genes));
    sinfo.push_str("\tGT");
    sinfo.push_str(&s_samples_gntp);
    sinfo.push('\n');
    Ok((sinfo, m_td_src, s_position_info))
}

/// The allele record tuple of `_cvt_one_rcd_to_pop_info`.
#[derive(Debug, Clone, Copy)]
pub struct AlleleRcd {
    pub an: i64,
    pub ac: i64,
    pub af: f64,
    pub n_het: i64,
    pub n_homo: i64,
    pub n_homo_ref: i64,
    pub het_freq: f64,
    pub homo_freq: f64,
    pub homo_ref_freq: f64,
    pub n_samples: i64,
}

/// `_prepare_one_record` (floats printed with Python's `str(float)`).
pub fn prepare_one_record(p: &str, r: &AlleleRcd) -> String {
    format!(
        "{p}AN={};{p}AC={};{p}AF={};{p}N_HET={};{p}N_HOMALT={};{p}N_HOMREF={};{p}FREQ_HET={};{p}FREQ_HOMALT={};{p}FREQ_HOMREF={};{p}N_BI_GENOS={};",
        r.an,
        r.ac,
        float_repr(r.af),
        r.n_het,
        r.n_homo,
        r.n_homo_ref,
        float_repr(r.het_freq),
        float_repr(r.homo_freq),
        float_repr(r.homo_ref_freq),
        r.n_samples,
        p = p
    )
}

/// `_gnrt_population_head_info_one_pop`
fn gnrt_population_head_info_one_pop(c1: &str, c2: &str, l_head_info: &mut Vec<String>) {
    l_head_info.push(format!("##INFO=<ID={c1}AC,Number=A,Type=Integer,Description=\"Number of non-reference {c2} alleles observed.\">\n"));
    l_head_info.push(format!("##INFO=<ID={c1}AF,Number=A,Type=Float,Description=\"{c2} allele frequency.\">\n"));
    l_head_info.push(format!("##INFO=<ID={c1}AN,Number=1,Type=Integer,Description=\"Total number of {c2} alleles genotyped.\">\n"));
    l_head_info.push(format!("##INFO=<ID={c1}FREQ_HET,Number=1,Type=Float,Description=\"{c2} heterozygous genotype frequency.\">\n"));
    l_head_info.push(format!("##INFO=<ID={c1}FREQ_HOMALT,Number=1,Type=Float,Description=\"{c2} homozygous alternate genotype frequency.\">\n"));
    l_head_info.push(format!("##INFO=<ID={c1}FREQ_HOMREF,Number=1,Type=Float,Description=\"{c2} homozygous reference genotype frequency.\">\n"));
    l_head_info.push(format!("##INFO=<ID={c1}N_BI_GENOS,Number=1,Type=Integer,Description=\"Total number of {c2} individuals with complete genotypes.\">\n"));
    l_head_info.push(format!("##INFO=<ID={c1}N_HET,Number=1,Type=Integer,Description=\"Number of {c2} individuals with heterozygous genotypes.\">\n"));
    l_head_info.push(format!("##INFO=<ID={c1}N_HOMALT,Number=1,Type=Integer,Description=\"Number of {c2} individuals with homozygous alternate genotypes.\">\n"));
    l_head_info.push(format!("##INFO=<ID={c1}N_HOMREF,Number=1,Type=Integer,Description=\"Number of {c2} individuals with homozygous reference genotypes.\">\n"));
}

/// `prepare_head_other`: ALT/FILTER/INFO/FORMAT header lines of the population VCF.
pub fn prepare_head_other(m_pop_cnt: &IndexMap<String, i64>) -> String {
    let mut sinfo = String::new();
    sinfo.push_str("##ALT=<ID=INS:ME:SVA,Description=\"Insertion of SVA element\">\n");
    sinfo.push_str("##ALT=<ID=INS:ME:LINE1,Description=\"Insertion of LINE1 element\">\n");
    sinfo.push_str("##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">\n");
    sinfo.push_str("##ALT=<ID=INS:ME:HERV-K,Description=\"Insertion of HERV-K element\">\n");
    sinfo.push_str("##ALT=<ID=INS:MT,Description=\"Mitochondrial insertion\">\n");
    sinfo.push_str("##ALT=<ID=INS:PSDGN,Description=\"Pseudogene insertion\">\n");
    sinfo.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
    sinfo.push_str("##FILTER=<ID=LowConfident,Description=\"All filters passed\">\n");
    let mut l_pop_head_info = Vec::new();
    gnrt_population_head_info_one_pop("", "", &mut l_pop_head_info); // total
    for s_code in m_pop_cnt.keys() {
        gnrt_population_head_info_one_pop(&format!("{}_", s_code), s_code, &mut l_pop_head_info);
    }
    for s in l_pop_head_info {
        sinfo.push_str(&s);
    }
    sinfo.push_str("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
    sinfo.push_str("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Insertion length\">\n");
    sinfo.push_str("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of this variant\">\n");
    sinfo.push_str("##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Insertion orientation (+/- for sense/antisense)\">\n");
    sinfo.push_str("##INFO=<ID=REF_REP,Number=1,Type=String,Description=\"Fall in reference repeat copy or not\">\n");
    sinfo.push_str("##INFO=<ID=GENE_INFO,Number=1,Type=String,Description=\"Fall in gene region\">\n");
    sinfo.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    sinfo
}

/// `load_in_sample_vcf`
pub fn load_in_sample_vcf(sample_id: &str, sf_vcf: &str, m_all_sites: &mut PopSites) -> Result<()> {
    if !is_file(sf_vcf) {
        println!("File {} doesn't exist!", sf_vcf);
        return Ok(());
    }
    for line in read_text_lines(sf_vcf)? {
        let s_rcd = py_rstrip(&line);
        let fields: Vec<&str> = s_rcd.split('\t').collect();
        if fields.len() < 3 || line.starts_with('#') {
            continue;
        }
        let ins_pos = py_int(fields[1])?;
        m_all_sites
            .entry(fields[0].to_string())
            .or_default()
            .entry(ins_pos)
            .or_default()
            .push((sample_id.to_string(), s_rcd.to_string()));
    }
    Ok(())
}

/// `call_peak_candidate_sites_vcf`: cluster the sorted positions of each chromosome (a gap
/// larger than `peak_window` starts a new cluster); every cluster is reported at its position
/// with the most records (first one in CPython set order on ties) and carries all records of
/// the cluster, in set order.
pub fn call_peak_candidate_sites_vcf(m_candidate_sites: &PopSites, peak_window: i64) -> PopSites {
    let mut m_peak_candidate_sites = PopSites::new();
    for (chrm, m_pos) in m_candidate_sites {
        let mut l_pos: Vec<i64> = m_pos.keys().copied().collect();
        l_pos.sort();
        let flush = |set_cluster: &PySet<i64>, out: &mut PopSites| {
            let mut max_cnt = 0usize;
            let mut tmp_candidate_pos: i64 = 0;
            let mut l_tmp_rcds: Vec<SampleRcd> = Vec::new();
            for &tmp_pos in set_cluster.iter() {
                let rcds = &m_pos[&tmp_pos];
                l_tmp_rcds.extend(rcds.iter().cloned());
                if max_cnt < rcds.len() {
                    tmp_candidate_pos = tmp_pos;
                    max_cnt = rcds.len();
                }
            }
            out.entry(chrm.clone()).or_default().entry(tmp_candidate_pos).or_default().extend(l_tmp_rcds);
        };
        let mut pre_pos: i64 = -1;
        let mut set_cluster: PySet<i64> = PySet::new();
        for pos in l_pos {
            if pre_pos == -1 {
                pre_pos = pos;
                set_cluster.add(pre_pos);
                continue;
            }
            if pos as i128 - pre_pos as i128 > peak_window as i128 {
                flush(&set_cluster, &mut m_peak_candidate_sites);
                set_cluster.clear();
            }
            pre_pos = pos;
            set_cluster.add(pre_pos);
        }
        // push out the last group
        flush(&set_cluster, &mut m_peak_candidate_sites);
    }
    m_peak_candidate_sites
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rcd(sample: &str, line: &str) -> SampleRcd {
        (sample.to_string(), line.to_string())
    }

    #[test]
    fn peak_calling_follows_cpython_set_order() {
        // Expected values printed by CPython 3.8 running PopVCFMerger.call_peak_candidate_sites_vcf
        // on the same input (positions inserted in this order, one record per (sample, pos)).
        let mut m: PopSites = PopSites::new();
        let ins = [(1015, "a"), (1007, "b"), (1008, "c"), (1008, "d"), (1007, "e"), (1015, "f"), (5000, "g"), (5040, "h"), (5100, "i"), (9, "j"), (1, "k")];
        for (p, s) in ins {
            m.entry("chr1".to_string()).or_default().entry(p).or_default().push(rcd(s, &format!("{}:{}", s, p)));
        }
        let out = call_peak_candidate_sites_vcf(&m, 50);
        let got: Vec<(i64, Vec<String>)> = out["chr1"].iter().map(|(p, v)| (*p, v.iter().map(|r| r.0.clone()).collect())).collect();
        let exp: Vec<(i64, Vec<String>)> = vec![
            (1, vec!["k".into(), "j".into()]),
            (1008, vec!["c".into(), "d".into(), "a".into(), "f".into(), "b".into(), "e".into()]),
            (5000, vec!["g".into(), "h".into()]),
            (5100, vec!["i".into()]),
        ];
        assert_eq!(got, exp);
    }

    #[test]
    fn negative_one_sentinel_quirk() {
        // a first position of -1 keeps `pre_pos == -1`, so the next position joins the cluster
        // whatever the distance; in set order 1000 comes before -1 (hash -2)
        // (CPython: {1000: [b, a]} for positions -1 and 1000 with window 5)
        let mut m: PopSites = PopSites::new();
        m.entry("c".into()).or_default().entry(1000).or_default().push(rcd("b", "x"));
        m.entry("c".into()).or_default().entry(-1).or_default().push(rcd("a", "y"));
        let out = call_peak_candidate_sites_vcf(&m, 5);
        let got: Vec<(i64, Vec<String>)> = out["c"].iter().map(|(p, v)| (*p, v.iter().map(|r| r.0.clone()).collect())).collect();
        assert_eq!(got, vec![(1000, vec!["b".to_string(), "a".to_string()])]);
    }

    #[test]
    fn allele_record_formatting() {
        let r = AlleleRcd { an: 6, ac: 1, af: 1.0 / 6.0, n_het: 1, n_homo: 0, n_homo_ref: 2, het_freq: 1.0 / 3.0, homo_freq: 0.0, homo_ref_freq: 2.0 / 3.0, n_samples: 3 };
        assert_eq!(
            prepare_one_record("EUR_", &r),
            "EUR_AN=6;EUR_AC=1;EUR_AF=0.16666666666666666;EUR_N_HET=1;EUR_N_HOMALT=0;EUR_N_HOMREF=2;\
             EUR_FREQ_HET=0.3333333333333333;EUR_FREQ_HOMALT=0.0;EUR_FREQ_HOMREF=0.6666666666666666;EUR_N_BI_GENOS=3;"
        );
    }

    #[test]
    fn conflicting_sites_keep_both_unless_replaced() {
        let mut sites = SampleSites::new();
        sites.entry("chr2".into()).or_default().set(140183, "L1 two_side_tprt_both".into());
        sites.entry("chr2".into()).or_default().set(5000, "L1 one_side".into());
        let dir = std::env::temp_dir().join(format!("xtea_vcf_merger_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let f = dir.join("sva.vcf");
        std::fs::write(&f, "chr2\t140195\tSVA two_side_tprt_both\nchr2\t5010\tSVA two_side_tprt_both\nchr2\t5049\tSVA one_side\n").unwrap();
        append_sites_from_vcf(f.to_str().unwrap(), 50, &mut sites).unwrap();
        let got: Vec<(i64, String)> = sites["chr2"].map.iter().map(|(p, r)| (*p, r.clone())).collect();
        assert_eq!(
            got,
            vec![
                (140183, "L1 two_side_tprt_both".to_string()),
                (5000, String::new()),
                (140195, "chr2\t140195\tSVA two_side_tprt_both".to_string()),
                (5010, "chr2\t5010\tSVA two_side_tprt_both".to_string()),
                (5049, "chr2\t5049\tSVA one_side".to_string()),
            ]
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn td_src_most_frequent_first_on_ties() {
        let mut m = TdSrc::new();
        m.entry("EUR".into()).or_default().insert("chr1:1-2".into(), 1);
        m.entry("EUR".into()).or_default().insert("chr2:3-4".into(), 1);
        m.entry("AFR".into()).or_default();
        m.entry("EAS".into()).or_default().insert("".into(), 2);
        m.entry("EAS".into()).or_default().insert("chr9:1".into(), 2);
        m.entry("AMR".into()).or_default().insert("x".into(), 1);
        m.entry("AMR".into()).or_default().insert("y".into(), 3);
        assert_eq!(cvt_td_src_rcd_to_str(&m), "\tEUR:chr1:1-2\tAMR:y");
    }
}
