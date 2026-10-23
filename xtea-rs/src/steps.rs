//! The individual steps of `x_TEA_main.py` (one function per mode).

use crate::cli::{automatic_gnrt_parameters, automatic_set_molecule_cutoff_for_10x_bam, py_f, set_globals_from_basic, Opts};
use crate::clip_disc_filter::XClipDiscFilter;
use crate::globals as gv;
use crate::intermediate_sites as xsites;
use crate::py::fmt::py_int;
use crate::tei_locator::{TeLocator, TeMultiLocator};
use crate::utils;
use anyhow::{anyhow, Result};

pub(crate) fn req<'a>(o: &'a Option<String>, name: &str) -> Result<&'a str> {
    o.as_deref().ok_or_else(|| anyhow!("option {} is required for this step", name))
}

pub(crate) fn opt_str(o: &Option<String>) -> String {
    o.clone().unwrap_or_default()
}

/// `-C`: candidate sites from clipped reads
pub fn run_clip(opts: &Opts, b_automatic: bool, b_tumor: bool, f_purity: f64, b_resume: bool, n_jobs: usize) -> Result<()> {
    println!("Working on \"clip\" step!");
    let sf_bam_list = opts.input.as_str();
    let s_working_folder = opts.wfolder.as_str();
    let sf_rep_cns = opt_str(&opts.cns);
    let sf_rep = opt_str(&opts.reference);
    let sf_annotation = opt_str(&opts.annotation);
    let sf_out = req(&opts.output, "-o")?;
    let b_se = opts.single;
    let sf_ref = opt_str(&opts.r#ref);
    let b_force = opts.force;
    let b_mosaic = opts.mosaic;
    if b_force {
        gv::set_force_clean();
    }
    gv::set_initial_min_clip_cutoff(opts.siteclip);
    let mut cutoff_left_clip = opts.lclip;
    let mut cutoff_right_clip = opts.rclip;
    let mut cutoff_clip_mate_in_rep = opts.cliprep;
    if !b_resume || !utils::is_file(sf_out) {
        if b_automatic {
            let (rcd, _basic) = automatic_gnrt_parameters(sf_bam_list, &sf_ref, s_working_folder, n_jobs, b_force, b_tumor, f_purity)?;
            cutoff_left_clip = rcd.0;
            cutoff_right_clip = rcd.0;
            cutoff_clip_mate_in_rep = rcd.2;
        }
        println!("Clip cutoff: {}, {}, {} are used!!!", cutoff_left_clip, cutoff_right_clip, cutoff_clip_mate_in_rep);
        let tem_locator = TeMultiLocator::new(sf_bam_list, s_working_folder, n_jobs, &sf_ref);
        let max_cov_cutoff = py_int(15.0 * opts.cov);
        let wfolder_pub_clip = opt_str(&opts.cwfolder);
        tem_locator.call_tei_candidate_sites_from_multiple_alignmts(
            &sf_annotation,
            &sf_rep_cns,
            &sf_rep,
            b_se,
            cutoff_left_clip,
            cutoff_right_clip,
            cutoff_clip_mate_in_rep,
            b_mosaic,
            &wfolder_pub_clip,
            b_force,
            max_cov_cutoff,
            sf_out,
        )?;
    }
    Ok(())
}

/// `-D`: filter the candidates with discordant pairs
pub fn run_discordant(opts: &Opts, b_automatic: bool, b_tumor: bool, f_purity: f64, b_resume: bool, n_jobs: usize) -> Result<()> {
    println!("Working on \"disc\" step!");
    let sf_bam_list = req(&opts.bam, "-b")?;
    let s_working_folder = opts.wfolder.as_str();
    let sf_annotation = opt_str(&opts.annotation);
    let sf_candidate_list = opts.input.as_str();
    let sf_out = req(&opts.output, "-o")?;
    let sf_ref = opt_str(&opts.r#ref);
    let mut peak_window = 100;
    if opts.post_f_mosaic || opts.somatic {
        peak_window = 30;
    }
    if !b_resume || !utils::is_file(sf_out) {
        let m_original_sites = xsites::load_in_candidate_list(sf_candidate_list)?;
        let sf_peak_sites = format!("{}clip_peak_candidate.list", s_working_folder);
        let m_sites_clip_peak = xsites::call_peak_candidate_sites_with_std_derivation(&m_original_sites, peak_window);
        xsites::output_candidate_sites(&m_sites_clip_peak, &sf_peak_sites)?;
        drop(m_original_sites);
        let (rcd, basic_rcd) = automatic_gnrt_parameters(sf_bam_list, &sf_ref, s_working_folder, n_jobs, true, b_tumor, f_purity)?;
        let rlth = basic_rcd.rlth;
        let mean_is = basic_rcd.mean_is;
        let std_var = basic_rcd.std_var;
        let max_is = py_int(mean_is + 3.0 * std_var) + py_int(rlth);
        let iextend = max_is;
        let i_is = 100000;
        let f_dev = std_var;
        let mut n_disc_cutoff = opts.ndisc;
        if b_automatic {
            n_disc_cutoff = rcd.1;
        }
        println!("Discordant cutoff: {} is used!!!", n_disc_cutoff);
        let sf_tmp = format!("{}disc_tmp.list", s_working_folder);
        let sf_raw_disc = format!("{}{}", sf_out, gv::RAW_DISC_TMP_SUFFIX);
        let tem_locator = TeMultiLocator::new(sf_bam_list, s_working_folder, n_jobs, &sf_ref);
        tem_locator.filter_candidate_sites_by_discordant_pairs_multi_alignmts(
            &m_sites_clip_peak,
            iextend,
            i_is,
            f_dev,
            n_disc_cutoff,
            &sf_annotation,
            &sf_tmp,
            &sf_raw_disc,
            b_tumor,
        )?;
        xsites::merge_clip_disc(&sf_tmp, sf_candidate_list, sf_out)?;
    }
    Ok(())
}

/// `-N`: filter with the reads re-aligned to the consensus
pub fn run_filter_csn(opts: &Opts, b_automatic: bool, b_tumor: bool, f_purity: f64, b_resume: bool, n_jobs: usize) -> Result<()> {
    println!("Working on \"clip-disc-filtering\" step!");
    let sf_bam_list = req(&opts.bam, "-b")?;
    let s_working_folder = opts.wfolder.as_str();
    println!("Current working folder is: {}\n", s_working_folder);
    let sf_ref = opt_str(&opts.r#ref);
    let sf_candidate_list = opts.input.as_str();
    let sf_raw_disc = opts.input2.as_str();
    let mut iextnd: i64 = 400;
    let bin_size: i64 = 50000000;
    let sf_cns = req(&opts.reference, "-r")?;
    let bmapped_cutoff = gv::MIN_CLIP_MAPPED_RATIO;
    let mut i_concord_dist: i64 = 550;
    let f_concord_ratio = 0.45;
    let sf_output = req(&opts.output, "-o")?;
    let sf_flank = opt_str(&opts.fflank);
    let i_flank_lenth = opts.flklen.unwrap_or(0);
    if !b_resume || !utils::is_file(sf_output) {
        let (rcd, basic_rcd) = automatic_gnrt_parameters(sf_bam_list, &sf_ref, s_working_folder, n_jobs, true, b_tumor, f_purity)?;
        let z = basic_rcd.zero;
        println!("Mean insert size is: {}\n", py_f(basic_rcd.mean_is, z));
        println!("Standard derivation is: {}\n", py_f(basic_rcd.std_var, z));
        let max_is = py_int(basic_rcd.mean_is + 3.0 * basic_rcd.std_var);
        if iextnd < max_is {
            iextnd = max_is;
        }
        if i_concord_dist < max_is {
            i_concord_dist = max_is;
        }
        set_globals_from_basic(&basic_rcd, max_is);
        println!("Read length is: {}\n", py_f(basic_rcd.rlth, z));
        println!("Maximum insert size is: {}\n", max_is);
        println!("Average coverage is: {}\n", basic_rcd.cov_py);
        let mut n_clip_cutoff = opts.cliprep;
        let mut n_disc_cutoff = opts.ndisc;
        if b_automatic {
            n_clip_cutoff = rcd.0;
            n_disc_cutoff = rcd.1;
        }
        println!("Filter (on cns) cutoff: {} and {} are used!!!", n_clip_cutoff, n_disc_cutoff);
        let x_cd_filter = XClipDiscFilter::new(sf_bam_list, s_working_folder, n_jobs, &sf_ref);
        x_cd_filter.call_meis_consensus(
            sf_candidate_list,
            sf_raw_disc,
            iextnd,
            bin_size,
            sf_cns,
            &sf_flank,
            i_flank_lenth,
            bmapped_cutoff,
            i_concord_dist,
            f_concord_ratio,
            n_clip_cutoff,
            n_disc_cutoff,
            sf_output,
        )?;
    }
    Ok(())
}

/// `--mosaic`: candidate sites from clipped reads (mosaic version, checks polyA)
pub fn run_mosaic_clip(opts: &Opts, b_automatic: bool, n_jobs: usize) -> Result<()> {
    println!("Working on mosaic \"clip\" step!");
    let sf_bam_list = opts.input.as_str();
    let s_working_folder = opts.wfolder.as_str();
    let sf_rep_cns = opt_str(&opts.cns);
    let sf_rep = opt_str(&opts.reference);
    let sf_annotation = opt_str(&opts.annotation);
    let sf_out = req(&opts.output, "-o")?;
    let b_se = opts.single;
    let sf_ref = opt_str(&opts.r#ref);
    let b_force = opts.force;
    if b_force {
        gv::set_force_clean();
    }
    gv::set_initial_min_clip_cutoff(opts.siteclip);
    let mut cutoff_left_clip = opts.lclip;
    let mut cutoff_right_clip = opts.rclip;
    let mut cutoff_clip_mate_in_rep = opts.cliprep;
    let cutoff_polya = 1;
    if b_automatic {
        let (rcd, _) = automatic_gnrt_parameters(sf_bam_list, &sf_ref, s_working_folder, n_jobs, false, false, 0.45)?;
        cutoff_left_clip = rcd.0;
        cutoff_right_clip = rcd.0;
        cutoff_clip_mate_in_rep = rcd.2;
    }
    println!("Clip cutoff: {}, {}, {} are used!!!", cutoff_left_clip, cutoff_right_clip, cutoff_clip_mate_in_rep);
    let tem_locator = TeMultiLocator::new(sf_bam_list, s_working_folder, n_jobs, &sf_ref);
    let max_cov_cutoff = py_int(0.5 * opts.cov);
    let wfolder_pub_clip = opt_str(&opts.cwfolder);
    tem_locator.call_tei_candidate_sites_from_multiple_alignmts_mosaic(
        &sf_annotation,
        &sf_rep_cns,
        &sf_rep,
        b_se,
        cutoff_left_clip,
        cutoff_right_clip,
        cutoff_clip_mate_in_rep,
        cutoff_polya,
        &wfolder_pub_clip,
        b_force,
        max_cov_cutoff,
        sf_out,
    )
}

/// `-B`: filter by barcode (molecule) coverage, 10X data
pub fn run_barcode(opts: &Opts, n_jobs: usize) -> Result<()> {
    let sf_ori_bam = req(&opts.bam, "-b")?;
    let sf_barcode_bam = opt_str(&opts.barcode_bam);
    let mut s_working_folder = opts.wfolder.clone();
    let sf_candidate_list = opts.input.as_str();
    let sf_out = req(&opts.output, "-o")?;
    let sf_ref = opt_str(&opts.r#ref);
    let m_input_sites = xsites::load_in_candidate_list(sf_candidate_list)?;
    let iextend = 2500;
    let i_cov_cutoff = automatic_set_molecule_cutoff_for_10x_bam(sf_ori_bam, &sf_ref, &s_working_folder, n_jobs)?;
    println!("Using molecule coverage cutoff: {}", i_cov_cutoff);
    if !s_working_folder.ends_with('/') {
        s_working_folder.push('/');
    }
    let sf_tmp = format!("{}barcode_tmp.list", s_working_folder);
    let caller = TeLocator::new(sf_ori_bam, &sf_barcode_bam, &s_working_folder, n_jobs, &sf_ref);
    let m_sites_barcode = caller.filter_candidate_sites_by_barcode_coverage(&m_input_sites, iextend, i_cov_cutoff)?;
    caller.output_candidate_sites(&m_sites_barcode, &sf_tmp)?;
    let sf_tmp2 = format!("{}barcode_tmp2.list", s_working_folder);
    xsites::merge_clip_disc_barcode(&sf_tmp, sf_candidate_list, &sf_tmp2)?;
    xsites::combine_closing_sites(&sf_tmp2, 200, sf_out)
}

fn not_ported(step: &str) -> Result<()> {
    Err(anyhow!("step `{}` is not available in this build of xtea-rs yet", step))
}

pub fn run_transduction(_opts: &Opts, _b_automatic: bool, _b_tumor: bool, _f_purity: f64, _b_resume: bool, _n_jobs: usize) -> Result<()> {
    not_ported("--transduction")
}
pub fn run_sibling(_opts: &Opts, _b_automatic: bool, _b_tumor: bool, _f_purity: f64, _n_jobs: usize) -> Result<()> {
    not_ported("--sibling")
}
pub fn run_case_control(_opts: &Opts, _b_automatic: bool, _b_tumor: bool, _n_jobs: usize) -> Result<()> {
    not_ported("--case_control")
}
pub fn run_post_filter(_opts: &Opts, _b_tumor: bool, _n_jobs: usize) -> Result<()> {
    not_ported("--postF")
}
pub fn run_gntp_classify(_opts: &Opts) -> Result<()> {
    not_ported("--gntp_classify")
}
pub fn run_gvcf(_opts: &Opts) -> Result<()> {
    not_ported("--gVCF")
}
pub fn run_joint(_opts: &Opts, _n_jobs: usize) -> Result<()> {
    not_ported("--joint")
}
pub fn run_igv(_opts: &Opts, _n_jobs: usize) -> Result<()> {
    not_ported("--igv")
}
pub fn run_collect_10x(_opts: &Opts, _n_jobs: usize) -> Result<()> {
    not_ported("-E")
}
pub fn run_mutation(_opts: &Opts, _n_jobs: usize) -> Result<()> {
    not_ported("-I")
}
pub fn run_gene(_opts: &Opts) -> Result<()> {
    not_ported("--gene")
}
pub fn run_assembly(_opts: &Opts, _n_jobs: usize) -> Result<()> {
    not_ported("-A")
}
pub fn run_map(_opts: &Opts, _n_jobs: usize) -> Result<()> {
    not_ported("-M")
}
pub fn run_flk_map(_opts: &Opts, _n_jobs: usize) -> Result<()> {
    not_ported("--flk_map")
}
pub fn run_filter_asm(_opts: &Opts, _n_jobs: usize) -> Result<()> {
    not_ported("-F")
}
