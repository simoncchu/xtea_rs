//! The individual steps of `x_TEA_main.py` (one function per mode).

use crate::cli::{automatic_gnrt_parameters, automatic_set_molecule_cutoff_for_10x_bam, py_f, set_globals_from_basic, Opts};
use crate::clip_disc_filter::XClipDiscFilter;
use crate::globals as gv;
use crate::intermediate_sites as xsites;
use crate::py::fmt::py_int;
use crate::orphan_transduction::XOrphanTransduction;
use crate::tei_locator::{TeLocator, TeMultiLocator};
use crate::transduction::XTransduction;
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

/// `--transduction`
pub fn run_transduction(opts: &Opts, b_automatic: bool, b_tumor: bool, f_purity: f64, b_resume: bool, n_jobs: usize) -> Result<()> {
    let sf_bam_list = req(&opts.bam, "-b")?;
    let sf_candidate_list = opts.input.as_str();
    let sf_raw_disc = opts.input2.as_str();
    let mut iextnd: i64 = 400;
    let bin_size: i64 = 50000000;
    let sf_cns = opt_str(&opts.reference);
    let s_working_folder = opts.wfolder.as_str();
    println!("Current working folder is: {}\n", s_working_folder);
    let sf_reference = opt_str(&opts.r#ref);
    let sf_flank = opt_str(&opts.fflank);
    let i_flank_lenth = opts.flklen.unwrap_or(0);
    let sf_output = req(&opts.output, "-o")?;
    let sf_rmsk = opt_str(&opts.annotation);
    let i_rep_type = opts.rep_type;
    if b_resume && utils::is_file(sf_output) {
        return Ok(());
    }
    if utils::is_file(&sf_flank) {
        let (rcd, basic_rcd) = automatic_gnrt_parameters(sf_bam_list, &sf_reference, s_working_folder, n_jobs, false, b_tumor, f_purity)?;
        let z = basic_rcd.zero;
        println!("Mean insert size is: {}\n", py_f(basic_rcd.mean_is, z));
        println!("Standard derivation is: {}\n", py_f(basic_rcd.std_var, z));
        let max_is = py_int(basic_rcd.mean_is + 3.0 * basic_rcd.std_var);
        if iextnd < max_is {
            iextnd = max_is;
        }
        let mut i_concord_dist: i64 = 550;
        let f_concord_ratio = 0.25;
        if i_concord_dist < max_is {
            i_concord_dist = max_is;
        }
        set_globals_from_basic(&basic_rcd, max_is);
        let ave_cov = basic_rcd.cov;
        let mut n_clip_cutoff = opts.cliprep;
        let mut n_disc_cutoff = opts.ndisc;
        if b_automatic {
            n_clip_cutoff = rcd.0;
            n_disc_cutoff = rcd.1;
        }
        let xtransduction = XTransduction::new(s_working_folder, n_jobs, &sf_reference);
        let sf_tmp_slct = format!("{}.slct", sf_raw_disc);
        let xannotation = xtransduction.prep_annotation_interval_tree(&sf_rmsk, 225)?;
        xtransduction.re_slct_with_clip_raw_disc_sites(sf_raw_disc, sf_candidate_list, n_disc_cutoff, &xannotation, i_rep_type, b_tumor, &sf_tmp_slct, true)?;
        let x_cd_filter = XClipDiscFilter::new(sf_bam_list, s_working_folder, n_jobs, &sf_reference);
        let i_max_cov = ave_cov * (gv::MAX_COV_TIMES + 1) as f64;
        let sf_output_tmp = format!("{}{}", sf_output, gv::TD_NON_SIBLING_SUFFIX);
        xtransduction.call_candidate_transduction_v3(
            &sf_tmp_slct,
            sf_candidate_list,
            &x_cd_filter,
            &sf_flank,
            &sf_cns,
            i_flank_lenth,
            iextnd,
            bin_size,
            n_clip_cutoff,
            n_disc_cutoff,
            i_concord_dist,
            f_concord_ratio,
            &xannotation,
            sf_bam_list,
            i_rep_type,
            i_max_cov,
            ave_cov,
            &sf_output_tmp,
        )?;
        let xorphan = XOrphanTransduction::new(s_working_folder, n_jobs, &sf_reference);
        let n_half_disc_cutoff = n_disc_cutoff as f64 / 2.0;
        let i_search_win = 2000;
        let sf_updated_cns = sf_output;
        let sf_sibling_td = format!("{}.sibling_transduction_from_existing_list", sf_output);
        xorphan.call_sibling_td_from_existing_list(
            &sf_output_tmp,
            sf_bam_list,
            iextnd,
            n_half_disc_cutoff,
            i_search_win,
            &xannotation,
            i_rep_type,
            i_max_cov,
            sf_updated_cns,
            &sf_sibling_td,
        )?;
        let sf_ori_hc = format!("{}{}", sf_candidate_list, gv::HIGH_CONFIDENT_SUFFIX);
        let sf_new_hc = format!("{}{}", sf_output, gv::HIGH_CONFIDENT_SUFFIX);
        xorphan.td.update_high_confident_callset(&sf_ori_hc, sf_updated_cns, &sf_new_hc)?;
    } else {
        utils::copyfile(sf_candidate_list, sf_output)?;
        let sf_ori_hc = format!("{}{}", sf_candidate_list, gv::HIGH_CONFIDENT_SUFFIX);
        let sf_new_hc = format!("{}{}", sf_output, gv::HIGH_CONFIDENT_SUFFIX);
        utils::copyfile(&sf_ori_hc, &sf_new_hc)?;
    }
    Ok(())
}

/// The common prologue of the `--sibling` branches.
fn sibling_params(opts: &Opts, b_automatic: bool, b_tumor: bool, f_purity: f64, n_jobs: usize) -> Result<(i64, i64, i64, f64)> {
    let sf_bam_list = req(&opts.bam, "-b")?;
    let sf_reference = opt_str(&opts.r#ref);
    let (rcd, basic_rcd) = automatic_gnrt_parameters(sf_bam_list, &sf_reference, &opts.wfolder, n_jobs, false, b_tumor, f_purity)?;
    let z = basic_rcd.zero;
    println!("Mean insert size is: {}\n", py_f(basic_rcd.mean_is, z));
    println!("Standard derivation is: {}\n", py_f(basic_rcd.std_var, z));
    let max_is = py_int(basic_rcd.mean_is + 3.0 * basic_rcd.std_var);
    let mut i_concord_dist: i64 = 550;
    if i_concord_dist < max_is {
        i_concord_dist = max_is;
    }
    gv::set_read_length(basic_rcd.rlth);
    gv::set_insert_size(basic_rcd.mean_is);
    gv::set_average_cov(basic_rcd.cov);
    gv::set_average_cov_is_int(matches!(basic_rcd.cov_py, crate::num::PyNum::Int(_)));
    let mut n_clip_cutoff = opts.cliprep;
    let mut n_disc_cutoff = opts.ndisc;
    if b_automatic {
        n_clip_cutoff = rcd.0;
        n_disc_cutoff = rcd.1;
    }
    Ok((i_concord_dist, n_clip_cutoff, n_disc_cutoff, basic_rcd.mean_is))
}

/// `--sibling`
pub fn run_sibling(opts: &Opts, b_automatic: bool, b_tumor: bool, f_purity: f64, n_jobs: usize) -> Result<()> {
    let sf_bam_list = req(&opts.bam, "-b")?;
    let sf_pre_step_out = opts.input.as_str();
    let sf_raw_disc = opts.input2.as_str();
    let s_working_folder = opts.wfolder.as_str();
    println!("Current working folder is: {}\n", s_working_folder);
    let sf_reference = opt_str(&opts.r#ref);
    let sf_flank = opt_str(&opts.fflank);
    let sf_output = req(&opts.output, "-o")?;
    let sf_rmsk = opt_str(&opts.annotation);
    let i_rep_type = opts.rep_type;
    let sf_black_list = opts.blacklist.as_str();
    if utils::is_file(&sf_flank) {
        let (i_concord_dist, n_clip_cutoff, n_disc_cutoff, mean_is) = sibling_params(opts, b_automatic, b_tumor, f_purity, n_jobs)?;
        let mut xorphan = XOrphanTransduction::new(s_working_folder, n_jobs, &sf_reference);
        xorphan.td.set_boundary_extend(mean_is);
        let xannotation = xorphan.td.prep_annotation_interval_tree(&sf_rmsk, 225)?;
        let sf_tmp_slct2 = format!("{}.slct2", sf_raw_disc);
        let sf_output_tmp = format!("{}{}", sf_pre_step_out, gv::TD_NON_SIBLING_SUFFIX);
        xorphan.td.re_slct_with_clip_raw_disc_sites(sf_raw_disc, &sf_output_tmp, n_disc_cutoff, &xannotation, i_rep_type, b_tumor, &sf_tmp_slct2, false)?;
        let m_failed_ori_td = xorphan.distinguish_source_from_insertion_for_td(sf_pre_step_out, sf_bam_list, i_concord_dist, sf_black_list, &sf_rmsk)?;
        let sf_sibling_td = sf_output;
        if !utils::is_file(sf_black_list) {
            println!("Blacklist file {} does not exist!", sf_black_list);
        }
        xorphan.call_novel_sibling_td_from_raw_list(&sf_tmp_slct2, sf_bam_list, i_concord_dist, n_clip_cutoff, n_disc_cutoff, sf_black_list, &sf_rmsk, sf_sibling_td)?;
        XOrphanTransduction::append_to_existing_list(sf_sibling_td, sf_pre_step_out, &m_failed_ori_td)?;
        XOrphanTransduction::append_to_existing_list(sf_sibling_td, &format!("{}{}", sf_pre_step_out, gv::HIGH_CONFIDENT_SUFFIX), &m_failed_ori_td)?;
    } else if utils::is_file(sf_pre_step_out) {
        let (i_concord_dist, _n_clip_cutoff, _n_disc_cutoff, _mean_is) = sibling_params(opts, b_automatic, b_tumor, f_purity, n_jobs)?;
        let xorphan = XOrphanTransduction::new(s_working_folder, n_jobs, &sf_reference);
        let m_failed_ori_td = xorphan.distinguish_source_from_insertion_for_td(sf_pre_step_out, sf_bam_list, i_concord_dist, sf_black_list, &sf_rmsk)?;
        XOrphanTransduction::update_existing_list_only(sf_pre_step_out, &m_failed_ori_td)?;
    }
    Ok(())
}
/// `--case_control`: somatic calls of the case filtered with the control BAMs
/// (`--somatic_hc`: keep the high-confidence calls that are in the raw somatic set).
pub fn run_case_control(opts: &Opts, b_automatic: bool, _b_tumor: bool, n_jobs: usize) -> Result<()> {
    let sf_candidate_list = opts.input.as_str();
    let s_working_folder = opts.wfolder.as_str();
    if !opts.somatic_hc {
        let sf_bam_list = req(&opts.bam, "-b")?;
        let sf_ref = opt_str(&opts.r#ref);
        let mut nclip_cutoff = opts.cliprep;
        let mut ndisc_cutoff = opts.ndisc;
        let sf_rep_cns = req(&opts.reference, "-r")?;
        let sf_flank = req(&opts.fflank, "--fflank")?;
        let b_force = true;
        let mut n_polya_cutoff = 0i64;
        let mut basic_rcd = None;
        if b_automatic {
            let (rcd, brcd) = crate::cli::automatic_gnrt_parameters_case_control(sf_bam_list, &sf_ref, s_working_folder, n_jobs, b_force)?;
            nclip_cutoff = rcd.0;
            ndisc_cutoff = rcd.1;
            n_polya_cutoff = rcd.2;
            basic_rcd = Some(brcd);
        }
        let ccm = crate::somatic::CaseControlMode::new(&sf_ref, s_working_folder, n_jobs);
        // without the automatic parameters the original subscripts None here
        let basic = basic_rcd.ok_or_else(|| anyhow!("TypeError: 'NoneType' object is not subscriptable"))?;
        let rlth = basic.rlth;
        let mean_is = basic.mean_is;
        gv::set_insert_size(py_int(mean_is) as f64);
        let std_var = basic.std_var;
        let max_is = py_int(mean_is + 3.0 * std_var) + py_int(rlth);
        let extnd = max_is;
        let bin_size = 50000000;
        println!("clip,disc,polyA-cutoff is ({}, {}, {})", nclip_cutoff, ndisc_cutoff, n_polya_cutoff);
        n_polya_cutoff = ndisc_cutoff;
        let sf_output = req(&opts.output, "-o")?;
        let i_flk_len = opts.flklen.ok_or_else(|| anyhow!("option --flklen is required for this step"))?;
        ccm.call_somatic_te_insertion(
            sf_bam_list,
            sf_candidate_list,
            extnd,
            nclip_cutoff,
            ndisc_cutoff,
            n_polya_cutoff,
            sf_rep_cns,
            sf_flank,
            i_flk_len,
            bin_size,
            sf_output,
        )
    } else {
        let sf_output = req(&opts.output, "-o")?;
        let ccm = crate::somatic::CaseControlMode::new("", s_working_folder, 1);
        ccm.parse_high_confident_somatic(sf_candidate_list, &opts.input2, sf_output)
    }
}
/// `--postF`: post filtering of the raw calls.
pub fn run_post_filter(opts: &Opts, b_tumor: bool, n_jobs: usize) -> Result<()> {
    let s_working_folder = opts.wfolder.as_str();
    let sf_xtea_rslt = opts.input.as_str();
    let sf_rmsk = req(&opts.annotation, "-a")?;
    let sf_new_out = req(&opts.output, "-o")?;
    let i_rep_type = opts.rep_type;
    let sf_black_list = opts.blacklist.as_str();
    let i_min_copy_len = 225;
    if opts.post_f_mosaic {
        let xpf_mosic = crate::mosaic::MosaicCaller::new(s_working_folder, n_jobs);
        return xpf_mosic.run_call_mosaic(sf_xtea_rslt, sf_rmsk, i_min_copy_len, i_rep_type, sf_black_list, sf_new_out);
    }
    let xbasic = crate::basic_info::XBasicInfo::new(s_working_folder, 1, "");
    let (f_cov, _f_ave_cov) = xbasic.load_get_cov()?;
    let xpost_filter = crate::post_filter::XPostFilter::new(s_working_folder, n_jobs);
    xpost_filter.run_post_filtering(sf_xtea_rslt, sf_rmsk, i_min_copy_len, i_rep_type, f_cov, sf_black_list, sf_new_out, b_tumor)
}

/// `--gntp_classify`: genotype the calls with the DF21 deep-forest model.
pub fn run_gntp_classify(opts: &Opts) -> Result<()> {
    if opts.train_gntp {
        return Err(anyhow!(
            "--train_gntp is not supported by xtea-rs: training used site lists hard-coded to the original authors' cluster and scikit-learn; train with the Python version"
        ));
    }
    let sf_xtea = opts.input.as_str();
    let sf_new = req(&opts.output, "-o")?;
    let mut gc = crate::genotype_classify::GntpClassifierDf21::new();
    gc.predict_for_site(&opts.model, sf_xtea, sf_new)
}

/// `--gVCF`: convert the final calls to VCF.
pub fn run_gvcf(opts: &Opts) -> Result<()> {
    run_gvcf_with_sample_id(opts, None)
}

/// `--gVCF`; the cloud entry point takes the sample id from `--id` instead of the BAM name.
pub fn run_gvcf_with_sample_id(opts: &Opts, sample_id: Option<&str>) -> Result<()> {
    let sf_ref = opts.r#ref.clone().unwrap_or_else(|| "None".to_string());
    let sf_bam_list = req(&opts.bam, "-b")?;
    let sf_prefix = req(&opts.output, "-o")?;
    crate::gvcf::run_gvcf_step(&sf_ref, sf_bam_list, &opts.input, opts.rep_type, sf_prefix, sample_id)
}
/// `--joint` (with `--postFmosaic`): keep the mosaic calls seen in a single sample.
pub fn run_joint(opts: &Opts, n_jobs: usize) -> Result<()> {
    let islack = 50;
    if opts.post_f_mosaic {
        let mzc_joint = crate::mosaic::MosaicJointCalling::new(&opts.wfolder, n_jobs)?;
        let sf_out = req(&opts.output, "-o")?;
        mzc_joint.call_mosaic_from_multi_samples(&opts.input, islack, sf_out)?;
    }
    Ok(())
}
/// `--igv`: IGV batch script for screenshots of the given sites.
pub fn run_igv(opts: &Opts, _n_jobs: usize) -> Result<()> {
    let mut sf_sites = opts.input.clone();
    let sf_bam_list = req(&opts.bam, "-b")?;
    let sf_out = req(&opts.output, "-o")?;
    let sf_gnm = opts.r#ref.clone().unwrap_or_else(|| "None".to_string());
    let x_igv = crate::igv::XIgv::new();
    if opts.single_sample {
        sf_sites = x_igv.gnrt_sites_single_sample(&sf_sites, sf_bam_list)?;
    }
    x_igv.prepare_igv_scripts_multi_bams(&sf_sites, sf_bam_list, &opts.wfolder, opts.extend, &sf_gnm, sf_out)
}
pub fn run_collect_10x(_opts: &Opts, _n_jobs: usize) -> Result<()> {
    not_ported("-E")
}
pub fn run_mutation(_opts: &Opts, _n_jobs: usize) -> Result<()> {
    not_ported("-I")
}
/// `--gene`: annotate the calls with the genes of a GFF3 file.
pub fn run_gene(opts: &Opts) -> Result<()> {
    run_gene_variant(opts, true)
}

/// `--gene`; the cloud entry point ignores `-e` (`b_user_extend` false).
pub fn run_gene_variant(opts: &Opts, b_user_extend: bool) -> Result<()> {
    let sf_gene_annotation = req(&opts.annotation, "-a")?;
    let sf_output = req(&opts.output, "-o")?;
    let mut gff = crate::gene_annotation::Gff3::new(sf_gene_annotation);
    let mut iextnd = gv::UP_DOWN_GENE;
    if b_user_extend && opts.extend > iextnd {
        iextnd = opts.extend;
    }
    gff.load_gene_annotation_with_extnd(iextnd)?;
    gff.index_gene_annotation_interval_tree();
    gff.annotate_results(&opts.input, sf_output)
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
