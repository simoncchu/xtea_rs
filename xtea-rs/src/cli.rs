//! Command line of `xtea_main` (port of `x_TEA_main.py`).  Option names, defaults and the
//! precedence of the modes (the `if/elif` chain of the original) are kept.

use crate::annotation::XAnnotation;
use crate::basic_info::{BasicRcd, XBasicInfo};
use crate::globals as gv;
use crate::num::PyNum;
use crate::parameter;
use anyhow::{anyhow, Result};
use clap::Parser;

#[derive(Parser, Debug, Clone, Default)]
#[command(
    name = "xtea_main",
    about = "xTea: transposable element insertion caller (Rust re-implementation of x_TEA_main.py)",
    infer_long_args = true,
    args_override_self = true,
    disable_help_flag = true
)]
pub struct Opts {
    #[arg(short = 'h', long = "help", help = "show this help message and exit")]
    pub help: bool,
    #[arg(short = 'P', long = "preprocess", help = "Preprocessing stpes")]
    pub preprocess: bool,
    #[arg(short = 'Q', long = "collectclip", help = "Call clipped reads from alignment")]
    pub collect_clip: bool,
    #[arg(short = 'C', long = "clip", help = "Call candidate TEI sites from clipped reads")]
    pub clip: bool,
    #[arg(short = 'S', long = "single", help = "Call clip positions from single-end reads")]
    pub single: bool,
    #[arg(short = 'D', long = "discordant", help = "Filter with discordant paired end reads")]
    pub discordant: bool,
    #[arg(short = 'N', long = "filter_csn", help = "Filter out candidate sites from map position on consensus")]
    pub filter_csn: bool,
    #[arg(long = "resume", help = "Resume the running, which will skip the step if output file already exists!")]
    pub resume: bool,
    #[arg(long = "transduction", help = "Call transduction for sites")]
    pub transduction: bool,
    #[arg(long = "sibling", help = "Call sibling transduction for sites")]
    pub sibling: bool,
    #[arg(long = "spectrum", help = "Spectrum analysis by tumor type")]
    pub spectrum: bool,
    #[arg(short = 'B', long = "barcode", help = "Indicate the input is 10X bam [10X; not supported]")]
    pub barcode: bool,
    #[arg(short = 'E', long = "collect", help = "Collect reads for candidate sites [10X; not supported]")]
    pub collect: bool,
    #[arg(short = 'I', long = "mutation", help = "Call internal mutation [10X; not supported]")]
    pub mutation: bool,
    #[arg(short = 'U', long = "collect_Illumina", help = "Collect reads for candidate sites from normal illumina alignment")]
    pub collect_illumina: bool,
    #[arg(short = 'F', long = "filter_asm", help = "Filter out candidate sites from assembly [10X; not supported]")]
    pub filter_asm: bool,
    #[arg(short = 'G', long = "contig_realign", help = "Filter out candidate sites from assembly")]
    pub contig_realign: bool,
    #[arg(short = 'T', long = "trace", help = "Trace the sources of TEIs")]
    pub trace: bool,
    #[arg(short = 'A', long = "assembly", help = "Do local assembly for collected reads [10X; not supported]")]
    pub assembly: bool,
    #[arg(short = 'L', long = "local", help = "Assemble the TEIs on local machine [10X; unused]")]
    pub local: bool,
    #[arg(short = 'M', long = "map", help = "map flank regions to the assembled contigs [10X; not supported]")]
    pub map: bool,
    #[arg(short = 'V', long = "visualization", help = "Show the heatmap figure of the selected regions")]
    pub visualization: bool,
    #[arg(short = 'K', long = "withflank", help = "Keep the flank regions with the repeat copies")]
    pub withflank: bool,
    #[arg(short = 'J', long = "joint", help = "Joint calling")]
    pub joint: bool,
    #[arg(long = "mit", help = "Indicate call mitochondrion insertion")]
    pub mit: bool,
    #[arg(long = "dna", help = "Not RNA mediated insertion (no polyA)")]
    pub dna: bool,
    #[arg(long = "cbs", help = "check by sample")]
    pub cbs: bool,
    #[arg(long = "sva", help = "For SVA insertion calling")]
    pub sva: bool,
    #[arg(long = "gntp_feature", help = "Collect genotyping features from bam")]
    pub gntp_feature: bool,
    #[arg(long = "postF", help = "Post filtering module")]
    pub post_f: bool,
    #[arg(long = "gntp_classify", help = "Train/predict genotpe classifier")]
    pub gntp_classify: bool,
    #[arg(long = "train_gntp", help = "Train the genotype classifer")]
    pub train_gntp: bool,
    #[arg(long = "postFmosaic", help = "Post filtering module for mosaic events")]
    pub post_f_mosaic: bool,
    #[arg(long = "igv", help = "Prepare screenshot command for given sites")]
    pub igv: bool,
    #[arg(long = "force", help = "Force to start from the very beginning")]
    pub force: bool,
    #[arg(long = "case_control", help = "case control mode")]
    pub case_control: bool,
    #[arg(long = "tumor", help = "Working on tumor samples")]
    pub tumor: bool,
    #[arg(long = "gVCF", help = "Generate the gVCF from xTEA raw output")]
    pub g_vcf: bool,
    #[arg(long = "bed", help = "Input annotation in bed format")]
    pub bed: bool,
    #[arg(long = "mosaic", help = "Call mosaic events")]
    pub mosaic: bool,
    #[arg(long = "flk_map", help = "Map flanks to contigs [10X; not supported]")]
    pub flk_map: bool,
    #[arg(long = "analysis", help = "Result analysis")]
    pub analysis: bool,
    #[arg(allow_hyphen_values = true, long = "flank", help = "flank regions")]
    pub flank: Option<String>,
    #[arg(long = "sv", help = "Call promoted SVs")]
    pub sv: bool,
    #[arg(long = "gene", help = "Check whether the insertion falls in genes")]
    pub gene: bool,
    #[arg(long = "somatic", help = "Only call somatic events from high coverage normal samples")]
    pub somatic: bool,
    #[arg(long = "somatic_hc", help = "Get high confident somatic events")]
    pub somatic_hc: bool,
    #[arg(long = "user", help = "User specific parameters, by default automatically calc the parameters")]
    pub user_specific: bool,
    #[arg(long = "single_sample", help = "For single sample (like igv screenshot)")]
    pub single_sample: bool,

    #[arg(allow_hyphen_values = true, short = 'i', long = "input", value_name = "FILE", default_value = "", help = "input file ")]
    pub input: String,
    #[arg(allow_hyphen_values = true, long = "input2", value_name = "FILE", default_value = "", help = "input file2 ")]
    pub input2: String,
    #[arg(allow_hyphen_values = true, short = 'r', long = "reference", value_name = "FILE", help = "The reference file ")]
    pub reference: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'a', long = "annotation", value_name = "FILE", help = "The annotation file ")]
    pub annotation: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'b', long = "bam", value_name = "FILE", help = "Input bam file")]
    pub bam: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'd', long = "barcode_bam", value_name = "FILE", help = "Input barcode indexed bam file [10X; unused]")]
    pub barcode_bam: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'o', long = "output", value_name = "FILE", help = "The output file")]
    pub output: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'p', long = "path", default_value = "./", help = "Working folder")]
    pub wfolder: String,
    #[arg(allow_hyphen_values = true, long = "cp", help = "Working folder for shared clipped reads")]
    pub cwfolder: Option<String>,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, short = 'n', long = "cores", default_value_t = 1, help = "number of cores")]
    pub cores: i64,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, short = 'e', long = "extend", default_value_t = 0, help = "extend length")]
    pub extend: i64,
    #[arg(allow_hyphen_values = true, short = 'u', long = "dup", value_name = "FILE", help = "duplication files")]
    pub duplication: Option<String>,
    #[arg(allow_hyphen_values = true, long = "fflank", value_name = "FILE", help = "flank region file")]
    pub fflank: Option<String>,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, long = "flklen", help = "flank region file")]
    pub flklen: Option<i64>,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_float, long = "purity", default_value_t = 0.45, help = "Tumor purity")]
    pub purity: f64,
    #[arg(allow_hyphen_values = true, long = "ref", value_name = "FILE", help = "repeat consensus/copies")]
    pub r#ref: Option<String>,
    #[arg(allow_hyphen_values = true, long = "cns", value_name = "FILE", help = "repeat consensus")]
    pub cns: Option<String>,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, long = "sc", default_value_t = 2, help = "cutoff of minimum # of clipped reads at the exact position, use larger value for 10X")]
    pub siteclip: i64,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, long = "lc", default_value_t = 3, help = "cutoff of minimum # of left clipped reads")]
    pub lclip: i64,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, long = "rc", default_value_t = 3, help = "cutoff of minimum # of rightt clipped reads")]
    pub rclip: i64,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, long = "cr", default_value_t = 1, help = "cutoff of minimum # of clipped parts fall in repeats")]
    pub cliprep: i64,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, long = "nd", default_value_t = 5, help = "cutoff of minimum # of discordant pair")]
    pub ndisc: i64,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, long = "nb", default_value_t = 500, help = "cutoff of maximum # of molecure coverage [10X; unused]")]
    pub nbarcode: i64,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, long = "teilen", help = "minimum length of the insertion for future analysis")]
    pub teilen: Option<i64>,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_float, long = "cov", default_value_t = 30.0, help = "approximate read depth")]
    pub cov: f64,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, long = "iniclip", default_value_t = 2, help = "initial minimum clip cutoff")]
    pub iniclip: i64,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_float, long = "af1", default_value_t = 0.005, help = "minimal allel fraction")]
    pub af1: f64,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_float, long = "af2", default_value_t = 0.45, help = "minimal allel fraction")]
    pub af2: f64,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, long = "rmsk_extnd", default_value_t = 100, help = "Length of the left extended region when loading the repeatmasker output")]
    pub rmsk_extnd: i64,
    #[arg(allow_hyphen_values = true, value_parser = crate::py::optparse::parse_int, long = "rtype", default_value_t = 1, help = "type of repeats: 1-L1, 2-Alu, 4-SVA, 8-HERV, 16-MIT, 32-MSTA")]
    pub rep_type: i64,
    #[arg(allow_hyphen_values = true, long = "blacklist", value_name = "FILE", default_value = "null", help = "Reference panel database for filtering, or a blacklist region")]
    pub blacklist: String,
    #[arg(allow_hyphen_values = true, long = "model", value_name = "FILE", default_value = "null", help = "Already trained model (.pkl file) for genotype classification")]
    pub model: String,
    #[arg(allow_hyphen_values = true, long = "seed", help = "Seed of the random site sampling (default: $XTEA_SEED if set, otherwise random as the original)")]
    pub seed: Option<u64>,
    /// Positional arguments are accepted and ignored (optparse collects them separately).
    #[arg(hide = true)]
    pub extra_args: Vec<String>,
}

fn req<'a>(o: &'a Option<String>, name: &str) -> Result<&'a str> {
    o.as_deref().ok_or_else(|| anyhow!("option {} is required for this step", name))
}

/// Python `str()` of a value computed as `x / n` (float) or the int 0 fallback.
pub fn py_f(x: f64, zero: bool) -> String {
    if zero {
        "0".to_string()
    } else {
        crate::py::fmt::float_repr(x)
    }
}

/// `automatic_gnrt_parameters` -> ((nclip, ndisc, nclip_disc), basic record)
pub fn automatic_gnrt_parameters(
    sf_bam_list: &str,
    sf_ref: &str,
    s_working_folder: &str,
    n_jobs: usize,
    b_force: bool,
    b_tumor: bool,
    f_purity: f64,
) -> Result<((i64, i64, i64), BasicRcd)> {
    let search_win = 500;
    let x_basic_info = XBasicInfo::new(s_working_folder, n_jobs, sf_ref);
    let rcd = x_basic_info.get_cov_is_rlth(sf_bam_list, sf_ref, search_win, b_force)?;
    let mut f_cov = rcd.cov_py;
    if b_tumor {
        f_cov = PyNum::Float(f_cov.f() * f_purity);
    }
    let par_rcd = parameter::parameters().get_par_by_cov(f_cov.f());
    println!(
        "Ave coverage is {}: automatic parameters (clip, disc, clip-disc) with value ({}, {} ,{})\n",
        f_cov, par_rcd.0, par_rcd.1, par_rcd.2
    );
    Ok((par_rcd, rcd))
}

/// `automatic_gnrt_parameters_case_control`
pub fn automatic_gnrt_parameters_case_control(sf_bam_list: &str, sf_ref: &str, s_working_folder: &str, n_jobs: usize, b_force: bool) -> Result<((i64, i64, i64), BasicRcd)> {
    let search_win = 500;
    let x_basic_info = XBasicInfo::new(s_working_folder, n_jobs, sf_ref);
    let rcd = x_basic_info.get_cov_is_rlth(sf_bam_list, sf_ref, search_win, b_force)?;
    let par_rcd = parameter::case_control_filter_pars().get_par_by_cov(rcd.cov);
    println!(
        "Ave coverage is {}: automatic parameters (clip, disc, clip-disc) with value ({}, {} ,{})\n",
        rcd.cov_py, par_rcd.0, par_rcd.1, par_rcd.2
    );
    Ok((par_rcd, rcd))
}

/// Set READ_LENGTH / insert size / AVE_COVERAGE from a basic record (the `global_values.set_*`
/// calls of the steps).
pub fn set_globals_from_basic(rcd: &BasicRcd, max_is: i64) {
    gv::set_read_length(rcd.rlth);
    gv::set_insert_size(max_is as f64);
    gv::set_average_cov(rcd.cov);
    gv::set_average_cov_is_int(matches!(rcd.cov_py, PyNum::Int(_)));
}

/// Which entry point is running: `x_TEA_main.py`, or `x_TEA_main_germline_cloud.py` (which
/// differs in `--gntp_feature`, `--gVCF`, `--gene`, has no `--igv` mode and takes a sample id).
#[derive(Clone, Debug)]
pub enum MainVariant {
    Local,
    GermlineCloud { sample_id: String },
}

/// Command line of `xtea_main_germline_cloud` (port of `x_TEA_main_germline_cloud.py`): the
/// options of `x_TEA_main.py` plus `--hard` and `--id`.
#[derive(Parser, Debug, Clone, Default)]
#[command(
    name = "xtea_main_germline_cloud",
    about = "xTea germline calling for the cloud workflow (Rust re-implementation of x_TEA_main_germline_cloud.py)",
    infer_long_args = true,
    args_override_self = true,
    disable_help_flag = true
)]
pub struct CloudOpts {
    #[command(flatten)]
    pub base: Opts,
    #[arg(long = "hard", help = "This is hard-cut for fitering out coverage abnormal candidates")]
    pub hard: bool,
    #[arg(allow_hyphen_values = true, long = "id", default_value = "sampleID", value_name = "SAMPLE_ID", help = "Sample ID")]
    pub sample_id: String,
}

/// `x_TEA_main.py`
pub fn run(opts: Opts) -> Result<()> {
    run_variant(opts, &MainVariant::Local)
}

/// `x_TEA_main_germline_cloud.py`
pub fn run_cloud(opts: CloudOpts) -> Result<()> {
    run_variant(opts.base, &MainVariant::GermlineCloud { sample_id: opts.sample_id })
}

pub fn run_variant(opts: Opts, variant: &MainVariant) -> Result<()> {
    let cloud = matches!(variant, MainVariant::GermlineCloud { .. });
    crate::rng::set_seed(opts.seed);
    if opts.mit {
        gv::turn_on_mit();
    }
    if opts.dna {
        gv::turn_off_rna_mediated();
    }
    if opts.cbs {
        gv::turn_on_check_by_sample();
    }
    if opts.sva {
        gv::turn_on_sva();
    }
    let b_automatic = !opts.user_specific;
    let b_tumor = opts.tumor;
    let f_purity = opts.purity;
    let b_resume = opts.resume;
    let n_jobs = opts.cores.max(1) as usize;

    if opts.preprocess {
        let sf_ref = req(&opts.reference, "-r")?;
        let sf_annotation = req(&opts.annotation, "-a")?;
        let sf_out_fa = req(&opts.output, "-o")?;
        let flank_lth = opts.extend;
        let b_with_flank = opts.withflank;
        let mut x_annotation = XAnnotation::new(sf_annotation);
        gv::set_load_rmsk_left_extnd(opts.rmsk_extnd);
        if opts.bed {
            x_annotation.collect_seqs_of_te_from_ref_bed_fmt(sf_ref, sf_out_fa, flank_lth)?;
        } else {
            let b_with_chr = x_annotation.is_ref_chrm_with_chr(sf_ref)?;
            x_annotation.set_with_chr(b_with_chr);
            x_annotation.load_rmsk_annotation()?;
            x_annotation.collect_seqs_of_te_from_ref(sf_ref, flank_lth, b_with_flank, sf_out_fa)?;
            x_annotation.bwa_index_te_seqs(sf_out_fa);
        }
    } else if opts.flank.as_deref().map_or(false, |s| !s.is_empty()) {
        let sf_ref = req(&opts.reference, "-r")?;
        let sf_annotation = req(&opts.annotation, "-a")?;
        let sf_out_fa = req(&opts.output, "-o")?;
        let flank_lth = opts.extend;
        let mut x_annotation = XAnnotation::new(sf_annotation);
        let b_with_chr = x_annotation.is_ref_chrm_with_chr(sf_ref)?;
        x_annotation.set_with_chr(b_with_chr);
        if opts.bed {
            println!("load from bed");
            x_annotation.load_annotation_no_extnd_from_bed()?;
        } else {
            println!("load from rmsk");
            x_annotation.load_rmsk_annotation_no_extnd()?;
        }
        x_annotation.collect_flank_regions_of_te_from_ref(sf_ref, flank_lth, sf_out_fa)?;
        x_annotation.bwa_index_te_seqs(sf_out_fa);
    } else if opts.clip {
        crate::steps::run_clip(&opts, b_automatic, b_tumor, f_purity, b_resume, n_jobs)?;
    } else if opts.discordant {
        crate::steps::run_discordant(&opts, b_automatic, b_tumor, f_purity, b_resume, n_jobs)?;
    } else if opts.filter_csn {
        crate::steps::run_filter_csn(&opts, b_automatic, b_tumor, f_purity, b_resume, n_jobs)?;
    } else if opts.transduction {
        crate::steps::run_transduction(&opts, b_automatic, b_tumor, f_purity, b_resume, n_jobs)?;
    } else if opts.sibling {
        crate::steps::run_sibling(&opts, b_automatic, b_tumor, f_purity, n_jobs)?;
    } else if opts.case_control {
        crate::steps::run_case_control(&opts, b_automatic, b_tumor, n_jobs)?;
    } else if opts.mosaic {
        crate::steps::run_mosaic_clip(&opts, b_automatic, n_jobs)?;
    } else if opts.barcode {
        crate::steps::run_x10_unsupported("-B/--barcode")?;
    } else if opts.post_f {
        crate::steps::run_post_filter(&opts, b_tumor, n_jobs)?;
    } else if opts.gntp_feature {
        let sf_bam_list = req(&opts.bam, "-b")?;
        let sf_ref = req(&opts.r#ref, "--ref")?;
        let sf_output = req(&opts.output, "-o")?;
        let extnd = if cloud {
            // the cloud version sizes the region from the estimated insert size
            let (_rcd, basic_rcd) = automatic_gnrt_parameters(sf_bam_list, sf_ref, &opts.wfolder, n_jobs, false, b_tumor, f_purity)?;
            let rlth = basic_rcd.rlth;
            let mean_is = basic_rcd.mean_is;
            gv::set_insert_size(crate::py::fmt::py_int(mean_is) as f64);
            let std_var = basic_rcd.std_var;
            crate::py::fmt::py_int(mean_is + 3.0 * std_var) + crate::py::fmt::py_int(rlth)
        } else {
            450
        };
        let x_gntper = crate::genotype_feature::XGenotyper::new(sf_ref, &opts.wfolder, n_jobs);
        x_gntper.call_genotype(sf_bam_list, &opts.input, extnd, sf_output)?;
    } else if opts.gntp_classify {
        crate::steps::run_gntp_classify(&opts)?;
    } else if opts.g_vcf {
        match variant {
            MainVariant::GermlineCloud { sample_id } => crate::steps::run_gvcf_with_sample_id(&opts, Some(sample_id))?,
            MainVariant::Local => crate::steps::run_gvcf(&opts)?,
        }
    } else if opts.joint {
        crate::steps::run_joint(&opts, n_jobs)?;
    } else if opts.igv && !cloud {
        crate::steps::run_igv(&opts, n_jobs)?;
    } else if opts.collect {
        crate::steps::run_x10_unsupported("-E/--collect")?;
    } else if opts.mutation {
        crate::steps::run_x10_unsupported("-I/--mutation")?;
    } else if opts.gene {
        crate::steps::run_gene_variant(&opts, !cloud)?;
    } else if opts.collect_illumina {
        // nothing to do (as in the original)
    } else if opts.assembly {
        crate::steps::run_x10_unsupported("-A/--assembly")?;
    } else if opts.map {
        crate::steps::run_x10_unsupported("-M/--map")?;
    } else if opts.flk_map {
        crate::steps::run_x10_unsupported("--flk_map")?;
    } else if opts.filter_asm {
        crate::steps::run_x10_unsupported("-F/--filter_asm")?;
    } else if opts.collect_clip {
        let sf_bam_list = &opts.input;
        let sf_ref = opts.r#ref.clone().unwrap_or_default();
        let sf_annotation = opts.annotation.clone().unwrap_or_default();
        let tem_locator = crate::tei_locator::TeMultiLocator::new(sf_bam_list, &opts.wfolder, n_jobs, &sf_ref);
        let wfolder_pub_clip = opts.cwfolder.clone().unwrap_or_default();
        tem_locator.collect_all_clipped_from_multiple_alignmts(&sf_annotation, opts.single, &opts.wfolder, &wfolder_pub_clip)?;
    }
    Ok(())
}
