//! `xtea_germline_cloud`: port of `gnrt_pipeline_germline_cloud_v38.py`, the single-sample
//! germline wrapper used by the cloud workflow (`xtea_germline_hg38.cwl`).  It indexes the BAM/CRAM
//! when no index exists (`pysam.index`, here htslib's index builder: `.bai` for BAM, `.crai` for
//! CRAM), unpacks the repeat library and the reference tarballs (`-D` defaults to True and cannot be
//! switched off), writes the configurations and pipeline scripts of the family, runs them with
//! `sh`, and finally copies the results into `<id>_results/` and packs them into
//! `<id>_<type>.tar.gz`.
//!
//! The generated commands call `${XTEA_MAIN}`, set to `xtea_main_germline_cloud` (the port of
//! `x_TEA_main_germline_cloud.py`) next to this executable or on `PATH`.
//!
//! Python behaviour kept: the repeat types are named `LINE1`, `ALU`, `SVA`, `HERV`,
//! `Mitochondria`, `Msta`, `Pseudogene` (folders and `<name>.config`) and all of them use the
//! standard (RNA) commands; HERV's configuration has REF before ANNOTATION1; `-n`, `-f` and `-y` have
//! no defaults (`-n` missing prints `-n None` in the commands, `-f`/`-y` missing raise); the gene
//! file is always `<work folder>/gencode.annotation.gff3` (`-g` is ignored); `print(l_sh)` prints
//! the Python list repr; every `<id>_<type>.tar.gz` archives the whole results folder.

use super::{
    annotation_steps, cp_file, decompress, exists, finish, gnrt_lib_config, gnrt_parameters, gnrt_pipelines, is_file, load_par_config,
    need_flag, need_str, parse_opts, popen, post_filter_steps_hg38, py_add_slash, py_int_or_none, py_list_repr, py_print, py_str_or_none,
    rna_steps, run_cmd, select_rep_types, type_error, write_py_file, xtea_main_path, CallArgs, RepType, PUB_CLIP, SR_ASM_INPUT_DEFAULT,
    XTEA_MAIN_GERMLINE_CLOUD_EXE,
};
use crate::py::optparse::parse_int;
use anyhow::{anyhow, Result};
use clap::Parser;

// these should be consistent with the definition within x_rep_type.py
pub const REP_TYPE_L1: &str = "LINE1";
pub const REP_TYPE_ALU: &str = "ALU";
pub const REP_TYPE_SVA: &str = "SVA";
pub const REP_TYPE_HERV: &str = "HERV";
pub const REP_TYPE_MIT: &str = "Mitochondria";
pub const REP_TYPE_MSTA: &str = "Msta";
pub const REP_TYPE_PSEUDOGENE: &str = "Pseudogene";

/// Repeat types and library files of `gnrt_lib_config` (hg38), in `-y` bit order.
pub static CLOUD_REP_TYPES: [RepType; 7] = [
    RepType {
        name: REP_TYPE_L1,
        bit: 1,
        anno: "LINE/hg38/hg38_L1_larger_500_with_all_L1HS.out",
        anno1: "LINE/hg38/hg38_L1.fa.out",
        copy_with_flank: "LINE/hg38/hg38_L1HS_copies_larger_5K_with_flank.fa",
        flank: Some("LINE/hg38/hg38_FL_L1_flanks_3k.fa"),
        cns: "consensus/LINE1.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_ALU,
        bit: 2,
        anno: "Alu/hg38/hg38_Alu.out",
        anno1: "Alu/hg38/hg38_Alu.out",
        copy_with_flank: "Alu/hg38/hg38_AluJabc_copies_with_flank.fa",
        flank: None,
        cns: "consensus/ALU.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_SVA,
        bit: 4,
        anno: "SVA/hg38/hg38_SVA.out",
        anno1: "SVA/hg38/hg38_SVA.out",
        copy_with_flank: "SVA/hg38/hg38_SVA_copies_with_flank.fa",
        flank: Some("SVA/hg38/hg38_FL_SVA_flanks_3k.fa"),
        cns: "consensus/SVA.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_HERV,
        bit: 8,
        anno: "HERV/hg38/hg38_HERV.out",
        anno1: "HERV/hg38/hg38_HERV.out",
        copy_with_flank: "HERV/hg38/hg38_HERV_copies_with_flank.fa",
        flank: None,
        cns: "consensus/HERV.fa",
        ref_before_anno1: true,
    },
    RepType {
        name: REP_TYPE_MIT,
        bit: 16,
        anno: "Mitochondrion/hg38/hg38_numtS.out",
        anno1: "Mitochondrion/hg38/hg38_numtS.out",
        copy_with_flank: "Mitochondrion/hg38/hg38_mitochondrion_copies_with_flank.fa",
        flank: None,
        cns: "consensus/mitochondrion.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_MSTA,
        bit: 32,
        anno: "MSTA/hg38/hg38_MSTA.out",
        anno1: "MSTA/hg38/hg38_MSTA.out",
        copy_with_flank: "MSTA/hg38/hg38_MSTA_copies_with_flank.fa",
        flank: None,
        cns: "consensus/MSTA.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_PSEUDOGENE,
        bit: 64,
        anno: "Pseudogene/hg38/gencode_v28_GRCh38_exon_annotation.out",
        anno1: "Pseudogene/hg38/gencode_v28_GRCh38_exon_annotation.out",
        copy_with_flank: "Pseudogene/hg38/gencode_v28_GRCh38_transcript_masked.fa",
        flank: None,
        cns: "consensus/gencode_v28_GRCh38_transcript_masked.fa",
        ref_before_anno1: false,
    },
];

/// Result files copied by `cp_compress_results` (followed by `<id>_<type>.vcf`).
const RESULT_FILES: [&str; 14] = [
    "candidate_disc_filtered_cns.txt",
    "candidate_list_from_clip.txt",
    "candidate_list_from_clip.txt_tmp",
    "candidate_list_from_disc.txt",
    "candidate_list_from_disc.txt.clip_sites_raw_disc.txt",
    "candidate_disc_filtered_cns.txt.before_filtering",
    "candidate_disc_filtered_cns.txt.gntp.features",
    "candidate_disc_filtered_cns.txt.high_confident",
    "candidate_disc_filtered_cns2.txt",
    "candidate_disc_filtered_cns2.txt.high_confident",
    "candidate_disc_filtered_cns_post_filtering.txt",
    "candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt",
    "candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt",
    "candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt",
];

/// `gnrt_calling_command` of the cloud script (steps 1, 2, 16, 256, 512, 1024 and 8192 only).
fn gnrt_calling_command(s_family_id: &str, a: &CallArgs, b_sva: bool) -> Result<String> {
    let s = rna_steps(a, b_sva, SR_ASM_INPUT_DEFAULT);
    let (sf_post_filter, sf_post_filter_hc) = post_filter_steps_hg38(a);
    let iflag = need_flag(a.iflag)?;
    let mut s_cmd = String::new();
    if (iflag & 1) == 1 {
        s_cmd += &s.sclip_step;
    }
    if (iflag & 2) == 2 {
        s_cmd += &s.sdisc_step;
    }
    if (iflag & 16) == 16 {
        s_cmd += &s.s_filter;
    }
    if (iflag & 256) == 256 {
        s_cmd += &s.sf_trsdct;
    }
    if (iflag & 512) == 512 {
        s_cmd += &sf_post_filter;
        s_cmd += &sf_post_filter_hc;
    }
    if (iflag & 1024) == 1024 && (iflag & 2048) != 2048 {
        s_cmd += &annotation_vcf(s_family_id, &a.ncores, a.i_rep_type);
    }
    if (iflag & 8192) == 8192 {
        s_cmd += &s.sf_sr_asm;
    }
    Ok(s_cmd)
}

/// `gnrt_calling_command_annotation_vcf` of the cloud script (gVCF with `--id`, no BamSnap).
fn annotation_vcf(s_family_id: &str, ncores: &str, i_rep_type: i64) -> String {
    let an = annotation_steps(ncores, i_rep_type);
    let sf_cvt_gvcf = format!(
        concat!(
            "${{XTEA_MAIN}} --gVCF -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt\" ",
            " -o ${{PREFIX}} -b ${{BAM_LIST}} --ref ${{REF}} --rtype {0} --id {1}\n"
        ),
        i_rep_type, s_family_id
    );
    [an.sf_gene, an.sf_gntp_classify, sf_cvt_gvcf].concat()
}

/// `pysam.index(sf_bam)` (`samtools index`): BAI for BAM, CRAI for CRAM, next to the file.
fn index_bam_cram(sf_bam: &str) -> Result<()> {
    rust_htslib::bam::index::build(sf_bam, None, rust_htslib::bam::index::Type::Bai, 0)
        .map_err(|e| anyhow!("SamtoolsError: samtools index: failed to create index for \"{}\": {}", sf_bam, e))
}

/// Values of the `__main__` scope used by `gnrt_running_shell`.
struct Ctx<'a> {
    o: &'a Opts,
    xtea_main: &'a str,
}

/// `gnrt_running_shell` of the cloud script.
#[allow(clippy::too_many_arguments)]
fn gnrt_running_shell(
    ctx: &Ctx,
    s_family_id: &str,
    sf_case_bam: &str,
    l_rep_type: &[&RepType],
    s_wfolder: &str,
    sf_folder_rep: Option<&str>,
    sf_ref: Option<&str>,
    sf_gene: &str,
    sf_black_list: &str,
    sf_folder_xtea: Option<&str>,
    ncores: Option<i64>,
    f_purity: f64,
    b_tumor: bool,
    sf_submit_sh: &str,
) -> Result<()> {
    let o = ctx.o;
    let (b_mosaic, b_user_par, b_force, b_hard_cut, b_resume) = (false, false, false, false, false);

    let s_wfolder = py_add_slash(Some(s_wfolder))?;
    if !exists(&s_wfolder) {
        py_print(&format!("Generate working folder: {}", s_wfolder));
        popen(&format!("mkdir {}", s_wfolder))?;
    }
    let sf_folder = format!("{}{}", s_wfolder, s_family_id); // first create folder
    if !exists(&sf_folder) {
        py_print(&format!("Generate family folder: {}", sf_folder));
        popen(&format!("mkdir {}", sf_folder))?;
        let sf_pub_clip = format!("{}/{}/", sf_folder, PUB_CLIP);
        popen(&format!("mkdir {}", sf_pub_clip))?;
    }

    // the sample and bam files
    for rep in l_rep_type {
        let s_wfolder_rep = format!("{}/{}", sf_folder, rep.name);
        if !exists(&s_wfolder_rep) {
            run_cmd(&format!("mkdir {}", s_wfolder_rep))?;
        }
        let sf_rep_sample_id = format!("{}/sample_id.txt", s_wfolder_rep);
        write_py_file(&sf_rep_sample_id, s_family_id)?;
        py_print(&format!("Sample id file: {}", sf_rep_sample_id));
        if is_file(sf_case_bam) {
            let sf_rep_bam = format!("{}/bam_list1.txt", s_wfolder_rep);
            write_py_file(&sf_rep_bam, &format!("{}\t{}\n", s_family_id, sf_case_bam))?;
        }
    }

    let mut l_sh: Vec<String> = Vec::new();
    let sf_sample_folder = format!("{}{}/", s_wfolder, s_family_id);
    let sf_pub_clip = format!("{}{}/", sf_sample_folder, PUB_CLIP);
    gnrt_lib_config(l_rep_type, sf_folder_rep, sf_ref, Some(sf_gene), sf_black_list, sf_folder_xtea, &sf_sample_folder, ctx.xtea_main)?;
    for rep in l_rep_type {
        let sf_config = format!("{}{}.config", sf_sample_folder, rep.name);
        let s_wfolder_rep = format!("{}{}", sf_sample_folder, rep.name);
        if !exists(&s_wfolder_rep) {
            run_cmd(&format!("mkdir {}", s_wfolder_rep))?;
        }
        let sf_rep_sample_id = format!("{}/sample_id.txt", s_wfolder_rep);
        let mut sf_rep_bam = format!("{}/bam_list1.txt", s_wfolder_rep);
        if !is_file(&sf_rep_bam) {
            sf_rep_bam = "null".to_string();
        }
        let sf_rep_x10_bam = "null";
        let s_head = "";
        let l_libs = load_par_config(&sf_config)?;
        let s_libs = gnrt_parameters(&l_libs);
        let a = CallArgs {
            iclip_c: o.nclip,
            iclip_rp: o.cliprep,
            idisc_c: o.ndisc,
            iflt_clip: o.nfilterclip,
            iflt_disc: o.nfilterdisc,
            ncores: py_int_or_none(ncores),
            iflk_len: o.flklen,
            min_tei_len: o.teilen,
            iflag: o.flag,
            b_mosaic,
            b_user_par,
            b_force,
            b_tumor,
            b_hard_cut,
            b_resume,
            f_purity,
            i_rep_type: rep.bit,
            s_cfolder: sf_pub_clip.clone(),
        };
        let mut s_calling_cmd = gnrt_calling_command(s_family_id, &a, false)?;
        if rep.name == REP_TYPE_SVA {
            s_calling_cmd = gnrt_calling_command(s_family_id, &a, true)?;
        }
        let l_tmp_sh =
            gnrt_pipelines(s_head, &s_libs, &s_calling_cmd, &sf_rep_sample_id, &sf_rep_bam, sf_rep_x10_bam, &sf_sample_folder, rep.name)?;
        l_sh.extend(l_tmp_sh);
    }
    py_print("List of command:\n");
    py_print(&py_list_repr(&l_sh));
    let mut s_submit = String::from("#!/bin/bash\n\n");
    for s_sh in &l_sh {
        s_submit += &format!("chmod +x {}\n", s_sh);
        s_submit += &format!("sh {}\n", s_sh);
    }
    write_py_file(sf_submit_sh, &s_submit)
}

/// `cp_compress_results` of the cloud script.
fn cp_compress_results(s_wfolder: &str, l_rep_type: &[&RepType], sample_id: &str) -> Result<()> {
    // create a "results" folder
    let sf_rslts = format!("{}{}_results/", s_wfolder, sample_id);
    if !exists(&sf_rslts) {
        run_cmd(&format!("mkdir {}", sf_rslts))?;
    }
    for rep in l_rep_type {
        let sf_rslts_rep_folder = format!("{}{}/", sf_rslts, rep.name);
        if !exists(&sf_rslts_rep_folder) {
            run_cmd(&format!("mkdir {}", sf_rslts_rep_folder))?;
        }
        let sf_samp_folder = format!("{}{}/", sf_rslts_rep_folder, sample_id);
        if !exists(&sf_samp_folder) {
            run_cmd(&format!("mkdir {}", sf_samp_folder))?;
        }
        let sf_source_folder = format!("{}{}/{}/", s_wfolder, sample_id, rep.name);
        for name in RESULT_FILES {
            cp_file(&format!("{}{}", sf_source_folder, name), &sf_samp_folder)?;
        }
        cp_file(&format!("{}{}_{}.vcf", sf_source_folder, sample_id, rep.name), &sf_samp_folder)?;
        // compress the results folder to one file
        let sf_compressed = format!("{}{}_{}.tar.gz", s_wfolder, sample_id, rep.name);
        run_cmd(&format!("tar -cvzf {} -C {} .", sf_compressed, sf_rslts))?;
    }
    Ok(())
}

/// The `__main__` block of the cloud script.
pub fn run(o: &Opts) -> Result<()> {
    let xtea_main = xtea_main_path(XTEA_MAIN_GERMLINE_CLOUD_EXE);
    let s_family_id = o.id.as_deref();
    let sf_case_bam = o.bam.as_deref(); // case bam file (not file list, but the exact bam/cram file)

    let mut s_wfolder = o.wfolder.clone().ok_or_else(|| anyhow!("AttributeError: 'NoneType' object has no attribute 'split'"))?;
    let wf_fields: Vec<&str> = s_wfolder.split('/').collect();
    if !wf_fields.is_empty() && wf_fields[0] == "." {
        let cwd = std::env::current_dir().map_err(|e| anyhow!("FileNotFoundError: {}", e))?;
        s_wfolder = format!("{}/{}", cwd.to_string_lossy(), wf_fields[1..].join("/"));
    }
    let s_wfolder = py_add_slash(Some(&s_wfolder))?;

    // index the bam/cram
    let sf_case_bam = sf_case_bam.ok_or_else(|| type_error("unsupported operand type(s) for +: 'NoneType' and 'str'"))?;
    if !is_file(&format!("{}.bai", sf_case_bam)) && !is_file(&format!("{}.crai", sf_case_bam)) {
        index_bam_cram(sf_case_bam)?;
    }
    py_print("Finished indexing bam files");

    let ncores = o.cores;
    let mut sf_folder_rep = o.lib.clone(); // the lib folder path
    let mut sf_ref = o.r#ref.clone(); // reference genome
    if o.decompress {
        decompress(py_str_or_none(o.lib.as_deref()), &s_wfolder)?;
        decompress(py_str_or_none(o.r#ref.as_deref()), &s_wfolder)?;
        sf_folder_rep = Some(format!("{}rep_lib_annotation/", s_wfolder));
        sf_ref = Some(format!("{}genome.fa", s_wfolder));
    }
    let sf_gene = format!("{}gencode.annotation.gff3", s_wfolder);

    let i_rep_type = need_flag(o.rep_type)?;
    let l_rep_type = select_rep_types(&CLOUD_REP_TYPES, i_rep_type);
    let s_family_id = need_str(s_family_id)?;
    let sf_sbatch_sh = format!("{}{}/run_xTea_germline_ins.sh", s_wfolder, s_family_id);
    // first, run the jobs separately for all the case and control samples
    py_print("Start command generating:");

    let f_purity = 1.0;
    let b_tumor = false;
    let ctx = Ctx { o, xtea_main: &xtea_main };
    gnrt_running_shell(
        &ctx,
        s_family_id,
        sf_case_bam,
        &l_rep_type,
        &s_wfolder,
        sf_folder_rep.as_deref(),
        sf_ref.as_deref(),
        &sf_gene,
        &o.blacklist,
        o.xtea.as_deref(),
        ncores,
        f_purity,
        b_tumor,
        &sf_sbatch_sh,
    )?;

    let cmd = format!("sh {}", sf_sbatch_sh);
    py_print("Start run the generated command:");
    run_cmd(&cmd)?;

    py_print("Finished running command, begin to compress files");
    cp_compress_results(&s_wfolder, &l_rep_type, s_family_id)?;
    py_print("Finished results compression");
    Ok(())
}

/// Options of `gnrt_pipeline_germline_cloud_v38.py` (no `-V`; `-x` is the xTea folder here).
#[derive(Parser, Debug, Clone)]
#[command(
    name = "xtea_germline_cloud",
    about = "Germline xTea wrapper for the cloud workflow on hg38 (Rust port of gnrt_pipeline_germline_cloud_v38.py)",
    infer_long_args = true,
    args_override_self = true,
    disable_help_flag = true
)]
pub struct Opts {
    #[arg(short = 'h', long = "help", help = "show this help message and exit")]
    pub help: bool,
    // optparse default is True, so the flag cannot be turned off
    #[arg(short = 'D', long = "decompress", default_value_t = true, help = "Decompress the rep lib and reference file")]
    pub decompress: bool,
    #[arg(allow_hyphen_values = true, short = 'i', long = "id", value_name = "FILE", help = "sample id list file ")]
    pub id: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'a', long = "par", value_name = "FILE", help = "parameter file ")]
    pub parameters: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'l', long = "lib", value_name = "FILE", help = "TE lib config file ")]
    pub lib: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'b', long = "bam", value_name = "FILE", help = "Input bam file")]
    pub bam: Option<String>,
    #[arg(allow_hyphen_values = true, long = "pa", default_value = "null", value_name = "FILE", help = "Input paternal bam file")]
    pub paternal: String,
    #[arg(allow_hyphen_values = true, long = "ma", default_value = "null", value_name = "FILE", help = "Input maternal bam file")]
    pub maternal: String,

    #[arg(allow_hyphen_values = true, short = 'p', long = "path", help = "Working folder")]
    pub wfolder: Option<String>,
    #[arg(allow_hyphen_values = true, value_parser = parse_int, short = 'n', long = "cores", help = "number of cores")]
    pub cores: Option<i64>,
    #[arg(allow_hyphen_values = true, short = 'r', long = "ref", help = "reference genome")]
    pub r#ref: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'g', long = "gene", help = "Gene annotation file")]
    pub gene: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'x', long = "xtea", help = "xTEA folder")]
    pub xtea: Option<String>,

    #[arg(
        allow_hyphen_values = true,
        value_parser = parse_int,
        short = 'f',
        long = "flag",
        help = "Flag indicates which step to run (1-clip, 2-disc, 4-barcode, 8-xfilter, 16-filter, 32-asm)"
    )]
    pub flag: Option<i64>,
    #[arg(
        allow_hyphen_values = true,
        value_parser = parse_int,
        short = 'y',
        long = "reptype",
        help = "Type of repeats working on: 1-L1, 2-Alu, 4-SVA, 8-HERV, 16-Mitochondrial"
    )]
    pub rep_type: Option<i64>,

    #[arg(allow_hyphen_values = true, value_parser = parse_int, long = "flklen", default_value_t = 3000, help = "flank region file")]
    pub flklen: i64,
    #[arg(allow_hyphen_values = true, value_parser = parse_int, long = "nclip", default_value_t = 3, help = "cutoff of minimum # of clipped reads")]
    pub nclip: i64,
    #[arg(
        allow_hyphen_values = true,
        value_parser = parse_int,
        long = "cr",
        default_value_t = 1,
        help = "cutoff of minimum # of clipped reads whose mates map in repetitive regions"
    )]
    pub cliprep: i64,
    #[arg(allow_hyphen_values = true, value_parser = parse_int, long = "nd", default_value_t = 5, help = "cutoff of minimum # of discordant pair")]
    pub ndisc: i64,
    #[arg(
        allow_hyphen_values = true,
        value_parser = parse_int,
        long = "nfclip",
        default_value_t = 3,
        help = "cutoff of minimum # of clipped reads in filtering step"
    )]
    pub nfilterclip: i64,
    #[arg(
        allow_hyphen_values = true,
        value_parser = parse_int,
        long = "nfdisc",
        default_value_t = 5,
        help = "cutoff of minimum # of discordant pair of each sample in filtering step"
    )]
    pub nfilterdisc: i64,
    #[arg(
        allow_hyphen_values = true,
        value_parser = parse_int,
        long = "teilen",
        default_value_t = 50,
        help = "minimum length of the insertion for future analysis"
    )]
    pub teilen: i64,

    #[arg(
        allow_hyphen_values = true,
        long = "blacklist",
        default_value = "null",
        value_name = "FILE",
        help = "Reference panel database for filtering, or a blacklist region"
    )]
    pub blacklist: String,
    #[arg(allow_hyphen_values = true, short = 'o', long = "output", value_name = "FILE", help = "The output file")]
    pub output: Option<String>,
    /// Positional arguments are accepted and ignored (optparse collects them separately).
    #[arg(hide = true)]
    pub extra_args: Vec<String>,
}

/// Entry point of the `xtea_germline_cloud` binary; returns the exit status.
pub fn main() -> i32 {
    let opts: Opts = parse_opts("xtea_germline_cloud");
    finish(run(&opts))
}
