//! `xtea_hg19`: port of `gnrt_pipeline_local.py` (identical to `bin/xtea_hg19`), the generator
//! of the calling pipelines for Illumina short reads on hg19.
//!
//! Differences to the hg38 script ([`super::local_v38`]) that are reproduced here: no `--hard`
//! option, no Pseudogene type (`-y 64` selects nothing), other option defaults (`-n 8`, `-m 25`,
//! `-q medium`, `-t 0-15:00`, `-f 19`, no default for `-g`), no BamSnap steps in the annotation
//! stage, `--postFmosaic` steps with `--blacklist`, a gene-annotation step (flag 512) instead of
//! post-filter + annotation for HERV, MSTA and Mitochondrion, `--ref "hg19"` in the igv steps,
//! no 32768 flag in de novo mode, the HERV configuration written in the normal order, the
//! `len(fields)<2` check of `split_bam_list`, and the `.case_list` file as BAM list in de novo
//! mode.

use super::{
    annotation_steps, denovo_from_trio, dna_steps, finish, need_flag, parse_opts, rna_steps, run_local, somatic_case_control, CallArgs,
    LocalOptions, LocalVariant, RepType, SR_ASM_INPUT_DEFAULT,
};
use crate::py::optparse::{parse_float, parse_int};
use anyhow::Result;
use clap::Parser;

pub const REP_TYPE_L1: &str = "L1";
pub const REP_TYPE_ALU: &str = "Alu";
pub const REP_TYPE_SVA: &str = "SVA";
pub const REP_TYPE_HERV: &str = "HERV";
pub const REP_TYPE_MIT: &str = "Mitochondrion";
pub const REP_TYPE_MSTA: &str = "MSTA";

/// Repeat types and library files of `gnrt_lib_config` (hg19), in `-y` bit order.
pub static HG19_REP_TYPES: [RepType; 6] = [
    RepType {
        name: REP_TYPE_L1,
        bit: 1,
        anno: "LINE/hg19/hg19_L1_larger500_with_all_L1HS.out",
        anno1: "LINE/hg19/hg19_L1.fa.out",
        copy_with_flank: "LINE/hg19/hg19_L1HS_copies_larger_5K_with_flank.fa",
        flank: Some("LINE/hg19/hg19_FL_L1_flanks_3k.fa"),
        cns: "consensus/LINE1.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_ALU,
        bit: 2,
        anno: "Alu/hg19/hg19_Alu.out",
        anno1: "Alu/hg19/hg19_Alu.out",
        copy_with_flank: "Alu/hg19/hg19_AluJabc_copies_with_flank.fa",
        flank: None,
        cns: "consensus/ALU.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_SVA,
        bit: 4,
        anno: "SVA/hg19/hg19_SVA.out",
        anno1: "SVA/hg19/hg19_SVA.out",
        copy_with_flank: "SVA/hg19/hg19_SVA_copies_with_flank.fa",
        flank: Some("SVA/hg19/hg19_FL_SVA_flanks_3k.fa"),
        cns: "consensus/SVA.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_HERV,
        bit: 8,
        anno: "HERV/hg19/hg19_HERV.out",
        anno1: "HERV/hg19/hg19_HERV.out",
        copy_with_flank: "HERV/hg19/hg19_HERV_copies_with_flank.fa",
        flank: None,
        cns: "consensus/HERV.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_MIT,
        bit: 16,
        anno: "Mitochondrion/hg19/hg19_numtS.out",
        anno1: "Mitochondrion/hg19/hg19_numtS.out",
        copy_with_flank: "Mitochondrion/hg19/hg19_mitochondrion_copies_with_flank.fa",
        flank: None,
        cns: "consensus/mitochondrion.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_MSTA,
        bit: 32,
        anno: "MSTA/hg19/hg19_MSTA.out",
        anno1: "MSTA/hg19/hg19_MSTA.out",
        copy_with_flank: "MSTA/hg19/hg19_MSTA_copies_with_flank.fa",
        flank: None,
        cns: "consensus/MSTA.fa",
        ref_before_anno1: false,
    },
];

/// The hg19 script.
pub struct Hg19;

/// `gnrt_calling_command`
fn gnrt_calling_command(a: &CallArgs, b_sva: bool) -> Result<String> {
    let s = rna_steps(a, b_sva, SR_ASM_INPUT_DEFAULT);
    let (i_rep_type, ncores, s_tumor) = (a.i_rep_type, &a.ncores, a.s_tumor());
    let mut sf_post_filter = format!(
        concat!(
            "${{XTEA_MAIN}} --postF --rtype {0} -p ${{TMP_CNS}} ",
            "-n {1} -i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\" ",
            "-a ${{ANNOTATION1}} {2} ",
            "-o ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\"\n"
        ),
        i_rep_type, ncores, s_tumor
    );
    if a.b_mosaic {
        sf_post_filter = format!(
            concat!(
                "${{XTEA_MAIN}} --postF --postFmosaic --blacklist ${{BLACK_LIST}} ",
                "--rtype {0} -p ${{TMP_CNS}} ",
                "-n {1} -i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt\" ",
                "-a ${{ANNOTATION1}} {2} ",
                "-o ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\"\n"
            ),
            i_rep_type, ncores, s_tumor
        );
    }
    let mut sf_post_filter_hc = format!(
        concat!(
            "${{XTEA_MAIN}} --postF --rtype {0} -p ${{TMP_CNS}} -n {1} ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt.high_confident\" -a ${{ANNOTATION1}} ",
            "--blacklist ${{BLACK_LIST}} {2} ",
            "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt\"\n"
        ),
        i_rep_type, ncores, s_tumor
    );
    if a.b_mosaic {
        sf_post_filter_hc = format!(
            concat!(
                "${{XTEA_MAIN}} --postF --postFmosaic --blacklist ${{BLACK_LIST}} ",
                "--rtype {0} -p ${{TMP_CNS}} -n {1} {2} ",
                "-i ${{PREFIX}}\"candidate_disc_filtered_cns2.txt.high_confident\" -a ${{ANNOTATION1}} ",
                "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt\"\n"
            ),
            i_rep_type, ncores, s_tumor
        );
    }
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
    if (iflag & 4096) == 4096 {
        s_cmd += &s.sf_sibling;
    }
    if (iflag & 512) == 512 {
        s_cmd += &sf_post_filter;
        s_cmd += &sf_post_filter_hc;
    }
    if (iflag & 1024) == 1024 && (iflag & 2048) != 2048 {
        let an = annotation_steps(ncores, i_rep_type);
        if a.b_tumor {
            // gnrt_calling_command_annotation_vcf_tumor_case
            s_cmd += &an.sf_gene;
            s_cmd += &an.sf_cvt_gvcf_tumor;
        } else {
            // gnrt_calling_command_annotation_vcf
            s_cmd += &an.sf_gene;
            s_cmd += &an.sf_gntp_classify;
            s_cmd += &an.sf_cvt_gvcf;
        }
    }
    if (iflag & 8192) == 8192 {
        s_cmd += &s.sf_sr_asm;
    }
    Ok(s_cmd)
}

/// `gnrt_calling_command_non_RNA` (`s_dna` = "--dna") and `gnrt_calling_command_MT`
/// (`s_dna` = "--dna --mit"); `sf_gene` is the gene-annotation step of flag 512.
fn gnrt_calling_command_dna(a: &CallArgs, s_dna: &str, sf_gene: &str) -> Result<String> {
    let d = dna_steps(a, s_dna);
    let iflag = need_flag(a.iflag)?;
    let mut s_cmd = String::new();
    if (iflag & 1) == 1 {
        s_cmd += &d.sclip_step;
    }
    if (iflag & 2) == 2 {
        s_cmd += &d.sdisc_step;
    }
    if (iflag & 16) == 16 {
        s_cmd += &d.s_filter;
    }
    if (iflag & 512) == 512 {
        s_cmd += sf_gene;
    }
    Ok(s_cmd)
}

fn gnrt_calling_command_non_rna(a: &CallArgs) -> Result<String> {
    let sf_gene = format!(
        concat!(
            "${{XTEA_MAIN}} --gene -a ${{GENE}} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" ",
            " -n {0} -o ${{PREFIX}}\"candidate_disc_filtered_cns_with_gene.txt\"\n"
        ),
        a.ncores
    );
    gnrt_calling_command_dna(a, "--dna", &sf_gene)
}

fn gnrt_calling_command_mt(a: &CallArgs) -> Result<String> {
    let sf_gene = format!(
        concat!(
            "${{XTEA_MAIN}} --gene --dna --mit -a ${{GENE}} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" ",
            " -n {0} -o ${{PREFIX}}\"candidate_disc_filtered_cns_with_gene.txt\"\n"
        ),
        a.ncores
    );
    gnrt_calling_command_dna(a, "--dna --mit", &sf_gene)
}

impl LocalVariant for Hg19 {
    fn genome(&self) -> &'static str {
        "hg19"
    }
    fn rep_types(&self) -> &'static [RepType] {
        &HG19_REP_TYPES
    }
    fn split_bam_list_skips_short_lines(&self) -> bool {
        true
    }
    fn denovo_uses_case_list(&self) -> bool {
        true
    }
    fn calling_command(&self, a: &CallArgs, rep: &RepType) -> Result<String> {
        let mut s_calling_cmd = gnrt_calling_command(a, false)?;
        if rep.name == REP_TYPE_SVA {
            s_calling_cmd = gnrt_calling_command(a, true)?;
        }
        if rep.name == REP_TYPE_MSTA || rep.name == REP_TYPE_HERV {
            s_calling_cmd = gnrt_calling_command_non_rna(a)?;
        } else if rep.name == REP_TYPE_MIT {
            s_calling_cmd = gnrt_calling_command_mt(a)?;
        }
        Ok(s_calling_cmd)
    }
    fn somatic_case_control(&self, a: &CallArgs, sf_ctr_list: &str, sf_case_control_bam_list: &str) -> Result<String> {
        somatic_case_control(a, sf_ctr_list, sf_case_control_bam_list, "hg19")
    }
    fn denovo_from_trio(&self, a: &CallArgs, sf_ctr_list: &str, sf_case_control_bam_list: &str) -> Result<String> {
        denovo_from_trio(a, sf_ctr_list, sf_case_control_bam_list, false)
    }
}

/// Options of `gnrt_pipeline_local.py`.
#[derive(Parser, Debug, Clone)]
#[command(
    name = "xtea_hg19",
    about = "Generate the xTea calling pipelines for Illumina short reads on hg19 (Rust port of gnrt_pipeline_local.py)",
    infer_long_args = true,
    args_override_self = true,
    disable_help_flag = true
)]
pub struct Opts {
    #[arg(short = 'h', long = "help", help = "show this help message and exit")]
    pub help: bool,
    #[arg(short = 'D', long = "decompress", help = "Decompress the rep lib and reference file")]
    pub decompress: bool,
    #[arg(short = 'M', long = "mosaic", help = "Calling mosaic events from high coverage data")]
    pub mosaic: bool,
    #[arg(short = 'C', long = "case_control", help = "Run in case control mode")]
    pub case_control: bool,
    #[arg(long = "denovo", help = "Run in de novo mode")]
    pub denovo: bool,
    #[arg(short = 'U', long = "user", help = "Use user specific parameters instead of automatically calculated ones")]
    pub user: bool,
    #[arg(long = "force", help = "Force to start from the very beginning")]
    pub force: bool,
    #[arg(long = "tumor", help = "Working on tumor samples")]
    pub tumor: bool,
    #[arg(allow_hyphen_values = true, value_parser = parse_float, long = "purity", default_value_t = 0.45, help = "Tumor purity")]
    pub purity: f64,
    #[arg(long = "lsf", help = "Indiates submit to LSF system")]
    pub lsf: bool,
    #[arg(long = "slurm", help = "Indiates submit to slurm system")]
    pub slurm: bool,
    #[arg(long = "resume", help = "Resume the running, which will skip the step if output file already exists!")]
    pub resume: bool,
    #[arg(short = 'V', long = "version", help = "Print xTea version")]
    pub version: bool,

    #[arg(allow_hyphen_values = true, short = 'i', long = "id", value_name = "FILE", help = "sample id list file ")]
    pub id: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'a', long = "par", value_name = "FILE", help = "parameter file ")]
    pub parameters: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'l', long = "lib", value_name = "FILE", help = "TE lib config file ")]
    pub lib: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'b', long = "bam", value_name = "FILE", help = "Input bam file")]
    pub bam: Option<String>,
    #[arg(
        allow_hyphen_values = true,
        short = 'x',
        long = "x10",
        value_name = "FILE",
        help = "Input 10X bam file [10X; not supported: only null is accepted]"
    )]
    pub x10: Option<String>,

    #[arg(allow_hyphen_values = true, value_parser = parse_int, short = 'n', long = "cores", default_value_t = 8, help = "number of cores")]
    pub cores: i64,
    #[arg(allow_hyphen_values = true, value_parser = parse_int, short = 'm', long = "memory", default_value_t = 25, help = "Memory limit in GB")]
    pub memory: i64,
    #[arg(allow_hyphen_values = true, short = 'q', long = "partition", default_value = "medium", help = "Which queue to run the job")]
    pub partition: String,
    #[arg(allow_hyphen_values = true, short = 't', long = "time", default_value = "0-15:00", help = "Time limit")]
    pub time: String,

    #[arg(allow_hyphen_values = true, short = 'p', long = "path", help = "Working folder")]
    pub wfolder: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'r', long = "ref", help = "reference genome")]
    pub r#ref: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'g', long = "gene", help = "Gene annotation file")]
    pub gene: Option<String>,
    #[arg(allow_hyphen_values = true, long = "xtea", help = "xTEA folder")]
    pub xtea: Option<String>,

    #[arg(
        allow_hyphen_values = true,
        value_parser = parse_int,
        short = 'f',
        long = "flag",
        default_value_t = 19,
        help = "Flag indicates which step to run (1-clip, 2-disc, 16-filter, 256-transduction, 512-post-filter, 1024-annotation, 4096-sibling; 5907 runs all; the 10X steps 4/8/32/64/128 are not supported)"
    )]
    pub flag: i64,
    #[arg(
        allow_hyphen_values = true,
        value_parser = parse_int,
        short = 'y',
        long = "reptype",
        default_value_t = 1,
        help = "Type of repeats working on: 1-L1, 2-Alu, 4-SVA, 8-HERV, 16-Mitochondrial"
    )]
    pub rep_type: i64,

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
        help = "minimum length of the insertion for future analysis [10X; unused]"
    )]
    pub teilen: i64,

    #[arg(
        allow_hyphen_values = true,
        short = 'o',
        long = "output",
        default_value = "submit_calling_jobs_for_samples.sh",
        value_name = "FILE",
        help = "The output file"
    )]
    pub output: String,
    #[arg(
        allow_hyphen_values = true,
        long = "blacklist",
        default_value = "null",
        value_name = "FILE",
        help = "Reference panel database for filtering, or a blacklist region"
    )]
    pub blacklist: String,
    /// Positional arguments are accepted and ignored (optparse collects them separately).
    #[arg(hide = true)]
    pub extra_args: Vec<String>,
}

impl From<Opts> for LocalOptions {
    fn from(o: Opts) -> Self {
        LocalOptions {
            decompress: o.decompress,
            mosaic: o.mosaic,
            case_control: o.case_control,
            denovo: o.denovo,
            user: o.user,
            force: o.force,
            hard: false,
            tumor: o.tumor,
            purity: o.purity,
            lsf: o.lsf,
            slurm: o.slurm,
            resume: o.resume,
            version: o.version,
            id: o.id,
            lib: o.lib,
            bam: o.bam,
            x10: o.x10,
            cores: o.cores,
            memory: o.memory,
            partition: Some(o.partition),
            time: o.time,
            wfolder: o.wfolder,
            r#ref: o.r#ref,
            gene: o.gene,
            xtea: o.xtea,
            flag: Some(o.flag),
            rep_type: o.rep_type,
            flklen: o.flklen,
            nclip: o.nclip,
            cliprep: o.cliprep,
            ndisc: o.ndisc,
            nfilterclip: o.nfilterclip,
            nfilterdisc: o.nfilterdisc,
            teilen: o.teilen,
            output: o.output,
            blacklist: o.blacklist,
        }
    }
}

/// Entry point of the `xtea_hg19` binary; returns the exit status.
pub fn main() -> i32 {
    let opts: Opts = parse_opts("xtea_hg19");
    finish(run_local(&Hg19, &opts.into()))
}
