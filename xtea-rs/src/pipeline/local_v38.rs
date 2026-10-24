//! `xtea`: port of `gnrt_pipeline_local_v38.py` (identical to `bin/xtea`), the generator of the
//! calling pipelines for Illumina short reads on hg38.  [`V38Family`] also serves the CHM13
//! copy of the script (see [`super::local_chm13`]).
//!
//! Differences to the hg19 script that are reproduced here: the `--hard` option (passed on to
//! the post-filter steps), the Pseudogene repeat type (`-y 64`), the BamSnap steps of the
//! annotation stage, post-filter + annotation steps (512/1024) for HERV, MSTA and
//! Mitochondrion, the `--postFmosaic` steps without `--blacklist`, the 32768 flag of the
//! de novo mode, the HERV configuration written with REF before ANNOTATION1, no
//! `len(fields)<2` check in `split_bam_list`, and the original BAM list in de novo mode.

use super::{
    annotation_steps, denovo_from_trio, dna_steps, finish, need_flag, parse_opts, post_filter_steps_hg38, rna_steps, run_local,
    somatic_case_control, CallArgs, LocalOptions, LocalVariant, RepType, SR_ASM_INPUT_DEFAULT,
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
pub const REP_TYPE_PSEUDOGENE: &str = "Pseudogene";

/// Repeat types and library files of `gnrt_lib_config` (hg38), in `-y` bit order.
pub static HG38_REP_TYPES: [RepType; 7] = [
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

/// The hg38 script and its CHM13 copy, which differs only in the library files, the reference
/// name of the igv steps, the input of the `--sr_asm` step and the `-V` text.
pub struct V38Family {
    pub genome: &'static str,
    pub igv_ref: &'static str,
    pub sr_asm_input: &'static str,
    pub rep_types: &'static [RepType],
}

pub static HG38: V38Family = V38Family { genome: "hg38", igv_ref: "hg38", sr_asm_input: SR_ASM_INPUT_DEFAULT, rep_types: &HG38_REP_TYPES };

impl V38Family {
    /// `gnrt_calling_command`
    fn gnrt_calling_command(&self, a: &CallArgs, b_sva: bool) -> Result<String> {
        let s = rna_steps(a, b_sva, self.sr_asm_input);
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
        if (iflag & 4096) == 4096 {
            s_cmd += &s.sf_sibling;
        }
        if (iflag & 512) == 512 {
            s_cmd += &sf_post_filter;
            s_cmd += &sf_post_filter_hc;
        }
        if (iflag & 1024) == 1024 && (iflag & 2048) != 2048 {
            if a.b_tumor {
                s_cmd += &annotation_vcf_tumor_case(&a.ncores, a.i_rep_type);
            } else {
                s_cmd += &annotation_vcf(&a.ncores, a.i_rep_type);
            }
        }
        if (iflag & 8192) == 8192 {
            s_cmd += &s.sf_sr_asm;
        }
        Ok(s_cmd)
    }
}

/// `gnrt_calling_command_non_RNA` (`s_dna` = "--dna") and `gnrt_calling_command_MT`
/// (`s_dna` = "--dna --mit"), which differ only by that option.
fn gnrt_calling_command_dna(a: &CallArgs, s_dna: &str) -> Result<String> {
    let d = dna_steps(a, s_dna);
    let (i_rep_type, ncores, s_hard_cut) = (a.i_rep_type, &a.ncores, a.s_hard_cut());
    let mut sf_post_filter = format!(
        concat!(
            "${{XTEA_MAIN}} --postF --rtype {0} -p ${{TMP_CNS}} ",
            "-n {1} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" ",
            "-a ${{ANNOTATION1}} {2} ",
            "-o ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\"\n"
        ),
        i_rep_type, ncores, s_hard_cut
    );
    if a.b_mosaic {
        sf_post_filter = format!(
            concat!(
                "${{XTEA_MAIN}} --postF --postFmosaic --rtype {0} -p ${{TMP_CNS}} ",
                "-n {1} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" ",
                "-a ${{ANNOTATION1}} ",
                "-o ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\"\n"
            ),
            i_rep_type, ncores
        );
    }
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
        s_cmd += &sf_post_filter;
    }
    if (iflag & 1024) == 1024 && (iflag & 2048) != 2048 {
        s_cmd += &annotation_vcf_non_rt(ncores, i_rep_type);
    }
    Ok(s_cmd)
}

/// `gnrt_calling_command_annotation_vcf`
fn annotation_vcf(ncores: &str, i_rep_type: i64) -> String {
    let s = annotation_steps(ncores, i_rep_type);
    let (s_bamsnap, s_bamsnap_hc) = bamsnap_steps(ncores);
    [s.sf_gene, s.sf_gntp_classify, s.sf_cvt_gvcf, s_bamsnap, s_bamsnap_hc].concat()
}

/// `gnrt_calling_command_annotation_vcf_tumor_case`
fn annotation_vcf_tumor_case(ncores: &str, i_rep_type: i64) -> String {
    let s = annotation_steps(ncores, i_rep_type);
    let (s_bamsnap, s_bamsnap_hc) = bamsnap_steps(ncores);
    [s.sf_gene, s.sf_cvt_gvcf_tumor, s_bamsnap, s_bamsnap_hc].concat()
}

/// The two BamSnap steps of the annotation stage.
fn bamsnap_steps(ncores: &str) -> (String, String) {
    let s_bamsnap = format!(
        concat!(
            "${{XTEA_MAIN}} --igv --bamsnap --single_sample -p ${{PREFIX}}\"tmp/igv\" -b ${{PREFIX}}\"bam_list1.txt\" ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" ",
            " --ref ${{REF}} -e {0} -n {1} ",
            "-o ${{PREFIX}}\"tmp/igv/bamsnap_screenshot.txt\"\n"
        ),
        1000, ncores
    );
    let s_bamsnap_hc = format!(
        concat!(
            "${{XTEA_MAIN}} --igv --bamsnap --single_sample -p ${{PREFIX}}\"tmp/igv\" -b ${{PREFIX}}\"bam_list1.txt\" ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt\" ",
            " --ref ${{REF}} -e {0} -n {1} ",
            "-o ${{PREFIX}}\"tmp/igv/bamsnap_screenshot_hc.txt\"\n"
        ),
        1000, ncores
    );
    (s_bamsnap, s_bamsnap_hc)
}

/// `gnrt_calling_command_annotation_vcf_non_RT`
fn annotation_vcf_non_rt(ncores: &str, i_rep_type: i64) -> String {
    let sf_gene = format!(
        concat!(
            "${{XTEA_MAIN}} --gene -a ${{GENE}} -i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\" ",
            " -n {0} -o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.post_filtering_with_gene.txt\"\n"
        ),
        ncores
    );
    let sf_gntp_classify = format!(
        concat!(
            "${{XTEA_MAIN}} --gntp_classify -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.post_filtering_with_gene.txt\" ",
            " -n {0} --model ${{XTEA_PATH}}\"genotyping/DF21_model_1_2\" ",
            " -o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.post_filtering_with_gene_gntp.txt\"\n"
        ),
        1
    );
    let sf_cvt_gvcf = format!(
        concat!(
            "${{XTEA_MAIN}} --gVCF -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.post_filtering_with_gene_gntp.txt\" ",
            " -o ${{PREFIX}} -b ${{BAM_LIST}} --ref ${{REF}} --rtype {0}\n"
        ),
        i_rep_type
    );
    let s_bamsnap = format!(
        concat!(
            "${{XTEA_MAIN}} --igv --bamsnap --single_sample -p ${{PREFIX}}\"tmp/igv\" -b ${{PREFIX}}\"bam_list1.txt\" ",
            "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" ",
            " --ref ${{REF}} -e {0} -n {1} ",
            "-o ${{PREFIX}}\"tmp/igv/bamsnap_screenshot.txt\"\n"
        ),
        1000, ncores
    );
    [sf_gene, sf_gntp_classify, sf_cvt_gvcf, s_bamsnap].concat()
}

impl LocalVariant for V38Family {
    fn genome(&self) -> &'static str {
        self.genome
    }
    fn rep_types(&self) -> &'static [RepType] {
        self.rep_types
    }
    fn split_bam_list_skips_short_lines(&self) -> bool {
        false
    }
    fn denovo_uses_case_list(&self) -> bool {
        false
    }
    /// The calling command of one repeat type as chosen in `gnrt_running_shell` (the default
    /// command is always built first, as in the original).
    fn calling_command(&self, a: &CallArgs, rep: &RepType) -> Result<String> {
        let mut s_calling_cmd = self.gnrt_calling_command(a, false)?;
        if rep.name == REP_TYPE_SVA {
            s_calling_cmd = self.gnrt_calling_command(a, true)?;
        }
        if rep.name == REP_TYPE_HERV || rep.name == REP_TYPE_MSTA {
            s_calling_cmd = gnrt_calling_command_dna(a, "--dna")?;
        }
        if rep.name == REP_TYPE_MIT {
            s_calling_cmd = gnrt_calling_command_dna(a, "--dna --mit")?;
        }
        Ok(s_calling_cmd)
    }
    fn somatic_case_control(&self, a: &CallArgs, sf_ctr_list: &str, sf_case_control_bam_list: &str) -> Result<String> {
        somatic_case_control(a, sf_ctr_list, sf_case_control_bam_list, self.igv_ref)
    }
    fn denovo_from_trio(&self, a: &CallArgs, sf_ctr_list: &str, sf_case_control_bam_list: &str) -> Result<String> {
        denovo_from_trio(a, sf_ctr_list, sf_case_control_bam_list, true)
    }
}

/// Options of `gnrt_pipeline_local_v38.py` (and of the CHM13 script, which has the same ones).
#[derive(Parser, Debug, Clone)]
#[command(
    name = "xtea",
    about = "Generate the xTea calling pipelines for Illumina short reads on hg38 (Rust port of gnrt_pipeline_local_v38.py)",
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
    #[arg(long = "hard", help = "This is hard-cut for fitering out coverage abnormal candidates")]
    pub hard: bool,
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

    #[arg(allow_hyphen_values = true, value_parser = parse_int, short = 'n', long = "cores", default_value_t = 1, help = "number of cores")]
    pub cores: i64,
    #[arg(allow_hyphen_values = true, value_parser = parse_int, short = 'm', long = "memory", default_value_t = 20, help = "Memory limit in GB")]
    pub memory: i64,
    #[arg(allow_hyphen_values = true, short = 'q', long = "partition", help = "Which queue to run the job")]
    pub partition: Option<String>,
    #[arg(allow_hyphen_values = true, short = 't', long = "time", default_value = "", help = "Time limit")]
    pub time: String,

    #[arg(allow_hyphen_values = true, short = 'p', long = "path", help = "Working folder")]
    pub wfolder: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'r', long = "ref", help = "reference genome")]
    pub r#ref: Option<String>,
    #[arg(allow_hyphen_values = true, short = 'g', long = "gene", default_value = "gene.aff3", help = "Gene annotation file")]
    pub gene: String,
    #[arg(allow_hyphen_values = true, long = "xtea", help = "xTEA folder")]
    pub xtea: Option<String>,

    #[arg(
        allow_hyphen_values = true,
        value_parser = parse_int,
        short = 'f',
        long = "flag",
        help = "Flag indicates which step to run (1-clip, 2-disc, 16-filter, 256-transduction, 512-post-filter, 1024-annotation, 4096-sibling; 5907 runs all; the 10X steps 4/8/32/64/128 are not supported)"
    )]
    pub flag: Option<i64>,
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
            hard: o.hard,
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
            partition: o.partition,
            time: o.time,
            wfolder: o.wfolder,
            r#ref: o.r#ref,
            gene: Some(o.gene),
            xtea: o.xtea,
            flag: o.flag,
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

/// Entry point of the `xtea` binary; returns the exit status.
pub fn main() -> i32 {
    let opts: Opts = parse_opts("xtea");
    finish(run_local(&HG38, &opts.into()))
}
