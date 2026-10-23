//! `xtea_chm13`: port of `gnrt_pipeline_local_chm13.py`, the CHM13 copy of the hg38 generator
//! ([`super::local_v38`]).  It has the same options and commands except:
//! * the library files (`<type>/chm13/chm13_*`; Pseudogene keeps the GRCh38 file names in a
//!   `Pseudogene/chm13/` folder),
//! * `--ref "chm13"` in the igv steps of the case-control mode,
//! * the `--sr_asm` step (flag 8192) reads `candidate_disc_filtered_cns_post_filtering_denovo.txt`
//!   for non-tumor samples,
//! * the `-V` text.

use super::local_v38::{
    Opts, V38Family, REP_TYPE_ALU, REP_TYPE_HERV, REP_TYPE_L1, REP_TYPE_MIT, REP_TYPE_MSTA, REP_TYPE_PSEUDOGENE, REP_TYPE_SVA,
};
use super::{finish, parse_opts, run_local, RepType};

/// Repeat types and library files of `gnrt_lib_config` (CHM13), in `-y` bit order.
pub static CHM13_REP_TYPES: [RepType; 7] = [
    RepType {
        name: REP_TYPE_L1,
        bit: 1,
        anno: "LINE/chm13/chm13_L1_larger_500_with_all_L1HS.out",
        anno1: "LINE/chm13/chm13_L1.fa.out",
        copy_with_flank: "LINE/chm13/chm13_L1HS_copies_larger_5K_with_flank.fa",
        flank: Some("LINE/chm13/chm13_FL_L1_flanks_3k.fa"),
        cns: "consensus/LINE1.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_ALU,
        bit: 2,
        anno: "Alu/chm13/chm13_Alu.out",
        anno1: "Alu/chm13/chm13_Alu.out",
        copy_with_flank: "Alu/chm13/chm13_AluJabc_copies_with_flank.fa",
        flank: None,
        cns: "consensus/ALU.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_SVA,
        bit: 4,
        anno: "SVA/chm13/chm13_SVA.out",
        anno1: "SVA/chm13/chm13_SVA.out",
        copy_with_flank: "SVA/chm13/chm13_SVA_copies_with_flank.fa",
        flank: Some("SVA/chm13/chm13_FL_SVA_flanks_3k.fa"),
        cns: "consensus/SVA.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_HERV,
        bit: 8,
        anno: "HERV/chm13/chm13_HERV.out",
        anno1: "HERV/chm13/chm13_HERV.out",
        copy_with_flank: "HERV/chm13/chm13_HERV_copies_with_flank.fa",
        flank: None,
        cns: "consensus/HERV.fa",
        ref_before_anno1: true,
    },
    RepType {
        name: REP_TYPE_MIT,
        bit: 16,
        anno: "Mitochondrion/chm13/chm13_numtS.out",
        anno1: "Mitochondrion/chm13/chm13_numtS.out",
        copy_with_flank: "Mitochondrion/chm13/chm13_mitochondrion_copies_with_flank.fa",
        flank: None,
        cns: "consensus/mitochondrion.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_MSTA,
        bit: 32,
        anno: "MSTA/chm13/chm13_MSTA.out",
        anno1: "MSTA/chm13/chm13_MSTA.out",
        copy_with_flank: "MSTA/chm13/chm13_MSTA_copies_with_flank.fa",
        flank: None,
        cns: "consensus/MSTA.fa",
        ref_before_anno1: false,
    },
    RepType {
        name: REP_TYPE_PSEUDOGENE,
        bit: 64,
        anno: "Pseudogene/chm13/gencode_v28_GRCh38_exon_annotation.out",
        anno1: "Pseudogene/chm13/gencode_v28_GRCh38_exon_annotation.out",
        copy_with_flank: "Pseudogene/chm13/gencode_v28_GRCh38_transcript_masked.fa",
        flank: None,
        cns: "consensus/gencode_v28_GRCh38_transcript_masked.fa",
        ref_before_anno1: false,
    },
];

pub static CHM13: V38Family = V38Family {
    genome: "chm13",
    igv_ref: "chm13",
    sr_asm_input: "candidate_disc_filtered_cns_post_filtering_denovo.txt",
    rep_types: &CHM13_REP_TYPES,
};

/// Entry point of the `xtea_chm13` binary; returns the exit status.
pub fn main() -> i32 {
    let opts: Opts = parse_opts("xtea_chm13");
    finish(run_local(&CHM13, &opts.into()))
}
