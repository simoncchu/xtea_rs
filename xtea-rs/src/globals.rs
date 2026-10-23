//! Port of `global_values.py`: constants and process-wide settings.
//!
//! The mutable settings are stored in atomics so that they can be changed at start-up
//! (mirroring the `global_values.set_*` functions) and read from worker threads.

use std::sync::atomic::{AtomicBool, AtomicI64, AtomicU64, Ordering};

macro_rules! global_i64 {
    ($name:ident, $getter:ident, $setter:ident, $default:expr) => {
        static $name: AtomicI64 = AtomicI64::new($default);
        #[inline]
        pub fn $getter() -> i64 {
            $name.load(Ordering::Relaxed)
        }
        #[inline]
        pub fn $setter(v: i64) {
            $name.store(v, Ordering::Relaxed)
        }
    };
}

macro_rules! global_bool {
    ($name:ident, $getter:ident, $setter:ident, $default:expr) => {
        static $name: AtomicBool = AtomicBool::new($default);
        #[inline]
        pub fn $getter() -> bool {
            $name.load(Ordering::Relaxed)
        }
        #[inline]
        pub fn $setter(v: bool) {
            $name.store(v, Ordering::Relaxed)
        }
    };
}

macro_rules! global_f64 {
    ($name:ident, $getter:ident, $setter:ident, $default:expr) => {
        static $name: AtomicU64 = AtomicU64::new(0);
        #[inline]
        pub fn $getter() -> f64 {
            let bits = $name.load(Ordering::Relaxed);
            if bits == 0 {
                $default
            } else {
                f64::from_bits(bits ^ 1)
            }
        }
        #[inline]
        pub fn $setter(v: f64) {
            // store with a tag bit flipped so that 0 means "unset" (default)
            $name.store(v.to_bits() ^ 1, Ordering::Relaxed)
        }
    };
}

// ---------------------------------------------------------------------------------------------
// originally used in x_TEA_main.py
global_bool!(GLOBAL_RNA_MEDIATED, rna_mediated, set_rna_mediated, true);
pub fn turn_off_rna_mediated() {
    set_rna_mediated(false)
}
global_bool!(GLOBAL_MITCHONDRION_SWITCH, mit_switch_on, set_mit_switch_on, false);
pub fn turn_on_mit() {
    set_mit_switch_on(true)
}

// originally used in clip_read.py
pub const INITIAL_MIN_CLIP_CUTOFF_ILLUMINA: i64 = 2;
pub const INITIAL_MIN_CLIP_CUTOFF_X10: i64 = 3;
global_i64!(INITIAL_MIN_CLIP_CUTOFF, initial_min_clip_cutoff, set_initial_min_clip_cutoff, 2);

pub const BWA_REALIGN_CUTOFF: i64 = 9;
pub const MINIMUM_POLYA_CLIP: i64 = 7;
pub const MAX_CLIP_CLIP_LEN: i64 = 8;
pub const CLIP_PHRED_SCORE_CUTOFF: i64 = 15;

global_i64!(MINIMUM_CLIP_MAPQ, minimum_clip_mapq, set_min_clip_mapq, 12);
pub const CLIP_LOW_MAPQ: i64 = 5;
pub const MAX_LOWQ_CLIP_RATIO: f64 = 0.65;
pub const NEARBY_LOW_Q_CLIP: i64 = 35;

global_i64!(MINIMUM_DISC_MAPQ, minimum_disc_mapq, set_min_disc_mapq, 20);

pub const MAX_BKGRND_LOW_MAPQ: i64 = 10;
pub const MAX_BKGRND_LOW_MAPQ_RATIO: f64 = 0.1;

pub const BWA_PATH: &str = "bwa";
pub const SAMTOOLS_PATH: &str = "samtools";
pub const CLIP_FQ_SUFFIX: &str = ".clipped.fq";
pub const CLIP_BAM_SUFFIX: &str = ".clipped.sam";
pub const CLIP_BAM_SUFFIX2: &str = ".clipped.bam";
pub const CLIP_POS_SUFFIX: &str = ".clip_pos";
pub const CLIP_RE_ALIGN_POS_SUFFIX: &str = ".clip_realign_pos";
pub const LCLIP_FA_SUFFIX: &str = ".lr_clipped.fa";
pub const FLAG_LEFT_CLIP: &str = "L";
pub const FLAG_RIGHT_CLIP: &str = "R";
pub const SEPERATOR: &str = "~";
pub const ALL_DISC_SUFFIX: &str = ".initial.all.disc";

pub const CK_POLYA_CLIP_WIN: i64 = 25;
pub const CK_POLYA_SEQ_MAX: i64 = 20;
pub const POLYA_RATIO: f64 = 0.4;

// originally used in x_alignments.py
pub const BARCODE_COV_CUTOFF: i64 = 600;
pub const MIN_RAW_DISC_CLUSTER_RATIO: f64 = 0.3;
pub const ILLUMINA: &str = "illumina";
pub const X10: &str = "10X";
pub const LONG_READ: &str = "long_read";

// originally used x_TEI_locator.py
pub const MIN_CLIP_FOR_CANDIDATE: i64 = 2;
pub const NEARBY_REGION: i64 = 50;
pub const CLIP_FREQ: i64 = 10;
pub const TRIM_CLIP_FREQ: i64 = 2;
pub const PEAK_WINDOW: i64 = 100;
pub const OUTPUT_BAM_SUFFIX: &str = ".out_bam";
pub const OUTPUT_BAM_HEADER: &str = ".bam_header.sam";
pub const CLIP_FOLDER: &str = "clip";
pub const DISC_FOLDER: &str = "disc";
pub const DISC_SUFFIX: &str = ".discord_pos.txt";
pub const DISC_SUFFIX_FILTER: &str = ".discdt";
pub const RAW_DISC_SUFFIX_FILTER: &str = ".raw.discdt";
pub const CLIP_TMP: &str = "clip_reads_tmp";
pub const DISC_TMP: &str = "discordant_reads_tmp";
pub const RAW_DISC_TMP: &str = "raw_discordant_reads_tmp";
pub const RAW_DISC_TMP_SUFFIX: &str = ".clip_sites_raw_disc.txt";

// originally used in x_clip_disc_filter.py
global_bool!(CHECK_BY_SAMPLE, check_by_sample, set_check_by_sample, false);
pub fn turn_on_check_by_sample() {
    set_check_by_sample(true)
}
global_bool!(IS_CALL_SVA, is_call_sva, set_is_call_sva, false);
pub fn turn_on_sva() {
    set_is_call_sva(true)
}
pub const SVA_ANNOTATION_EXTND: i64 = 200;

pub const DISC_THRESHOLD: i64 = 2000;
pub const TSD_CUTOFF: i64 = 100;
pub const TRANSDCT_UNIQ_MAPQ: i64 = 50;
pub const MINIMAL_TRANSDUCT_MAPQ: i64 = 10;
pub const N_MIN_A_T: i64 = 5;
pub const NEARBY_CLIP: i64 = 50;
pub const CLIP_SEARCH_WINDOW: i64 = 15;
pub const CLIP_CONSISTENT_RATIO: f64 = 0.4;
pub const MAX_COV_TIMES: i64 = 3;
pub const COV_SEARCH_WINDOW: i64 = 1000;
pub const COV_ISD_CHK_WIN: i64 = 900;
pub const LOCAL_COV_WIN: i64 = 200;
pub const MIN_CLIP_MAPPED_RATIO: f64 = 0.65;
pub const MIN_DISC_MAPPED_RATIO: f64 = 0.7;
pub const TWO_SIDE_CLIP_MIN_LEN: i64 = 8;

pub const DISC_NAME_SUFFIX: &str = ".disc_names";
pub const DISC_POS_SUFFIX: &str = ".disc_pos";
pub const GNTP_FEATURE_SUFFIX: &str = ".gntp_features";
pub const ALLELE_FREQUENCY_SUFFIX: &str = ".af";
pub const NOT_TRANSDUCTION: &str = "not_transduction";
pub const ONE_SIDE_FLANKING: &str = "may_two_side_with_one_side_on_flanking";
pub const TWO_SIDE: &str = "two_side";
pub const TWO_SIDE_TPRT_BOTH: &str = "two_side_tprt_both";
pub const TWO_SIDE_TPRT: &str = "two_side_tprt";
pub const ONE_HALF_SIDE: &str = "one_half_side";
pub const ONE_HALF_SIDE_TRPT_BOTH: &str = "one_half_side_tprt_both";
pub const ONE_HALF_SIDE_TRPT: &str = "one_half_side_tprt";
pub const ONE_HALF_SIDE_POLYA_DOMINANT: &str = "polyA_dominant_one-half-side_may_low_confident";
pub const ONE_SIDE: &str = "one_side";
pub const ONE_SIDE_COVERAGE_CONFLICT: &str = "one_side_coverage_conflict";
pub const ONE_SIDE_TRSDCT: &str = "one_side_and_half_transduction";
pub const ONE_SIDE_WEAK: &str = "one_side_weak";
pub const ONE_SIDE_OTHER: &str = "one_side_other";
pub const ONE_SIDE_SV: &str = "one_side_sv";
pub const ONE_SIDE_POLYA_DOMINANT: &str = "polyA_dominant_one-side_may_low_confident";

pub const ORPHAN_TRANSDUCTION: &str = "orphan_or_sibling_transduction";
pub const TWO_SIDE_POLYA_DOMINANT: &str = "polyA_dominant_both-side_may_low_confident";
pub const HIGH_COV_ISD: &str = "high_coverage_island_low_confident";
pub const OTHER_TYPE: &str = "other_type_low_confident";
pub const HIT_END_OF_CNS: &str = "hit_end_of_consensus";
pub const NOT_HIT_END_OF_CNS: &str = "not_hit_end_of_consensus";
pub const BOTH_END_CONSISTNT: &str = "both_end_consistent";
pub const ONE_END_CONSISTNT: &str = "one_end_consistent";
pub const MAX_POLYA_RATIO: f64 = 0.85;
pub const RAW_DISC_FA_SUFFIX: &str = ".raw_disc.fa";
pub const RAW_CLIP_FQ_SUFFIX: &str = ".raw_clip.fq";
pub const HIGH_CONFIDENT_SUFFIX: &str = ".high_confident";
pub const FIVE_PRIME_INVERSION: &str = "5-prime-inversion";
pub const NOT_FIVE_PRIME_INV: &str = "Not-5prime-inversion";
pub const SIBLING_LABEL: &str = "sibling";
pub const ORPHAN_LABEL: &str = "orphan";

// originally used in x_transduction.py
pub const MIN_POLYMORPHIC_SOURCE_DIST: i64 = 1000;
pub const TRANSDCT_MULTI_SOURCE_MIN_RATIO: f64 = 0.45;
pub const TRANSDCT_OTHER_SIDE_CLUSTER_MIN_RATIO: f64 = 0.1;
pub const TRANSDCT_REGION_CLUSTER_MIN_RATIO: f64 = 0.67;
pub const TD_CLIP_QLFD_RATIO: f64 = 0.8;
pub const TD_DECOY_LINE: &str = "LINE1";
pub const TD_DECOY_SVA: &str = "SVA";
pub const TD_DECOY_ALU: &str = "ALU";
pub const TD_DECOY_HERV: &str = "HERV";
pub const S_POLYMORPHIC: &str = "polymerphic";
pub const F_MIN_TRSDCT_DISC_MAP_RATION: f64 = 0.65;
pub const F_MIN_RSC_DISC_MAP_RATION: f64 = 0.6;
pub const TWO_SITES_CLOSE_DIST: i64 = 100;
pub const MIN_HALF_SIBLING_DISC_RATIO: f64 = 0.85;
pub const TD_NON_SIBLING_SUFFIX: &str = ".all_non_sibling_td.txt";
pub const TD_NEW_SITES_SUFFIX: &str = ".new_sites";
// originally used in x_orphan_transdcution.py
pub const MIN_SIBLING_FLANK_LEN: i64 = 15;
pub const MAX_SIBLING_FLANK_LEN: i64 = 5000;
pub const MIN_ORPHAN_DISC_PAIR_INSERT_SIZE: i64 = 100000;
pub const MIN_CLUSTER_RC_RATIO: f64 = 0.6;
pub const MIN_POLYA_CLIP_RATIO: f64 = 0.1;
pub const DOMINANT_POLYA_T_MIN_RATIO: f64 = 0.75;

// originally defined in x_annotation.py, x_reference.py
pub const S_DELIM: &str = "~";
pub const LEFT_FLANK: &str = "left";
pub const RIGHT_FLANK: &str = "right";
global_i64!(LOAD_RMSK_LEFT_EXTND, load_rmsk_left_extnd, set_load_rmsk_left_extnd, 100);

// originally used in x_reads_collection.py
pub const READS_FOLDER: &str = "reads_fa";
pub const HAP1: &str = "hap1";
pub const HAP2: &str = "hap2";
pub const HAP_UNKNOWN: &str = "hap_unknown";
pub const ALL_HAP: &str = "all_hap";
pub const HAP_DISCORD: &str = "hap_discord";
pub const HAP1_SLCT: &str = "hap1_slct";
pub const HAP2_SLCT: &str = "hap2_slct";
pub const HAP_UNKNOWN_SLCT: &str = "hap_unknown_slct";
pub const ALL_HAP_SLCT: &str = "all_hap_slct";
pub const DISC_HAP_SLCT: &str = "disc_slct";
pub const UNIQUE_MAP_CUTOFF: i64 = 60;

// originally used in x_contig.py and x_local_assembly.py
pub const IDBA_UD: &str = "idba_ud";
pub const MINIMAP2: &str = "minimap2";
pub const BWA: &str = "bwa";
pub const CMD_FOLDER: &str = "command";
pub const ASM_FOLDER: &str = "asm_contig";
pub const FLANK_FOLDER: &str = "flank_regions";
pub const MAP_FOLDER: &str = "contig_alignmt";
pub const FILTER_FOLDER: &str = "filter_out";
pub const TEI_SEQ_FOLDER: &str = "tei_seq";
pub const TEI_ON_CNS_CLIP_LENTH: i64 = 100;
pub const MIN_CONTIG_CLIP_LENTH: i64 = 30;

// originally used in x_mutation.py
pub const MUTATION_FOLDER: &str = "mutations";
pub const BCFTOOLS_PATH: &str = "bcftools";
pub const TABIX_PATH: &str = "tabix";

// originally used in x_gene_annotation.py
pub const UP_STREAM_REGION: &str = "up_stream";
pub const DOWN_STREAM_REGION: &str = "down_stream";
pub const UP_DOWN_GENE: i64 = 1500;
pub const NON_GENE: &str = "not_gene_region";

global_bool!(FORCE_CLEAN, force_clean, set_force_clean_value, false);
pub fn set_force_clean() {
    set_force_clean_value(true)
}

// originally used in l_asm.py (long reads; kept for completeness)
pub const WTDBG2: &str = "wtdbg2";
pub const RACON: &str = "racon";
pub const WTPOA: &str = "wtpoa-cns";
pub const FAIL_CLP: &str = "COLLECT_CLIP_FAIL";
pub const FAIL_ASM: &str = "ASM_FAIL";
pub const FAIL_CNS: &str = "POLISH_FAIL";
pub const SUCCD: &str = "succeed";
pub const LRD_CNS_SUFFIX: &str = "_wtdbg2.ctg.lay.fa";
pub const LRD_ALNMT_SUFFIX: &str = "_aln_flank_2_ctg.bam";
pub const LRD_CLIP_SEARCH_WIN: i64 = 20000;
pub const LRD_MIN_CLIP_LTH: i64 = 50;
pub const BRKPNT_CK_WIN: i64 = 500;
pub const LRD_MIN_INS_LTH: i64 = 40;
pub const LRD_MAX_BRKPNT_CHK_REGION: i64 = 300;
pub const LRD_MIN_CLIP_COLLECT_LTH: i64 = 200;
pub const LRD_MIN_CTN_COLLECT_LTH: i64 = 50;
global_i64!(LRD_EXTND_LEN, lrd_extnd_len, set_lrd_extnd_len, 2000);

// originally used in x_coverage.py
pub const N_RANDOM_SITES: i64 = 3000;
pub const MIN_COV_RANDOM_SITE: f64 = 5.0;

global_f64!(AVE_COVERAGE, average_cov, set_average_cov, 30.0);
// whether AVE_COVERAGE holds a Python int (the default 30, or the (0, 0, 0, 0) fallback)
global_bool!(AVE_COVERAGE_IS_INT, average_cov_is_int, set_average_cov_is_int, true);
/// `str(global_values.AVE_COVERAGE)`
pub fn average_cov_str() -> String {
    if average_cov_is_int() {
        format!("{}", average_cov() as i64)
    } else {
        crate::py::fmt::float_repr(average_cov())
    }
}

// x_genotype.py
pub const BWA_HALF_READ_MIN_SCORE: i64 = 45;
pub const CLIP_EXACT_CLIP_SLACK: i64 = 3;
pub const LARGE_INDEL_IN_READ: i64 = 3;
global_f64!(DFT_IS, dft_is, set_insert_size, 550.0);
global_f64!(READ_LENGTH, read_length, set_read_length, 100.0);

// originally used in x_basic_info.py
pub const BASIC_INFO_FILE: &str = "basic_cov_is_rlth_info.txt";
pub const MAX_NORMAL_INSERT_SIZE: i64 = 2000;

// originally used in x_insert_size.py
pub const ISIZE_MEAN_SUFFIX: &str = ".mean_temp_is";

// originally used in x_tprt_filter.py
pub const ONE_SIDE_POLYA_CUTOFF: f64 = 0.75;

pub const LRD_MIN_MAPQ: i64 = 20;
pub const LRD_MIN_MAP_LEN: i64 = 1000;
pub const LRD_MIN_FOCUS_INS_LEN: i64 = 50;
pub const LRD_BRKPNT_FOCAL_REGIN: i64 = 75;
pub const LRD_BRKPNT_FOCAL_CLIP_RATIO: i64 = 1;
pub const LRD_BRKPNT_MAXIMUM_STD: i64 = 250;
pub const LRD_POLYMORPHIC_MIN_CLIP_LEN: i64 = 1500;
global_i64!(LRD_POLYM_BRK_CHK_WIN, lrd_polym_brk_chk_win, set_polymorphic_brk_chk_win, 250);
pub const LRD_PRI_SUP_MAX_OVRLAP: i64 = 100;
pub const LRD_PRI_SUP_MAX_GAP_IN_SEQ: i64 = 75;
pub const LRD_PRI_SUP_DOMINANT_POLYA_RATIO: f64 = 0.45;
pub const LRD_PRI_SUP_MAX_GAP_MAP_RATIO: f64 = 0.25;
pub const LRD_PRI_SUP_MAX_CLIP: i64 = 45;
pub const LRD_MIN_INTERNAL_DEL: i64 = 50;
pub const LRD_TRSDCT_MIN_MAP_RATION: f64 = 0.8;
pub const LRD_MAX_TSD_LEN: i64 = 80;
pub const LRD_MIN_TSD_LEN: i64 = 4;
pub const LRD_MIN_TSD_MATCH_RATIO: f64 = 0.85;
pub const LRD_PSEUDOGENE_INPUT_SUF: &str = ".for_pseudogene.fa";

// originally used in x_post_filter.py
global_i64!(TWO_CLIP_CLUSTER_DIFF_CUTOFF, two_clip_cluster_diff_cutoff, set_two_clip_cluster_diff_cutoff, 300);

pub const REP_DIVERGENT_CUTOFF: i64 = 15;
pub const REP_LOW_DIVERGENT_CUTOFF: i64 = 7;
pub const TD_REP_DIVERGENT_CUTOFF: i64 = 5;
pub const ALPHA_SAT_RMSK: &str = "ALR/Alpha";
pub const HSATII_SAT_RMSK: &str = "HSATII";
