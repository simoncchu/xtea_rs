//! Port of `x_rep_type.py`.

pub const REP_TYPE_L1: &str = "LINE1";
pub const REP_TYPE_ALU: &str = "ALU";
pub const REP_TYPE_SVA: &str = "SVA";
pub const REP_TYPE_HERV: &str = "HERV";
pub const REP_TYPE_MIT: &str = "Mitochondria";
pub const REP_TYPE_MSTA: &str = "Msta";
pub const PSUDOGENE: &str = "Pseudogene";

pub const REP_TYPE_L1_VCF: &str = "INS:ME:LINE1";
pub const REP_TYPE_ALU_VCF: &str = "INS:ME:ALU";
pub const REP_TYPE_SVA_VCF: &str = "INS:ME:SVA";
pub const REP_TYPE_HERV_VCF: &str = "INS:ME:HERV-K";
pub const REP_TYPE_MIT_VCF: &str = "INS:MT";
pub const REP_TYPE_MSTA_VCF: &str = "INS:Msta";
pub const PSUDOGENE_VCF: &str = "INS:PSDGN";

const TYPES: [(i64, &str, &str); 7] = [
    (1, REP_TYPE_L1, REP_TYPE_L1_VCF),
    (2, REP_TYPE_ALU, REP_TYPE_ALU_VCF),
    (4, REP_TYPE_SVA, REP_TYPE_SVA_VCF),
    (8, REP_TYPE_HERV, REP_TYPE_HERV_VCF),
    (16, REP_TYPE_MIT, REP_TYPE_MIT_VCF),
    (32, REP_TYPE_MSTA, REP_TYPE_MSTA_VCF),
    (64, PSUDOGENE, PSUDOGENE_VCF),
];

/// `get_rep_type`
pub fn get_rep_type(i_rep_type: i64) -> &'static str {
    for (bit, name, _) in TYPES.iter() {
        if i_rep_type & bit != 0 {
            return name;
        }
    }
    "unkown"
}

/// `get_rep_vcf_alt_type`
pub fn get_rep_vcf_alt_type(i_rep_type: i64) -> &'static str {
    for (bit, _, vcf) in TYPES.iter() {
        if i_rep_type & bit != 0 {
            return vcf;
        }
    }
    "unkown"
}

/// `get_all_rep_types` / `parse_out_all_rep_type` (ordered list)
pub fn parse_out_all_rep_type(i_rep_type: i64) -> Vec<&'static str> {
    TYPES.iter().filter(|(bit, _, _)| i_rep_type & bit != 0).map(|(_, n, _)| *n).collect()
}
