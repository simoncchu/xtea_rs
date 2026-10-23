//! Conversion of the final xTea calls to VCF (port of `x_gvcf.py`, `x_TSD.py` and
//! `XReference.get_ref_seqs_of_sites`).

use crate::annotation::process_chrm_name;
use crate::bam::AlignmentFile;
use crate::fasta::FastaFile;
use crate::py::fmt::{float_repr, parse_py_float, parse_py_int};
use crate::rep_type;
use crate::utils::{self, split_ws};
use anyhow::{anyhow, bail, Result};
use std::collections::HashSet;
use std::io::Write;

fn pint(s: &str) -> Result<i64> {
    parse_py_int(s).ok_or_else(|| anyhow!("invalid literal for int() with base 10: '{}'", s))
}

fn pfloat(s: &str) -> Result<f64> {
    parse_py_float(s).ok_or_else(|| anyhow!("could not convert string to float: '{}'", s))
}

/// `XTSD.get_TSD_seqs`: "+seq" for a target site duplication (left clip < right clip),
/// "-seq" for a target site deletion, "NULL" when a position is unknown.
pub fn get_tsd_seqs(sf_ref: &str, l_sites: &[(String, i64, i64)]) -> Result<Vec<String>> {
    let f_fa = FastaFile::open(sf_ref)?;
    let b_with_chr = f_fa.has_reference("chr1");
    let mut l_tsd = Vec::with_capacity(l_sites.len());
    for (ins_chrm, lpos, rpos) in l_sites {
        let (lpos, rpos) = (*lpos, *rpos);
        if lpos <= 0 || rpos <= 0 {
            l_tsd.push("NULL".to_string());
            continue;
        }
        let ref_chrm = process_chrm_name(ins_chrm, b_with_chr);
        let (mut istart, mut iend, mut s_signal) = (lpos, rpos, "+");
        if lpos > rpos {
            istart = rpos;
            iend = lpos;
            s_signal = "-";
        }
        let s_tsd = f_fa.fetch(&ref_chrm, istart, iend)?;
        l_tsd.push(format!("{}{}", s_signal, s_tsd));
    }
    Ok(l_tsd)
}

/// `XReference.get_ref_seqs_of_sites`: reference sequence of each [start, end) site ("." when
/// unknown or when the chromosome is missing from the reference).
pub fn get_ref_seqs_of_sites(sf_ref: &str, l_sites: &[(String, i64, i64)]) -> Result<Vec<String>> {
    let f_fa = FastaFile::open(sf_ref)?;
    let refs: HashSet<&str> = f_fa.references().iter().map(|s| s.as_str()).collect();
    let b_with_chr = refs.contains("chr1");
    let mut l_seqs = Vec::with_capacity(l_sites.len());
    for (ins_chrm, lpos, rpos) in l_sites {
        let (lpos, rpos) = (*lpos, *rpos);
        if lpos <= 0 || rpos <= 0 {
            l_seqs.push(".".to_string());
            continue;
        }
        let ref_chrm = process_chrm_name(ins_chrm, b_with_chr);
        let (istart, iend) = if lpos > rpos { (rpos, lpos) } else { (lpos, rpos) };
        if !refs.contains(ref_chrm.as_str()) {
            println!("{} doesn't exist in reference!", ref_chrm);
            l_seqs.push(".".to_string());
            continue;
        }
        l_seqs.push(f_fa.fetch(&ref_chrm, istart, iend)?);
    }
    Ok(l_seqs)
}

/// `gVCF`
#[derive(Default)]
pub struct GVcf;

const HEAD_OTHER: &str = "##ALT=<ID=INS:ME:SVA,Description=\"Insertion of SVA element\">\n\
##ALT=<ID=INS:ME:LINE1,Description=\"Insertion of LINE1 element\">\n\
##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">\n\
##ALT=<ID=INS:ME:HERV-K,Description=\"Insertion of HERV-K element\">\n\
##ALT=<ID=INS:MT,Description=\"Mitochondrial insertion\">\n\
##ALT=<ID=INS:PSDGN,Description=\"Pseudogene insertion\">\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FILTER=<ID=LowConfident,Description=\"All filters passed\">\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n\
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Insertion length\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of this variant\">\n\
##INFO=<ID=TSD,Number=1,Type=String,Description=\"Precise Target Site Duplication for bases, if unknown, value will be NULL\">\n\
##INFO=<ID=TSDLEN,Number=1,Type=Integer,Description=\"TSD (dup/del) length, if unknown, then -1\">\n\
##INFO=<ID=SUBTYPE,Number=1,Type=String,Description=\"Subtype based on TSD/PolyA/side-info\">\n\
##INFO=<ID=TD_SRC,Number=1,Type=String,Description=\"Transduction source\">\n\
##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Insertion orientation (+/- for sense/antisense)\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1)\">\n\
##INFO=<ID=LCLIP,Number=1,Type=Integer,Description=\"Number of left clipped reads (aligned to copies)\">\n\
##INFO=<ID=RCLIP,Number=1,Type=Integer,Description=\"Number of right clipped reads (aligned to copies)\">\n\
##INFO=<ID=LDISC,Number=1,Type=Integer,Description=\"Number of left discordant reads (aligned to copies)\">\n\
##INFO=<ID=RDISC,Number=1,Type=Integer,Description=\"Number of right discordant reads (aligned to copies)\">\n\
##INFO=<ID=LPOLYA,Number=1,Type=Integer,Description=\"Number of left polyA reads\">\n\
##INFO=<ID=RPOLYA,Number=1,Type=Integer,Description=\"Number of right polyA reads\">\n\
##INFO=<ID=LRAWCLIP,Number=1,Type=Integer,Description=\"Number of left raw clip reads\">\n\
##INFO=<ID=RRAWCLIP,Number=1,Type=Integer,Description=\"Number of right raw clip reads\">\n\
##INFO=<ID=AF_CLIP,Number=1,Type=Integer,Description=\"Number of effective clip reads (for AF)\">\n\
##INFO=<ID=AF_FMAP,Number=1,Type=Integer,Description=\"Number of effective fully mapped reads (for AF)\">\n\
##INFO=<ID=AF_DISC,Number=1,Type=Integer,Description=\"Number of effective discordant pairs (for AF)\">\n\
##INFO=<ID=AF_CONCORDNT,Number=1,Type=Integer,Description=\"Number of effective concordant pairs (for AF)\">\n\
##INFO=<ID=LDRC,Number=1,Type=Integer,Description=\"Number of left discordant reads are reverse-complementary\">\n\
##INFO=<ID=LDNRC,Number=1,Type=Integer,Description=\"Number of left discordant reads are not reverse-complementary\">\n\
##INFO=<ID=RDRC,Number=1,Type=Integer,Description=\"Number of right discordant reads are reverse-complementary\">\n\
##INFO=<ID=RDNRC,Number=1,Type=Integer,Description=\"Number of right discordant reads are not reverse-complementary\">\n\
##INFO=<ID=LCOV,Number=1,Type=Float,Description=\"Left focal coverage\">\n\
##INFO=<ID=RCOV,Number=1,Type=Float,Description=\"Right focal coverage\">\n\
##INFO=<ID=LD_AKR_RC,Number=1,Type=Integer,Description=\"Number of left discordant anchor reads are reverse-complementary\">\n\
##INFO=<ID=LD_AKR_NRC,Number=1,Type=Integer,Description=\"Number of left discordant anchor reads are not reverse-complementary\">\n\
##INFO=<ID=RD_AKR_RC,Number=1,Type=Integer,Description=\"Number of right discordant anchor reads are reverse-complementary\">\n\
##INFO=<ID=RD_AKR_NRC,Number=1,Type=Integer,Description=\"Number of right discordant anchor reads are not reverse-complementary\">\n\
##INFO=<ID=LC_CLUSTER,Number=1,Type=String,Description=\"Left clip formed cluster in consensus\">\n\
##INFO=<ID=RC_CLUSTER,Number=1,Type=String,Description=\"Right clip formed cluster in consensus\">\n\
##INFO=<ID=LD_CLUSTER,Number=1,Type=String,Description=\"Left discordant formed cluster in consensus\">\n\
##INFO=<ID=RD_CLUSTER,Number=1,Type=String,Description=\"Right discordant formed cluster in consensus\">\n\
##INFO=<ID=NINDEL,Number=1,Type=Integer,Description=\"Number of large indel (or with several mismatches) reads\">\n\
##INFO=<ID=CLIP_LEN,Number=1,Type=String,Description=\"Clip reads length\">\n\
##INFO=<ID=INS_INV,Number=1,Type=String,Description=\"5 prime inversion\">\n\
##INFO=<ID=REF_REP,Number=1,Type=String,Description=\"Fall in reference repeat copy or not\">\n\
##INFO=<ID=GENE_INFO,Number=1,Type=String,Description=\"Fall in gene region\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";

const HG19: [(&str, i64); 24] = [
    ("1", 249250621), ("2", 243199373), ("3", 198022430), ("4", 191154276), ("5", 180915260),
    ("6", 171115067), ("7", 159138663), ("8", 146364022), ("9", 141213431), ("10", 135534747),
    ("11", 135006516), ("12", 133851895), ("13", 115169878), ("14", 107349540), ("15", 102531392),
    ("16", 90354753), ("17", 81195210), ("18", 78077248), ("19", 59128983), ("20", 63025520),
    ("21", 48129895), ("22", 51304566), ("X", 155270560), ("Y", 59373566),
];

const HG38: [(&str, i64); 24] = [
    ("1", 248956422), ("2", 242193529), ("3", 198295559), ("4", 190214555), ("5", 181538259),
    ("6", 170805979), ("7", 159345973), ("8", 145138636), ("9", 138394717), ("10", 133797422),
    ("11", 135086622), ("12", 133275309), ("13", 114364328), ("14", 107043718), ("15", 101991189),
    ("16", 90338345), ("17", 83257441), ("18", 80373285), ("19", 58617616), ("20", 64444167),
    ("21", 46709983), ("22", 50818468), ("X", 156040895), ("Y", 57227415),
];

fn contig_lines(tbl: &[(&str, i64)], b_chr: bool) -> String {
    let s_chr = if b_chr { "chr" } else { "" };
    tbl.iter().map(|(c, l)| format!("##contig=<ID={}{},length={}>\n", s_chr, c, l)).collect()
}

impl GVcf {
    pub fn new() -> Self {
        GVcf
    }

    /// `prepare_head_hg19`
    pub fn prepare_head_hg19(&self, b_chr: bool) -> String {
        contig_lines(&HG19, b_chr)
    }

    /// `prepare_head_hg38`
    pub fn prepare_head_hg38(&self, b_chr: bool) -> String {
        contig_lines(&HG38, b_chr)
    }

    /// `prepare_head_other`
    pub fn prepare_head_other(&self) -> String {
        HEAD_OTHER.to_string()
    }

    /// `gnrt_vcf_head`: VCF header with the contigs of the (first) BAM.
    pub fn gnrt_vcf_head(&self, s_sample_id: &str, sf_bam: &str) -> Result<String> {
        let samfile = AlignmentFile::open(sf_bam, "")?;
        let mut shead = String::from("##fileformat=VCFv4.2\n");
        match samfile.header_sq_records() {
            None => {
                println!("No head information in bam file: {}", sf_bam);
                // the Python code calls an undefined module-level `prepare_head_hg38` here
                bail!("name 'prepare_head_hg38' is not defined");
            }
            Some(recs) => {
                let mut sf_ref = String::new();
                let mut l_chrm = Vec::new();
                for rcd in recs.iter() {
                    let get = |k: &str| rcd.iter().find(|(tk, _)| tk == k).map(|(_, v)| v.clone());
                    sf_ref = get("UR").unwrap_or_else(|| "NULL".to_string());
                    let i_chrm_len = get("LN").unwrap_or_else(|| "-1".to_string());
                    let s_chrm = get("SN").unwrap_or_else(|| "NULL".to_string());
                    l_chrm.push((s_chrm, i_chrm_len));
                }
                let last = sf_ref.rsplit('/').next().unwrap_or("");
                shead.push_str(&format!("##reference={}\n", last));
                for (c, l) in l_chrm.iter() {
                    if c.contains("HLA") || c.contains("Un") || c.contains('_') {
                        continue;
                    }
                    shead.push_str(&format!("##contig=<ID={},length={}>\n", c, l));
                }
            }
        }
        shead.push_str(&self.prepare_head_other());
        shead.push_str(&format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n", s_sample_id));
        Ok(shead)
    }

    /// `sort_raw_rslt`
    pub fn sort_raw_rslt(&self, sf_raw_rslt: &str, sf_sorted_rslt: &str) {
        let cmd = format!("sort -k1,1V -k2,2n -o {} {}", sf_sorted_rslt, sf_raw_rslt);
        utils::run_cmd_small_output(&cmd);
    }

    /// `cvt_raw_rslt_to_gvcf`
    pub fn cvt_raw_rslt_to_gvcf(&self, s_sample_id: &str, sf_bam: &str, sf_raw_rslt: &str, i_rep_type: i64, sf_ref: &str, sf_vcf: &str) -> Result<()> {
        let sf_sorted = format!("{}.sorted", sf_raw_rslt);
        self.sort_raw_rslt(sf_raw_rslt, &sf_sorted);
        let mut l_rslts: Vec<Vec<String>> = Vec::new();
        let mut l_sites = Vec::new();
        let mut l_pos = Vec::new();
        for line in utils::read_raw_lines(&sf_sorted)? {
            let fields: Vec<String> = split_ws(&line).into_iter().map(|s| s.to_string()).collect();
            let g = |i: usize| -> Result<&str> {
                fields.get(i).map(|s| s.as_str()).ok_or_else(|| anyhow!("list index out of range"))
            };
            let lclip_pos = pint(g(2)?)?;
            let rclip_pos = pint(g(3)?)?;
            l_sites.push((g(0)?.to_string(), lclip_pos, rclip_pos));
            let (mut i_begin, i_end) = if lclip_pos > rclip_pos { (rclip_pos, lclip_pos) } else { (lclip_pos, rclip_pos) };
            if i_begin <= 0 {
                i_begin = i_end;
            }
            l_pos.push((g(0)?.to_string(), i_begin, i_begin + 1));
            l_rslts.push(fields);
        }
        let s_rep_type = rep_type::get_rep_vcf_alt_type(i_rep_type);
        let l_stsd = get_tsd_seqs(sf_ref, &l_sites)?;
        let l_ref_seqs = get_ref_seqs_of_sites(sf_ref, &l_pos)?;
        let mut w = utils::create(sf_vcf)?;
        let s_head = self.gnrt_vcf_head(s_sample_id, sf_bam)?;
        w.write_all(s_head.as_bytes())?;
        if l_rslts.len() != l_sites.len() || l_sites.len() != l_ref_seqs.len() || l_rslts.len() != l_ref_seqs.len() {
            println!("[Error in x_vcf.py]: Arrzy size is different!\n");
        }
        for ((l_fields, s_tsd), s_site_seq) in l_rslts.iter().zip(l_stsd.iter()).zip(l_ref_seqs.iter()) {
            let s_rcd = self.cvt_to_vcf_rcd(l_fields, s_rep_type, s_tsd, ".", s_site_seq)?;
            w.write_all(s_rcd.as_bytes())?;
        }
        w.flush()?;
        Ok(())
    }

    /// `_cvt_to_vcf_rcd`
    fn cvt_to_vcf_rcd(&self, l_fields: &[String], s_type: &str, s_tsd: &str, s_id: &str, s_ref_seq: &str) -> Result<String> {
        let g = |i: usize| -> Result<&str> {
            l_fields.get(i).map(|s| s.as_str()).ok_or_else(|| anyhow!("list index out of range"))
        };
        let gi = |i: usize| -> Result<i64> { pint(g(i)?) };
        let ins_chrm = g(0)?;
        let lclip_pos = gi(2)?;
        let rclip_pos = gi(3)?;
        let i_tsd_len = gi(4)?;
        let n_lclip_cns = gi(5)?;
        let n_rclip_cns = gi(6)?;
        let n_ldisc_cns = gi(7)?;
        let n_rdisc_cns = gi(8)?;
        let n_lpolya = gi(9)?;
        let n_rpolya = gi(10)?;
        let f_lcov = pfloat(g(11)?)?;
        let f_rcov = pfloat(g(12)?)?;
        let slc_cluster = g(19)?;
        let src_cluster = g(20)?;
        let sld_cluster = g(21)?;
        let srd_cluster = g(22)?;
        let s_transduction = g(23)?;
        let n_lpe_same = gi(27)?;
        let n_lpe_diff = gi(28)?;
        let n_rpe_same = gi(29)?;
        let n_rpe_diff = gi(30)?;
        let s_inv = g(31)?;
        let s_sub_type = g(32)?;
        let n_af_clip = gi(35)?;
        let n_af_full = gi(36)?;
        // (sic) the condition tests the clip count twice
        let s_af = if n_af_clip + n_af_clip > 0 {
            let den = n_af_clip + n_af_full;
            if den == 0 {
                bail!("float division by zero");
            }
            float_repr(n_af_clip as f64 / den as f64)
        } else {
            "0".to_string()
        };
        let n_raw_lclip = gi(37)?;
        let n_raw_rclip = gi(38)?;
        let n_af_disc = gi(39)?;
        let n_af_concod = gi(40)?;
        let n_indel = gi(41)?;
        let s_clip_lens = g(42)?;
        let n_lanchor_rc = gi(43)?;
        let n_lanchor_nrc = gi(44)?;
        let n_ranchor_rc = gi(45)?;
        let n_ranchor_nrc = gi(46)?;
        let i_ins_len = gi(47)?;
        let s_fall_in_ref_copy = g(48)?;
        let s_gene_info = if l_fields.len() > 49 { g(49)? } else { "." };
        let s_gntp = if l_fields.len() > 50 { g(50)? } else { "./." };
        let s_strand = if n_lpolya > 0 && n_rpolya == 0 {
            "+"
        } else if n_lpolya == 0 && n_rpolya > 0 {
            "-"
        } else {
            "."
        };
        let (mut i_begin, i_end) = if lclip_pos > rclip_pos { (rclip_pos, lclip_pos) } else { (lclip_pos, rclip_pos) };
        if i_begin <= 0 {
            i_begin = i_end;
        }
        let mut s = String::new();
        s.push_str(&format!("{}\t{}\t{}\t{}\t<{}>\t.\tPASS\t", ins_chrm, i_begin, s_id, s_ref_seq, s_type));
        s.push_str(&format!("SVTYPE={};", s_type));
        s.push_str(&format!("SVLEN={};", i_ins_len.abs()));
        s.push_str(&format!("END={};", i_end));
        s.push_str(&format!("TSD={};", s_tsd));
        s.push_str(&format!("TSDLEN={};", i_tsd_len));
        s.push_str(&format!("SUBTYPE={};", s_sub_type));
        s.push_str(&format!("TD_SRC={};", s_transduction));
        s.push_str(&format!("STRAND={};", s_strand));
        s.push_str(&format!("AF={};", s_af));
        s.push_str(&format!("LCLIP={};", n_lclip_cns));
        s.push_str(&format!("RCLIP={};", n_rclip_cns));
        s.push_str(&format!("LDISC={};", n_ldisc_cns));
        s.push_str(&format!("RDISC={};", n_rdisc_cns));
        s.push_str(&format!("LPOLYA={};", n_lpolya));
        s.push_str(&format!("RPOLYA={};", n_rpolya));
        s.push_str(&format!("LRAWCLIP={};", n_raw_lclip));
        s.push_str(&format!("RRAWCLIP={};", n_raw_rclip));
        s.push_str(&format!("AF_CLIP={};", n_af_clip));
        s.push_str(&format!("AF_FMAP={};", n_af_full));
        s.push_str(&format!("AF_DISC={};", n_af_disc));
        s.push_str(&format!("AF_CONCORDNT={};", n_af_concod));
        s.push_str(&format!("LDRC={};", n_lpe_same));
        s.push_str(&format!("LDNRC={};", n_lpe_diff));
        s.push_str(&format!("RDRC={};", n_rpe_same));
        s.push_str(&format!("RDNRC={};", n_rpe_diff));
        s.push_str(&format!("LCOV={};", float_repr(f_lcov)));
        s.push_str(&format!("RCOV={};", float_repr(f_rcov)));
        s.push_str(&format!("LD_AKR_RC={};", n_lanchor_rc));
        s.push_str(&format!("LD_AKR_NRC={};", n_lanchor_nrc));
        s.push_str(&format!("RD_AKR_RC={};", n_ranchor_rc));
        s.push_str(&format!("RD_AKR_NRC={};", n_ranchor_nrc));
        s.push_str(&format!("LC_CLUSTER={};", slc_cluster));
        s.push_str(&format!("RC_CLUSTER={};", src_cluster));
        s.push_str(&format!("LD_CLUSTER={};", sld_cluster));
        s.push_str(&format!("RD_CLUSTER={};", srd_cluster));
        s.push_str(&format!("NINDEL={};", n_indel));
        s.push_str(&format!("CLIP_LEN={};", s_clip_lens));
        s.push_str(&format!("INS_INV={};", s_inv));
        s.push_str(&format!("REF_REP={};", s_fall_in_ref_copy));
        s.push_str(&format!("GENE_INFO={}\t", s_gene_info));
        s.push_str(&format!("GT\t{}\n", s_gntp));
        Ok(s)
    }
}

/// The `--gVCF` step of `x_TEA_main.py`: sample id from the first BAM of the list, output
/// `<prefix>/<sample>_<type>.vcf`.
/// `fixed_sample_id`: the cloud entry point's `--id` (otherwise the id is the BAM file name
/// without its extension).
pub fn run_gvcf_step(sf_ref: &str, sf_bam_list: &str, sf_raw_rslt: &str, i_rep_type: i64, sf_prefix: &str, fixed_sample_id: Option<&str>) -> Result<()> {
    let mut sf_bam = String::new();
    let mut s_sample_id = fixed_sample_id.unwrap_or("null").to_string();
    for line in utils::read_raw_lines(sf_bam_list)? {
        let fields = split_ws(&line);
        if fields.is_empty() {
            break;
        }
        sf_bam = fields[0].to_string();
        if fixed_sample_id.is_none() {
            let s_sample_id2 = sf_bam.rsplit('/').next().unwrap_or("").to_string();
            let tmp: Vec<&str> = s_sample_id2.split('.').collect();
            s_sample_id = tmp[..tmp.len() - 1].join(".");
        }
        break;
    }
    if utils::is_file(&sf_bam) {
        let gvcf = GVcf::new();
        let mut prefix = sf_prefix.to_string();
        if !prefix.ends_with('/') {
            prefix.push('/');
        }
        let s_rep_type = rep_type::get_rep_type(i_rep_type);
        let sf_vcf = format!("{}{}_{}.vcf", prefix, s_sample_id, s_rep_type);
        gvcf.cvt_raw_rslt_to_gvcf(&s_sample_id, &sf_bam, sf_raw_rslt, i_rep_type, sf_ref, &sf_vcf)?;
    } else {
        println!("Wrong bam file: {}", sf_bam);
    }
    Ok(())
}
