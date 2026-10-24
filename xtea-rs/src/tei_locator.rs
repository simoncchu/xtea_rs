//! Port of `x_TEI_locator.py` (TE_Multi_Locator and TELocator).

use crate::alignments::{cnt_discordant_pairs, BamInfo};
use crate::annotation::XAnnotation;
use crate::bam::AlignmentFile;
use crate::bwa_align::BwAlign;
use crate::chrom::is_decoy_contig_chrms;
use crate::clip_read::ClipReadInfo;
use crate::globals as gv;
use crate::intermediate_sites::{self as xis, SiteMap};
use crate::num::PyVal;
use crate::utils::{self, split_ws};
use anyhow::{anyhow, Result};
use indexmap::IndexMap;
use std::io::Write;
use std::sync::Arc;

/// The Rust version calls insertions from Illumina short reads only: a BAM list entry with the
/// 10X platform tag (the Python version lowers the initial clip cutoff for those) is rejected
/// before any BAM is processed.
fn reject_x10_bams(sf_list: &str, l_lines: &[String]) -> Result<()> {
    for line in l_lines {
        let fields = split_ws(line);
        if fields.len() >= 2 && utils::rstrip(fields[1]) == gv::X10 {
            return Err(anyhow!(
                "{} is listed as a 10X linked-read bam in {}; the Rust version of xTea supports Illumina short reads only (use the Python xTea for 10X data)",
                fields[0],
                sf_list
            ));
        }
    }
    Ok(())
}

pub struct TeMultiLocator {
    pub sf_list: String,
    pub working_folder: String,
    pub n_jobs: usize,
    pub sf_ref: String,
}

impl TeMultiLocator {
    pub fn new(sf_list: &str, working_folder: &str, n_jobs: usize, sf_ref: &str) -> Self {
        TeMultiLocator { sf_list: sf_list.to_string(), working_folder: working_folder.to_string(), n_jobs, sf_ref: sf_ref.to_string() }
    }

    /// `collect_all_clipped_from_multiple_alignmts`
    pub fn collect_all_clipped_from_multiple_alignmts(&self, sf_annotation: &str, b_se: bool, s_clip_wfolder: &str, wfolder_pub_clip: &str) -> Result<()> {
        for line in utils::read_lines(&self.sf_list)? {
            let fields = split_ws(&line);
            if fields.is_empty() {
                continue;
            }
            let sf_ori_bam = fields[0];
            if sf_ori_bam.len() <= 1 {
                continue;
            }
            let caller = TeLocator::new(sf_ori_bam, &self.working_folder, self.n_jobs, &self.sf_ref);
            caller.collect_all_clipped_reads_only(sf_annotation, b_se, s_clip_wfolder, wfolder_pub_clip)?;
        }
        Ok(())
    }

    /// Shared driver of `call_TEI_candidate_sites_from_multiple_alignmts(_mosaic)`.
    #[allow(clippy::too_many_arguments)]
    fn call_sites_impl(
        &self,
        sf_annotation: &str,
        sf_rep_cns: &str,
        sf_ref: &str,
        b_se: bool,
        cutoff_left_clip: i64,
        cutoff_right_clip: i64,
        cutoff_clip_mate_in_rep: i64,
        b_mosaic: bool,
        mosaic_mode: bool,
        cutoff_polya: i64,
        sf_clip_folder: &str,
        b_force: bool,
        max_cov: i64,
        sf_out: &str,
    ) -> Result<()> {
        let l_lines = utils::read_lines(&self.sf_list)?;
        reject_x10_bams(&self.sf_list, &l_lines)?;
        let mut cnt = 0usize;
        let mut s_sample_bam = String::new();
        let mut b_set = false;
        let mut i_idx_bam = 0usize;
        for line in l_lines {
            let fields = split_ws(&line);
            if fields.len() < 2 {
                println!("Error input bam file {}!!!", utils::rstrip(&line));
                continue;
            }
            let sf_ori_bam = fields[0];
            let s_read_type = utils::rstrip(fields[1]);
            println!("Input bam {} is sequenced from {} platform!", sf_ori_bam, s_read_type);
            if mosaic_mode {
                if cutoff_left_clip <= 3 {
                    println!("Clip cutoff is small (<=3), keep all the clipped reads (initial cutoff set as 1)!!!");
                    gv::set_initial_min_clip_cutoff(1);
                } else {
                    gv::set_initial_min_clip_cutoff(gv::INITIAL_MIN_CLIP_CUTOFF_ILLUMINA);
                }
            } else if cutoff_left_clip <= 2 && b_mosaic {
                println!("Clip cutoff is small (<=2) , we are using 1 for initial cutoff");
                gv::set_initial_min_clip_cutoff(1);
            } else {
                gv::set_initial_min_clip_cutoff(gv::INITIAL_MIN_CLIP_CUTOFF_ILLUMINA);
            }
            if sf_ori_bam.len() <= 1 {
                continue;
            }
            if !b_set {
                s_sample_bam = sf_ori_bam.to_string();
                b_set = true;
            }
            let b_cutoff = true;
            let cutoff_hit_rep_copy = gv::initial_min_clip_cutoff();
            let sf_out_tmp = format!("{}{}{}", self.working_folder, gv::CLIP_TMP, cnt);
            cnt += 1;
            let caller = TeLocator::new(sf_ori_bam, &self.working_folder, self.n_jobs, &self.sf_ref);
            if sf_clip_folder.is_empty() {
                println!("public folder is null!!!!");
                continue;
            }
            let sf_new_pub = if sf_clip_folder.ends_with('/') { format!("{}{}/", sf_clip_folder, i_idx_bam) } else { format!("{}/{}/", sf_clip_folder, i_idx_bam) };
            if mosaic_mode {
                caller.call_tei_candidate_sites_from_clip_reads_v2_mosaic(
                    sf_annotation, sf_rep_cns, sf_ref, b_se, cutoff_hit_rep_copy, cutoff_hit_rep_copy, b_cutoff, &sf_new_pub, i_idx_bam, b_force, max_cov, &sf_out_tmp,
                )?;
            } else {
                caller.call_tei_candidate_sites_from_clip_reads_v2(
                    sf_annotation, sf_rep_cns, sf_ref, b_se, cutoff_hit_rep_copy, cutoff_hit_rep_copy, b_cutoff, &sf_new_pub, i_idx_bam, b_force, max_cov, &sf_out_tmp,
                )?;
            }
            i_idx_bam += 1;
        }
        let bam_info = BamInfo::new(&s_sample_bam, &self.sf_ref);
        let b_with_chr = bam_info.is_chrm_contain_chr()?;
        let m_chrms = bam_info.get_all_reference_names()?;
        let sf_out_merged = format!("{}_tmp", sf_out);
        let mut fout_merged = utils::create(&sf_out_merged)?;
        let mut fout = utils::create(sf_out)?;
        // read the per-bam outputs once
        let mut tmp_lines: Vec<Vec<String>> = Vec::new();
        for i in 0..cnt {
            let sf_tmp = format!("{}{}{}", self.working_folder, gv::CLIP_TMP, i);
            if !utils::is_file(&sf_tmp) {
                tmp_lines.push(Vec::new());
                continue;
            }
            tmp_lines.push(utils::read_lines(&sf_tmp)?);
        }
        for chrm in m_chrms.iter() {
            if is_decoy_contig_chrms(chrm) {
                continue;
            }
            let mut m_sites_chrm: IndexMap<i64, Vec<i64>> = IndexMap::new();
            for i in 0..cnt {
                let sf_tmp = format!("{}{}{}", self.working_folder, gv::CLIP_TMP, i);
                if !utils::is_file(&sf_tmp) {
                    println!("Errors happen, file {} doens't exist!", sf_tmp);
                    continue;
                }
                for line in tmp_lines[i].iter() {
                    let f = split_ws(line);
                    if f.is_empty() {
                        continue;
                    }
                    let tmp_chrm = crate::annotation::process_chrm_name(f[0], b_with_chr);
                    if &tmp_chrm != chrm {
                        continue;
                    }
                    let pos = utils::pint(f[1])?;
                    match m_sites_chrm.get_mut(&pos) {
                        None => {
                            let mut v = Vec::new();
                            for x in &f[2..] {
                                v.push(utils::pint(x)?);
                            }
                            m_sites_chrm.insert(pos, v);
                        }
                        Some(v) => {
                            for (k, x) in f[2..].iter().enumerate() {
                                v[k] += utils::pint(x)?;
                            }
                        }
                    }
                }
            }
            for (pos, v) in m_sites_chrm.iter() {
                write!(fout_merged, "{}\t{}\t", chrm, pos)?;
                for x in v {
                    write!(fout_merged, "{}\t", x)?;
                }
                writeln!(fout_merged)?;
            }
            let filtered = if mosaic_mode {
                xis::parse_sites_with_clip_cutoff_for_chrm_with_polya(&m_sites_chrm, cutoff_left_clip, cutoff_right_clip, cutoff_clip_mate_in_rep, cutoff_polya)
            } else {
                xis::parse_sites_with_clip_cutoff_for_chrm(&m_sites_chrm, cutoff_left_clip, cutoff_right_clip, cutoff_clip_mate_in_rep)
            };
            for (pos, v) in filtered.iter() {
                write!(fout, "{}\t{}\t", chrm, pos)?;
                for x in v {
                    write!(fout, "{}\t", x.py_str())?;
                }
                writeln!(fout)?;
            }
        }
        fout_merged.flush()?;
        fout.flush()?;
        Ok(())
    }

    /// `call_TEI_candidate_sites_from_multiple_alignmts` (step -C)
    #[allow(clippy::too_many_arguments)]
    pub fn call_tei_candidate_sites_from_multiple_alignmts(
        &self,
        sf_annotation: &str,
        sf_rep_cns: &str,
        sf_ref: &str,
        b_se: bool,
        cutoff_left_clip: i64,
        cutoff_right_clip: i64,
        cutoff_clip_mate_in_rep: i64,
        b_mosaic: bool,
        sf_clip_folder: &str,
        b_force: bool,
        max_cov: i64,
        sf_out: &str,
    ) -> Result<()> {
        self.call_sites_impl(
            sf_annotation, sf_rep_cns, sf_ref, b_se, cutoff_left_clip, cutoff_right_clip, cutoff_clip_mate_in_rep, b_mosaic, false, 0, sf_clip_folder, b_force, max_cov, sf_out,
        )
    }

    /// `call_TEI_candidate_sites_from_multiple_alignmts_mosaic` (step --mosaic)
    #[allow(clippy::too_many_arguments)]
    pub fn call_tei_candidate_sites_from_multiple_alignmts_mosaic(
        &self,
        sf_annotation: &str,
        sf_rep_cns: &str,
        sf_ref: &str,
        b_se: bool,
        cutoff_left_clip: i64,
        cutoff_right_clip: i64,
        cutoff_clip_mate_in_rep: i64,
        cutoff_polya: i64,
        sf_clip_folder: &str,
        b_force: bool,
        max_cov: i64,
        sf_out: &str,
    ) -> Result<()> {
        self.call_sites_impl(
            sf_annotation, sf_rep_cns, sf_ref, b_se, cutoff_left_clip, cutoff_right_clip, cutoff_clip_mate_in_rep, false, true, cutoff_polya, sf_clip_folder, b_force, max_cov, sf_out,
        )
    }

    /// `filter_candidate_sites_by_discordant_pairs_multi_alignmts` (step -D)
    #[allow(clippy::too_many_arguments)]
    pub fn filter_candidate_sites_by_discordant_pairs_multi_alignmts(
        &self,
        m_sites: &SiteMap<Vec<PyVal>>,
        iext: i64,
        i_is: i64,
        f_dev: f64,
        cutoff: i64,
        sf_annotation: &str,
        sf_out: &str,
        sf_raw_disc: &str,
        b_tumor: bool,
    ) -> Result<()> {
        let mut cnt = 0usize;
        for line in utils::read_lines(&self.sf_list)? {
            let fields = split_ws(&line);
            if fields.is_empty() {
                continue;
            }
            let sf_bam = fields[0];
            let caller = TeLocator::new(sf_bam, &self.working_folder, self.n_jobs, &self.sf_ref);
            let tmp_cutoff = 1;
            let (m_disc, m_raw) = caller.filter_candidate_sites_by_discordant_pairs_non_barcode(m_sites, iext, i_is, f_dev, sf_annotation, tmp_cutoff)?;
            let sf_out_tmp = format!("{}{}{}", self.working_folder, gv::DISC_TMP, cnt);
            xis::output_candidate_sites(&to_pyval_map(&m_disc), &sf_out_tmp)?;
            let sf_raw_tmp = format!("{}{}{}", self.working_folder, gv::RAW_DISC_TMP, cnt);
            xis::output_candidate_sites(&m_raw, &sf_raw_tmp)?;
            cnt += 1;
        }
        let mut merged: IndexMap<String, IndexMap<i64, [i64; 2]>> = IndexMap::new();
        for i in 0..cnt {
            let sf_tmp = format!("{}{}{}", self.working_folder, gv::DISC_TMP, i);
            for line in utils::read_lines(&sf_tmp)? {
                let f = split_ws(&line);
                let pos = utils::pint(f[1])?;
                let l = utils::pint(f[2])?;
                let r = utils::pint(f[3])?;
                let e = merged.entry(f[0].to_string()).or_default().entry(pos).or_insert([0, 0]);
                e[0] += l;
                e[1] += r;
            }
        }
        {
            let mut fout = utils::create(sf_out)?;
            let n_half_cutoff = cutoff as f64 / 2.0;
            for (chrm, m) in merged.iter() {
                for (pos, v) in m.iter() {
                    let (nl, nr) = (v[0], v[1]);
                    let keep = if b_tumor {
                        nl + nr >= cutoff
                    } else {
                        ((nl as f64) > n_half_cutoff && (nr as f64) > n_half_cutoff) || nl > cutoff || nr > cutoff
                    };
                    if keep {
                        writeln!(fout, "{}\t{}\t{}\t{}", chrm, pos, nl, nr)?;
                    }
                }
            }
            fout.flush()?;
        }
        if sf_raw_disc.is_empty() {
            return Ok(());
        }
        struct Raw {
            v: [i64; 4],
            lcl: Vec<String>,
            rcl: Vec<String>,
        }
        let mut merged_raw: IndexMap<String, IndexMap<i64, Raw>> = IndexMap::new();
        for i in 0..cnt {
            let sf_tmp = format!("{}{}{}", self.working_folder, gv::RAW_DISC_TMP, i);
            for line in utils::read_lines(&sf_tmp)? {
                let f = split_ws(&line);
                let pos = utils::pint(f[1])?;
                let vals = [utils::pint(f[2])?, utils::pint(f[3])?, utils::pint(f[4])?, utils::pint(f[5])?];
                let m = merged_raw.entry(f[0].to_string()).or_default();
                match m.get_mut(&pos) {
                    None => {
                        m.insert(pos, Raw { v: vals, lcl: vec![f[6].to_string()], rcl: vec![f[9].to_string()] });
                    }
                    Some(r) => {
                        for k in 0..4 {
                            r.v[k] += vals[k];
                        }
                        r.lcl.push(f[6].to_string());
                        r.rcl.push(f[9].to_string());
                    }
                }
            }
        }
        let mut i_half_cutoff = cutoff as f64 / 2.0;
        if !b_tumor {
            i_half_cutoff = (crate::py::fmt::py_int(cutoff as f64 * 3.0 / 4.0) + 1) as f64;
        }
        if i_half_cutoff < 1.0 {
            i_half_cutoff = 1.0;
        }
        let mut fout = utils::create(sf_raw_disc)?;
        for (chrm, m) in merged_raw.iter() {
            for (pos, r) in m.iter() {
                let bl = (r.v[0] as f64) >= i_half_cutoff && r.lcl.iter().any(|s| s == "1");
                let br = (r.v[1] as f64) >= i_half_cutoff && r.rcl.iter().any(|s| s == "1");
                if bl || br {
                    writeln!(fout, "{}\t{}\t{}\t{}\t{}\t{}", chrm, pos, r.v[0], r.v[1], r.v[2], r.v[3])?;
                }
            }
        }
        fout.flush()?;
        Ok(())
    }
}

fn to_pyval_map(m: &SiteMap<Vec<i64>>) -> SiteMap<Vec<PyVal>> {
    let mut out: SiteMap<Vec<PyVal>> = IndexMap::new();
    for (c, mp) in m.iter() {
        let e = out.entry(c.clone()).or_default();
        for (p, v) in mp.iter() {
            e.insert(*p, v.iter().map(|x| PyVal::Int(*x)).collect());
        }
    }
    out
}

pub struct TeLocator {
    pub sf_bam: String,
    pub working_folder: String,
    pub n_jobs: usize,
    pub sf_reference: String,
}

impl TeLocator {
    pub fn new(sf_bam: &str, working_folder: &str, n_jobs: usize, sf_ref: &str) -> Self {
        TeLocator {
            sf_bam: sf_bam.to_string(),
            working_folder: working_folder.to_string(),
            n_jobs,
            sf_reference: sf_ref.to_string(),
        }
    }

    fn pool(&self) -> rayon::ThreadPool {
        rayon::ThreadPoolBuilder::new().num_threads(self.n_jobs.max(1)).build().unwrap()
    }

    /// `call_TEI_candidate_sites_from_clip_reads_v2`
    #[allow(clippy::too_many_arguments)]
    pub fn call_tei_candidate_sites_from_clip_reads_v2(
        &self,
        sf_annotation: &str,
        sf_rep_cns: &str,
        sf_ref: &str,
        b_se: bool,
        cutoff_left_clip: i64,
        cutoff_right_clip: i64,
        b_cutoff: bool,
        sf_pub_folder0: &str,
        idx_bam: usize,
        b_force: bool,
        max_cov_cutoff: i64,
        sf_out: &str,
    ) -> Result<()> {
        let sf_pub_folder = utils::with_slash(sf_pub_folder0);
        if !utils::exists(&sf_pub_folder) {
            utils::run_cmd_small_output(&format!("mkdir -p {}", sf_pub_folder));
        }
        let sf_clip_working_folder = format!("{}{}/{}/", self.working_folder, gv::CLIP_FOLDER, idx_bam);
        if !utils::exists(&sf_clip_working_folder) {
            utils::run_cmd_small_output(&format!("mkdir -p {}", sf_clip_working_folder));
        }
        let mut clip_info = ClipReadInfo::new(&self.sf_bam, self.n_jobs, &self.sf_reference);
        let sf_bam_name = utils::basename(&self.sf_bam);
        let sf_all_clip_fq = format!("{}{}{}", sf_pub_folder, sf_bam_name, gv::CLIP_FQ_SUFFIX);
        clip_info.set_working_folder(&sf_clip_working_folder);
        if !utils::is_link(&sf_all_clip_fq) || b_force {
            println!("Collected clipped reads file {} doesn't exist. Generate it now!", sf_all_clip_fq);
            let initial = gv::initial_min_clip_cutoff();
            println!("Initial minimum clip cutoff is {}", initial);
            clip_info.collect_clip_positions(sf_annotation, initial, b_se, &sf_pub_folder)?;
            println!("({}, {})", crate::py::fmt::py_str_repr("Output info: Collect clipped parts for file "), crate::py::fmt::py_str_repr(&self.sf_bam));
            let sf_all_clip_fq_ori = format!("{}{}{}", sf_clip_working_folder, sf_bam_name, gv::CLIP_FQ_SUFFIX);
            clip_info.collect_clipped_parts(&sf_all_clip_fq_ori)?;
            if utils::is_file(&sf_all_clip_fq) || utils::is_link(&sf_all_clip_fq) {
                utils::remove_file(&sf_all_clip_fq);
            }
            utils::run_cmd_small_output(&format!("ln -s {} {}", sf_all_clip_fq_ori, sf_all_clip_fq));
        } else {
            println!("Collected clipped reads file {} already exist!", sf_all_clip_fq);
        }
        let sf_algnmt = format!("{}{}{}", self.working_folder, sf_bam_name, gv::CLIP_BAM_SUFFIX);
        println!("({}, {})", crate::py::fmt::py_str_repr("Output info: Re-align clipped parts for file "), crate::py::fmt::py_str_repr(&self.sf_bam));
        let bwa = BwAlign::new(gv::BWA_PATH, gv::BWA_REALIGN_CUTOFF, self.n_jobs as i64);
        bwa.two_stage_realign(sf_rep_cns, sf_ref, &sf_all_clip_fq, &sf_algnmt)?;
        clip_info.cnt_clip_part_aligned_to_rep(&sf_algnmt)?;
        if !b_cutoff {
            clip_info.merge_clip_positions(&sf_pub_folder, sf_out)?;
        } else {
            clip_info.merge_clip_positions_with_cutoff(cutoff_left_clip, cutoff_right_clip, max_cov_cutoff, &sf_pub_folder, sf_out)?;
        }
        if utils::is_file(&sf_algnmt) {
            utils::remove_file(&sf_algnmt);
        }
        Ok(())
    }

    /// `call_TEI_candidate_sites_from_clip_reads_v2_mosaic`
    #[allow(clippy::too_many_arguments)]
    pub fn call_tei_candidate_sites_from_clip_reads_v2_mosaic(
        &self,
        sf_annotation: &str,
        sf_rep_cns: &str,
        sf_ref: &str,
        b_se: bool,
        cutoff_left_clip: i64,
        cutoff_right_clip: i64,
        b_cutoff: bool,
        sf_pub_folder0: &str,
        idx_bam: usize,
        _b_force: bool,
        max_cov_cutoff: i64,
        sf_out: &str,
    ) -> Result<()> {
        let sf_pub_folder = utils::with_slash(sf_pub_folder0);
        if !utils::exists(&sf_pub_folder) {
            utils::run_cmd_small_output(&format!("mkdir -p {}", sf_pub_folder));
        }
        let sf_clip_working_folder = format!("{}{}/{}/", self.working_folder, gv::CLIP_FOLDER, idx_bam);
        if !utils::exists(&sf_clip_working_folder) {
            utils::run_cmd_small_output(&format!("mkdir -p {}", sf_clip_working_folder));
        }
        let mut clip_info = ClipReadInfo::new(&self.sf_bam, self.n_jobs, &self.sf_reference);
        let sf_bam_name = utils::basename(&self.sf_bam);
        let sf_all_clip_fq = format!("{}{}{}", sf_pub_folder, sf_bam_name, gv::CLIP_FQ_SUFFIX);
        clip_info.set_working_folder(&sf_clip_working_folder);
        let initial = gv::initial_min_clip_cutoff();
        println!("Initial minimum clip cutoff is {}", initial);
        clip_info.collect_clip_positions(sf_annotation, initial, b_se, &sf_pub_folder)?;
        println!("({}, {})", crate::py::fmt::py_str_repr("Output info: Collect clipped parts for file "), crate::py::fmt::py_str_repr(&self.sf_bam));
        let sf_all_clip_fq_ori = format!("{}{}{}", sf_clip_working_folder, sf_bam_name, gv::CLIP_FQ_SUFFIX);
        clip_info.collect_clipped_parts(&sf_all_clip_fq_ori)?;
        if utils::is_file(&sf_all_clip_fq) {
            utils::remove_file(&sf_all_clip_fq);
        }
        utils::run_cmd_small_output(&format!("ln -s {} {}", sf_all_clip_fq_ori, sf_all_clip_fq));
        let sf_algnmt = format!("{}{}{}", self.working_folder, sf_bam_name, gv::CLIP_BAM_SUFFIX);
        println!("({}, {})", crate::py::fmt::py_str_repr("Output info: Re-align clipped parts for file "), crate::py::fmt::py_str_repr(&self.sf_bam));
        let bwa = BwAlign::new(gv::BWA_PATH, gv::BWA_REALIGN_CUTOFF, self.n_jobs as i64);
        bwa.two_stage_realign(sf_rep_cns, sf_ref, &sf_all_clip_fq, &sf_algnmt)?;
        clip_info.cnt_clip_part_aligned_to_rep_mosaic(&sf_algnmt)?;
        if !b_cutoff {
            clip_info.merge_clip_positions(&sf_pub_folder, sf_out)?;
        } else {
            clip_info.merge_clip_positions_with_cutoff_polya(cutoff_left_clip, cutoff_right_clip, max_cov_cutoff, &sf_pub_folder, sf_out)?;
        }
        if utils::is_file(&sf_algnmt) {
            utils::remove_file(&sf_algnmt);
        }
        Ok(())
    }

    /// `collect_all_clipped_reads_only`
    pub fn collect_all_clipped_reads_only(&self, sf_annotation: &str, b_se: bool, s_working_folder: &str, wfolder_pub_clip: &str) -> Result<()> {
        let mut wf = format!("{}{}/", s_working_folder, gv::CLIP_FOLDER);
        if wf.len() > 1 && !wf.ends_with('/') {
            wf.push('/');
        }
        if !utils::exists(&wf) {
            utils::run_cmd_small_output(&format!("mkdir {}", wf));
        }
        let mut clip_info = ClipReadInfo::new(&self.sf_bam, self.n_jobs, &self.sf_reference);
        clip_info.set_working_folder(&wf);
        println!("({}, {})", crate::py::fmt::py_str_repr("Output info: Collect clip positions for file "), crate::py::fmt::py_str_repr(&self.sf_bam));
        clip_info.collect_clip_positions(sf_annotation, 2, b_se, wfolder_pub_clip)?;
        let sf_bam_name = utils::basename(&self.sf_bam);
        let sf_all_clip_fq = format!("{}{}{}", wf, sf_bam_name, gv::CLIP_FQ_SUFFIX);
        println!("({}, {})", crate::py::fmt::py_str_repr("Output info: Collect clipped parts for file "), crate::py::fmt::py_str_repr(&self.sf_bam));
        clip_info.collect_clipped_parts(&sf_all_clip_fq)?;
        Ok(())
    }

    /// `output_candidate_sites_by_chrm`
    pub fn output_candidate_sites_by_chrm(&self, m: &SiteMap<Vec<PyVal>>, sf_folder: &str, s_suffix: &str) -> Result<()> {
        for (chrm, mp) in m.iter() {
            let mut f = utils::create(&format!("{}{}{}", sf_folder, chrm, s_suffix))?;
            for (pos, v) in mp.iter() {
                write!(f, "{}\t{}\t", chrm, pos)?;
                for x in v {
                    write!(f, "{}\t", x.py_str())?;
                }
                writeln!(f)?;
            }
            f.flush()?;
        }
        Ok(())
    }

    /// `run_filter_by_discordant_pair_by_chrom_non_barcode`
    #[allow(clippy::too_many_arguments)]
    fn run_filter_by_discordant_pair_by_chrom_non_barcode(
        &self,
        site_chrm1: &str,
        iextend: i64,
        i_is: i64,
        f_dev: i64,
        xannotation: &XAnnotation,
        sf_disc_working_folder: &str,
        s_suffix: &str,
    ) -> Result<()> {
        let sf_candidate_list = format!("{}{}{}", sf_disc_working_folder, site_chrm1, s_suffix);
        if !utils::exists(&sf_candidate_list) {
            return Ok(());
        }
        let mut m_candidate_pos: IndexMap<i64, String> = IndexMap::new();
        for line in utils::read_lines(&sf_candidate_list)? {
            let f = split_ws(&line);
            let pos = utils::pint(f[1])?;
            m_candidate_pos.insert(pos, f[2..].join("\t"));
        }
        let bam_info = BamInfo::new(&self.sf_bam, &self.sf_reference);
        let b_with_chr = bam_info.is_chrm_contain_chr()?;
        let m_chrms = bam_info.get_all_reference_names()?;
        let site_chrm = crate::annotation::process_chrm_name(site_chrm1, b_with_chr);
        if !m_chrms.contains(&site_chrm) {
            return Ok(());
        }
        let mut bamfile = AlignmentFile::open(&self.sf_bam, &self.sf_reference)?;
        let m_chrm_ids = bamfile.chrm_id_name();
        let mut m_new: IndexMap<i64, [i64; 2]> = IndexMap::new();
        let mut m_raw: IndexMap<i64, Vec<String>> = IndexMap::new();
        for (&site_pos, _) in m_candidate_pos.iter() {
            if site_pos < iextend {
                continue;
            }
            let (nl, nlraw, lcl) = cnt_discordant_pairs(&mut bamfile, &m_chrm_ids, &site_chrm, site_pos - iextend, site_pos, i_is, f_dev, xannotation)?;
            let (nr, nrraw, rcl) = cnt_discordant_pairs(&mut bamfile, &m_chrm_ids, &site_chrm, site_pos + 1, site_pos + iextend, i_is, f_dev, xannotation)?;
            m_new.insert(site_pos, [nl, nr]);
            let (mut s_lc, mut s_lchrm, mut s_lpos) = ("0".to_string(), "-1".to_string(), "-1".to_string());
            if lcl.0 {
                s_lc = "1".to_string();
                s_lchrm = lcl.1.clone().unwrap_or_default();
                s_lpos = lcl.2.unwrap_or(-1).to_string();
            }
            let (mut s_rc, mut s_rchrm, mut s_rpos) = ("0".to_string(), "-1".to_string(), "-1".to_string());
            if rcl.0 {
                s_rc = "1".to_string();
                s_rchrm = rcl.1.clone().unwrap_or_default();
                s_rpos = rcl.2.unwrap_or(-1).to_string();
            }
            m_raw.insert(
                site_pos,
                vec![nlraw.to_string(), nrraw.to_string(), nl.to_string(), nr.to_string(), s_lc, s_lchrm, s_lpos, s_rc, s_rchrm, s_rpos],
            );
        }
        {
            let mut f = utils::create(&format!("{}{}", sf_candidate_list, gv::DISC_SUFFIX_FILTER))?;
            for (pos, v) in m_new.iter() {
                write!(f, "{}\t{}\t{}\t\n", pos, v[0], v[1])?;
            }
            f.flush()?;
        }
        {
            let mut f = utils::create(&format!("{}{}", sf_candidate_list, gv::RAW_DISC_SUFFIX_FILTER))?;
            for (pos, v) in m_raw.iter() {
                write!(f, "{}\t", pos)?;
                for x in v {
                    write!(f, "{}\t", x)?;
                }
                writeln!(f)?;
            }
            f.flush()?;
        }
        Ok(())
    }

    /// `filter_candidate_sites_by_discordant_pairs_non_barcode`
    #[allow(clippy::type_complexity)]
    pub fn filter_candidate_sites_by_discordant_pairs_non_barcode(
        &self,
        m_candidate_sites: &SiteMap<Vec<PyVal>>,
        iextend: i64,
        i_is: i64,
        f_dev: f64,
        sf_annotation: &str,
        n_discordant_cutoff: i64,
    ) -> Result<(SiteMap<Vec<i64>>, SiteMap<Vec<PyVal>>)> {
        let mut folder = format!("{}{}", self.working_folder, gv::DISC_FOLDER);
        if !utils::exists(&folder) {
            utils::run_cmd_small_output(&format!("mkdir {}", folder));
        }
        folder.push('/');
        self.output_candidate_sites_by_chrm(m_candidate_sites, &folder, gv::DISC_SUFFIX)?;
        let chrms: Vec<String> = m_candidate_sites.keys().filter(|c| !is_decoy_contig_chrms(c)).cloned().collect();
        // annotation (identical for all chromosomes; loaded once instead of per worker)
        let bam_info = BamInfo::new(&self.sf_bam, &self.sf_reference);
        let b_with_chr = bam_info.is_chrm_contain_chr()?;
        let mut xann = XAnnotation::new(sf_annotation);
        xann.set_with_chr(b_with_chr);
        let boundary_extnd = if gv::is_call_sva() { gv::SVA_ANNOTATION_EXTND } else { 0 };
        xann.load_rmsk_annotation_with_extnd_with_lenth_cutoff(boundary_extnd, 0)?;
        xann.index_rmsk_annotation_interval_tree();
        let xann = Arc::new(xann);
        // record: int(record[2]) etc.
        let iext = iextend;
        let iis = i_is;
        let idev = crate::py::fmt::py_int(f_dev);
        let results: Vec<Result<()>> = self.pool().install(|| {
            utils::par_map(chrms.len(), |i| self.run_filter_by_discordant_pair_by_chrom_non_barcode(&chrms[i], iext, iis, idev, &xann, &folder, gv::DISC_SUFFIX))
        });
        for r in results {
            r?;
        }
        let mut m_new: SiteMap<Vec<i64>> = IndexMap::new();
        for chrm in m_candidate_sites.keys() {
            let sf = format!("{}{}{}{}", folder, chrm, gv::DISC_SUFFIX, gv::DISC_SUFFIX_FILTER);
            if !utils::exists(&sf) {
                continue;
            }
            for line in utils::read_lines(&sf)? {
                let f = split_ws(&line);
                let pos = utils::pint(f[0])?;
                let l = utils::pint(f[1])?;
                let r = utils::pint(f[2])?;
                if l + r < n_discordant_cutoff {
                    continue;
                }
                let e = m_new.entry(chrm.clone()).or_default();
                e.entry(pos).or_insert_with(|| vec![l, r]);
            }
        }
        let mut m_raw: SiteMap<Vec<PyVal>> = IndexMap::new();
        for chrm in m_candidate_sites.keys() {
            let sf = format!("{}{}{}{}", folder, chrm, gv::DISC_SUFFIX, gv::RAW_DISC_SUFFIX_FILTER);
            if !utils::exists(&sf) {
                continue;
            }
            for line in utils::read_lines(&sf)? {
                let f = split_ws(&line);
                let pos = utils::pint(f[0])?;
                let mut v: Vec<PyVal> = Vec::new();
                for k in 1..5 {
                    v.push(PyVal::Int(utils::pint(f[k])?));
                }
                for k in 5..11 {
                    v.push(PyVal::Str(f[k].to_string()));
                }
                let e = m_raw.entry(chrm.clone()).or_default();
                e.entry(pos).or_insert(v);
            }
        }
        Ok((m_new, m_raw))
    }

}
