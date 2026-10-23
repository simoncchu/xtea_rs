//! Case-control (somatic) calling: port of `x_somatic_calling.py` (`CaseControlMode`, and the
//! parts of `SomaticMEICaller` used by the mosaic caller).

use crate::clip_disc_filter::{ClipChk, XClipDiscFilter};
use crate::coverage::ReadDepth;
use crate::genotype_feature::{GntpFeat, XGenotyper};
use crate::globals as gv;
use crate::intermediate_sites::SiteMap;
use crate::num::PyNum;
use crate::orphan_transduction::XOrphanTransduction;
use crate::py::text::{py_float, py_int};
use crate::transduction::XTransduction;
use crate::utils::{self, split_ws};
use anyhow::{anyhow, Result};
use indexmap::IndexMap;
use std::collections::HashSet;
use std::io::Write;

fn fld<'a>(f: &[&'a str], i: usize) -> Result<&'a str> {
    f.get(i).copied().ok_or_else(|| anyhow!("IndexError: list index out of range"))
}

/// `int(float(s))`
fn int_float(s: &str) -> Result<i64> {
    let v = py_float(s)?;
    if !v.is_finite() {
        return Err(anyhow!("cannot convert float {} to integer", s));
    }
    Ok(v.trunc() as i64)
}

/// Python float division (ZeroDivisionError -> Err).
fn fdiv(a: f64, b: f64) -> Result<f64> {
    if b == 0.0 {
        return Err(anyhow!("ZeroDivisionError: float division by zero"));
    }
    Ok(a / b)
}

fn site_has<T>(m: &SiteMap<T>, chrm: &str, pos: i64) -> bool {
    m.get(chrm).map_or(false, |mm| mm.contains_key(&pos))
}

/// The record of `load_TEI_info_from_file`.
type TeiInfo = (i64, i64, i64, i64, i64, i64, String, String, String, String);

/// `CaseControlMode`
pub struct CaseControlMode {
    pub s_wfolder: String,
    pub sf_ref: String,
    pub n_jobs: usize,
    pub iextnd: i64,
    pub bin_size: i64,
    pub bmapped_cutoff: f64,
    pub i_concord_dist: i64,
    pub f_concord_ratio: f64,
    clip_consist_dist: i64,
    disc_consist_dist: i64,
}

impl CaseControlMode {
    pub fn new(sf_ref: &str, s_wfolder: &str, n_jobs: usize) -> Self {
        // (sic) a missing trailing "/" is appended to the local variable only
        let wf = if utils::exists(s_wfolder) { s_wfolder.to_string() } else { "./".to_string() };
        CaseControlMode {
            s_wfolder: wf,
            sf_ref: sf_ref.to_string(),
            n_jobs,
            iextnd: 400,
            bin_size: 50000000,
            bmapped_cutoff: 0.65,
            i_concord_dist: 550,
            f_concord_ratio: 0.45,
            clip_consist_dist: 35,
            disc_consist_dist: 50,
        }
    }

    /// `load_sites`: chrm -> pos -> the remaining fields (first occurrence kept).
    pub fn load_sites(sf_ce: &str) -> Result<SiteMap<String>> {
        let mut m: SiteMap<String> = IndexMap::new();
        for line in utils::read_raw_lines(sf_ce)? {
            let fields = split_ws(&line);
            let chrm = fld(&fields, 0)?;
            let pos = py_int(fld(&fields, 1)?)?;
            let sinfo = fields[2..].join("\t");
            m.entry(chrm.to_string()).or_default().entry(pos).or_insert(sinfo);
        }
        Ok(m)
    }

    /// `extract_sites_raw_gntp_feature_depth`: depth and genotype features of the sites in the
    /// control BAMs.
    pub fn extract_sites_raw_gntp_feature_depth(&self, sf_candidate_list: &str, sf_bam_list: &str) -> Result<(SiteMap<[PyNum; 4]>, SiteMap<GntpFeat>)> {
        let rd = ReadDepth::new(&self.s_wfolder, self.n_jobs, &self.sf_ref);
        let m_sites = Self::load_sites(sf_candidate_list)?;
        let l_sites: Vec<(String, i64)> =
            m_sites.iter().flat_map(|(c, m)| m.keys().map(move |&p| (c.clone(), p))).collect();
        let m_read_depth =
            rd.calc_coverage_of_two_regions(&l_sites, sf_bam_list, gv::COV_SEARCH_WINDOW, gv::LOCAL_COV_WIN, gv::COV_ISD_CHK_WIN)?;
        ReadDepth::dump_coverage_info(&m_read_depth, &format!("{}.read_depth", sf_candidate_list))?;
        let x_gntper = XGenotyper::new(&self.sf_ref, &self.s_wfolder, self.n_jobs);
        let sf_gntp_feature = format!("{}.gntp.features", sf_candidate_list);
        x_gntper.call_genotype_f(sf_bam_list, sf_candidate_list, gv::dft_is(), &sf_gntp_feature)?;
        let i_max_cov = gv::MAX_COV_TIMES as f64 * gv::average_cov();
        let m_gntp_info = XGenotyper::load_in_features_from_file_with_cov_cutoff(&sf_gntp_feature, i_max_cov)?;
        Ok((m_read_depth, m_gntp_info))
    }

    /// `call_somatic_TE_insertion`: filter the case calls with the control BAMs.
    #[allow(clippy::too_many_arguments)]
    pub fn call_somatic_te_insertion(
        &self,
        sf_bam_list: &str,
        sf_case_candidates: &str,
        extnd: i64,
        nclip_cutoff: i64,
        ndisc_cutoff: i64,
        npolya_cutoff: i64,
        sf_rep_cns: &str,
        sf_flank: &str,
        i_flk_len: i64,
        bin_size: i64,
        sf_out: &str,
    ) -> Result<()> {
        let xclip_disc = XClipDiscFilter::new(sf_bam_list, &self.s_wfolder, self.n_jobs, &self.sf_ref);
        let sf_non_td = format!("{}.non_td.tmp", sf_case_candidates);
        let sf_td = format!("{}.td.tmp", sf_case_candidates);
        let sf_td_sibling = format!("{}.td_sibling.tmp", sf_case_candidates);
        let sf_orphan = format!("{}.orphan.tmp", sf_case_candidates);
        XClipDiscFilter::sprt_tei_to_td_orphan_non_td(sf_case_candidates, &sf_non_td, &sf_td, &sf_td_sibling, &sf_orphan)?;

        let sf_tmp_cluster = format!("{}.candidate_somatic_cluster_from_ctrl.txt", sf_out);
        let m_clip_cluster_non_td = xclip_disc.call_clip_disc_cluster(
            &sf_non_td,
            self.iextnd,
            self.bin_size,
            sf_rep_cns,
            self.bmapped_cutoff,
            self.i_concord_dist,
            self.f_concord_ratio,
            nclip_cutoff as f64,
            ndisc_cutoff as f64,
            &self.s_wfolder,
            &sf_tmp_cluster,
        )?;
        let sf_tmp_td_cluster = format!("{}.candidate_somatic_td_cluster_from_ctrl.txt", sf_out);
        let m_clip_cluster_td = xclip_disc.call_clip_disc_cluster(
            &sf_td,
            self.iextnd,
            self.bin_size,
            sf_rep_cns,
            self.bmapped_cutoff,
            self.i_concord_dist,
            self.f_concord_ratio,
            nclip_cutoff as f64,
            ndisc_cutoff as f64,
            &self.s_wfolder,
            &sf_tmp_td_cluster,
        )?;

        let xtd = XTransduction::new(&self.s_wfolder, self.n_jobs, &self.sf_ref);
        let sf_tmp_cluster_td = format!("{}.candidate_somatic_cluster_from_ctrl_td.txt", sf_out);
        let (m_td_cluster, m_td_polya) = xtd.collect_realign_reads(
            &sf_td,
            sf_case_candidates,
            &xclip_disc,
            sf_flank,
            i_flk_len,
            extnd,
            bin_size,
            ndisc_cutoff,
            sf_rep_cns,
            &sf_tmp_cluster_td,
        )?;
        {
            let mut fout = utils::create(&format!("{}.tmp_ctrl_clip_polyA", sf_out))?;
            for (ins_chrm, m) in m_td_polya.iter() {
                for (ins_pos, rcd) in m.iter() {
                    writeln!(fout, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}", ins_chrm, ins_pos, rcd[0], rcd[1], rcd[2], rcd[3], rcd[4], rcd[5])?;
                }
            }
            fout.flush()?;
        }

        let xorphan = XOrphanTransduction::new(&self.s_wfolder, self.n_jobs, &self.sf_ref);
        let m_need_filter_out = xorphan.check_features_for_given_sites(&sf_orphan, sf_bam_list, nclip_cutoff, ndisc_cutoff)?;
        let m_filter_out_sibling = xorphan.check_features_for_given_sites(&sf_td_sibling, sf_bam_list, nclip_cutoff, ndisc_cutoff)?;

        let m_cluster_info = XClipDiscFilter::load_tei_info_from_file(&sf_tmp_cluster)?;
        let m_cluster_info_td = XClipDiscFilter::load_tei_info_from_file(&sf_tmp_td_cluster)?;
        let sf_before_final = format!("{}.before_final_filter", sf_out);
        {
            let mut m_unique_tmp: HashSet<(String, i64)> = HashSet::new();
            let mut fout = utils::create(&sf_before_final)?;
            for line in utils::read_raw_lines(&sf_non_td)? {
                let fields = split_ws(&line);
                let chrm = fld(&fields, 0)?;
                let pos = py_int(fld(&fields, 1)?)?;
                if !m_unique_tmp.insert((chrm.to_string(), pos)) {
                    continue;
                }
                if self.pass_cns_chk(&m_cluster_info, &m_clip_cluster_non_td, &fields, npolya_cutoff)? {
                    fout.write_all(line.as_bytes())?;
                }
            }
            for line in utils::read_raw_lines(&sf_td)? {
                let l_fields = split_ws(&line);
                let chrm = fld(&l_fields, 0)?;
                let pos = py_int(fld(&l_fields, 1)?)?;
                if !m_unique_tmp.insert((chrm.to_string(), pos)) {
                    continue;
                }
                if self.pass_td_chk(&m_td_cluster, &m_td_polya, &l_fields, npolya_cutoff)?
                    && self.pass_cns_chk(&m_cluster_info_td, &m_clip_cluster_td, &l_fields, npolya_cutoff)?
                {
                    if site_has(&m_filter_out_sibling, chrm, pos) {
                        continue;
                    }
                    fout.write_all(line.as_bytes())?;
                }
            }
            for line in utils::read_raw_lines(&sf_orphan)? {
                let fields = split_ws(&line);
                let chrm = fld(&fields, 0)?;
                let pos = py_int(fld(&fields, 1)?)?;
                if !m_unique_tmp.insert((chrm.to_string(), pos)) {
                    continue;
                }
                if site_has(&m_need_filter_out, chrm, pos) {
                    continue;
                }
                fout.write_all(line.as_bytes())?;
            }
            fout.flush()?;
        }

        let f_disc_ef_ratio = 0.1;
        let f_ef_disc_cutoff = 0.05;
        let f_disc_cutoff = 0.15;
        let f_case_ctrl_cutoff = 0.35;
        let f_ctrl_clip_ratio = 0.75;
        let f_ctrl_raw_clip_flex_ratio = 0.35;
        let f_ctrl_polya_cutoff = 0.3;
        let f_case_raw_clip_cutoff = 0.05;
        self.filter_by_case_ctrl_disc_depth_ratio(
            &sf_before_final,
            sf_bam_list,
            f_disc_ef_ratio,
            f_ef_disc_cutoff,
            f_disc_cutoff,
            f_case_raw_clip_cutoff,
            f_case_ctrl_cutoff,
            f_ctrl_clip_ratio,
            f_ctrl_raw_clip_flex_ratio,
            f_ctrl_polya_cutoff,
            sf_out,
        )
    }

    /// `filter_by_case_ctrl_disc_depth_ratio`: drop germline events by comparing the case and
    /// control discordant / clip / polyA ratios.
    #[allow(clippy::too_many_arguments)]
    pub fn filter_by_case_ctrl_disc_depth_ratio(
        &self,
        sf_candidate: &str,
        sf_ctrl_bam: &str,
        f_disc_ef_ratio: f64,
        f_ef_disc_depth_cutoff: f64,
        f_disc_cutoff: f64,
        f_case_raw_clip_cutoff: f64,
        f_case_ctrl_cutoff: f64,
        f_ctrl_clip_ratio: f64,
        f_ctrl_raw_clip_flex_ratio: f64,
        f_ctrl_polya_cutoff: f64,
        sf_out: &str,
    ) -> Result<()> {
        let (m_depth_ctrl, m_feature_ctrl) = self.extract_sites_raw_gntp_feature_depth(sf_candidate, sf_ctrl_bam)?;
        let lines = utils::read_raw_lines(sf_candidate)?;
        let mut fout_new = utils::create(sf_out)?;
        let mut fout_log = utils::create(&format!("{}.filter_log", sf_out))?;
        for line in lines.iter() {
            let b_trsdct = !line.contains(gv::NOT_TRANSDUCTION);
            let b_orphan = line.contains(gv::ORPHAN_LABEL);
            let b_sibling = line.contains(gv::SIBLING_LABEL);
            let fields = split_ws(line);
            let gi = |i: usize| -> Result<i64> { py_int(fld(&fields, i)?) };
            let gf = |i: usize| -> Result<f64> { py_float(fld(&fields, i)?) };
            let ins_chrm = fld(&fields, 0)?;
            let ins_pos = gi(1)?;
            let n_l_ef_disc = gi(7)?;
            let n_r_ef_disc = gi(8)?;
            let n_lpolya = gi(9)?;
            let n_rpolya = gi(10)?;
            let f_lcov = gf(11)?;
            let f_rcov = gf(12)?;
            let n_ef_clip = gi(5)? + gi(6)?;
            let n_ef_disc = gi(7)? + gi(8)?;
            let _n_clip = gi(35)?;
            let _n_full_map = gi(36)?;
            let n_raw_lclip_case = gi(37)?;
            let n_raw_rclip_case = gi(38)?;
            let n_disc = gi(39)?;
            let _n_concod = gi(40)?;

            if f_lcov + f_rcov <= 0.0000000001 {
                continue;
            }
            let f_ef_disc_ratio_case = n_ef_disc as f64 * 2.0 / (f_lcov + f_rcov);
            if f_ef_disc_ratio_case < f_ef_disc_depth_cutoff {
                writeln!(fout_log, "{}:{} is filtered out as ef_disc_ratio_case is small!", ins_chrm, ins_pos)?;
                continue;
            }
            if b_orphan
                && ((f_lcov > 0.0 && gf(7)? / f_lcov < f_ef_disc_depth_cutoff)
                    || (f_rcov > 0.0 && gf(8)? / f_rcov < f_ef_disc_depth_cutoff))
            {
                writeln!(fout_log, "{}:{} (orpha) is filtered out as ef_disc_ratio_case is small!", ins_chrm, ins_pos)?;
                continue;
            } else if b_orphan
                && ((f_lcov > 0.0 && n_raw_lclip_case as f64 / f_lcov < f_case_raw_clip_cutoff)
                    || (f_rcov > 0.0 && n_raw_rclip_case as f64 / f_rcov < f_case_raw_clip_cutoff))
            {
                writeln!(fout_log, "{}:{} (orpha) is filtered out as raw clip ratio is small!", ins_chrm, ins_pos)?;
                continue;
            }

            let f_disc_ratio_case = n_disc as f64 * 2.0 / (f_lcov + f_rcov);
            if f_disc_ratio_case < f_disc_cutoff {
                writeln!(fout_log, "{}:{} is filtered out as disc_ratio_case is small!", ins_chrm, ins_pos)?;
                continue;
            }
            if n_disc <= 0 {
                writeln!(fout_log, "{}:{} is filtered out as no discordant reads!", ins_chrm, ins_pos)?;
                continue;
            }
            let f_ef_disc_raw_disc_ratio = (n_l_ef_disc + n_r_ef_disc) as f64 / n_disc as f64;
            if f_ef_disc_raw_disc_ratio < f_disc_ef_ratio && !b_orphan && !b_sibling {
                writeln!(fout_log, "{}:{} is filtered out as ef_disc_raw_disc_ratio is small!", ins_chrm, ins_pos)?;
                continue;
            }
            let f_clip_ratio_case = (n_raw_lclip_case + n_raw_rclip_case) as f64 * 2.0 / (f_lcov + f_rcov);
            let f_ef_clip_ratio_case = n_ef_clip as f64 * 2.0 / (f_lcov + f_rcov);

            let feat = match m_feature_ctrl.get(ins_chrm).and_then(|m| m.get(&ins_pos)) {
                Some(f) => f,
                None => {
                    writeln!(fout_log, "{}:{} is filtered out as site is not in m_feature_ctrl!", ins_chrm, ins_pos)?;
                    continue;
                }
            };
            let n_af_clip_ctrl = feat.n_af_clip;
            let n_full_map_ctrl = feat.n_full_map;
            let n_raw_lclip_ctrl = feat.n_l_raw_clip;
            let n_raw_rclip_ctrl = feat.n_r_raw_clip;
            let n_disc_pairs_ctrl = feat.n_disc_pairs;
            let n_polya_ctrl = feat.n_polya;
            let depth = match m_depth_ctrl.get(ins_chrm).and_then(|m| m.get(&ins_pos)) {
                Some(d) => d,
                None => {
                    writeln!(fout_log, "{}:{} is filtered out as site is not in m_depth_ctrl!", ins_chrm, ins_pos)?;
                    continue;
                }
            };
            let f_lcov_ctrl = depth[0].f();
            let f_rcov_ctrl = depth[1].f();
            if f_lcov_ctrl + f_rcov_ctrl <= 0.0000000001 {
                writeln!(fout_log, "{}:{} is filtered out as ctrl depth is 0!", ins_chrm, ins_pos)?;
                continue;
            }
            let f_disc_ratio_ctrl = n_disc_pairs_ctrl as f64 * 2.0 / (f_lcov_ctrl + f_rcov_ctrl);
            let f_clip_ratio_ctrl = (n_raw_lclip_ctrl + n_raw_rclip_ctrl) as f64 * 2.0 / (f_lcov_ctrl + f_rcov_ctrl);
            let f_ef_clip_ratio_ctrl = n_af_clip_ctrl as f64 * 2.0 / (f_lcov_ctrl + f_rcov_ctrl);
            let f_polya_ratio_ctrl = n_polya_ctrl as f64 / (f_lcov_ctrl + f_rcov_ctrl);

            if b_trsdct {
                if fdiv(f_disc_ratio_ctrl, f_disc_ratio_case)? > f_case_ctrl_cutoff
                    && fdiv(f_clip_ratio_ctrl, f_clip_ratio_case)? > f_case_ctrl_cutoff
                {
                    writeln!(fout_log, "{}:{} (transduction) is filtered out as event have disc and clip in ctrl!", ins_chrm, ins_pos)?;
                    continue;
                }
            } else if fdiv(f_disc_ratio_ctrl, f_disc_ratio_case)? > f_case_ctrl_cutoff
                && fdiv(f_ef_clip_ratio_ctrl, f_ef_clip_ratio_case)? > f_case_ctrl_cutoff
            {
                writeln!(fout_log, "{}:{} is filtered out as event have disc and clip in ctrl!", ins_chrm, ins_pos)?;
                continue;
            }
            if f_polya_ratio_ctrl > f_ctrl_polya_cutoff {
                writeln!(fout_log, "{}:{} is filtered out at polyA ratio step!", ins_chrm, ins_pos)?;
                continue;
            }
            if ((n_af_clip_ctrl + n_full_map_ctrl) == 0
                || n_af_clip_ctrl as f64 / (n_af_clip_ctrl + n_full_map_ctrl) as f64 > f_ctrl_clip_ratio)
                && f_disc_ratio_ctrl > f_disc_cutoff
            {
                writeln!(fout_log, "{}:{} is filtered out at ctrl clip/full-map ratio step!", ins_chrm, ins_pos)?;
                continue;
            }
            let n_raw_clip_ctrl = n_raw_lclip_ctrl + n_raw_rclip_ctrl;
            let mut f_raw_clip_ratio_ctrl = 0.0;
            if n_raw_clip_ctrl + n_full_map_ctrl > 0 {
                f_raw_clip_ratio_ctrl = n_raw_clip_ctrl as f64 / (n_raw_clip_ctrl + n_full_map_ctrl) as f64;
            }
            let mut f_ctrl_polya_clip_ratio = 0.0;
            if n_raw_clip_ctrl > 0 {
                f_ctrl_polya_clip_ratio = n_polya_ctrl as f64 / n_raw_clip_ctrl as f64;
            }
            let b_ctrl_clip_polya_dominant = f_ctrl_polya_clip_ratio > f_ctrl_clip_ratio / 2.0;
            if ((n_raw_clip_ctrl + n_full_map_ctrl) == 0 || f_raw_clip_ratio_ctrl > f_ctrl_clip_ratio)
                && (f_disc_ratio_ctrl > f_disc_cutoff || (b_ctrl_clip_polya_dominant && n_disc_pairs_ctrl > 0))
            {
                writeln!(fout_log, "{}:{} is filtered out at ctrl raw clip step!", ins_chrm, ins_pos)?;
                continue;
            }
            if n_lpolya + n_rpolya == 0 && f_raw_clip_ratio_ctrl > f_ctrl_raw_clip_flex_ratio {
                writeln!(fout_log, "{}:{} is filtered out at case_polyA and ctrl raw clip step2!", ins_chrm, ins_pos)?;
                continue;
            }
            writeln!(fout_new, "{}", utils::rstrip(line))?;
        }
        fout_new.flush()?;
        fout_log.flush()?;
        Ok(())
    }

    /// `pass_cns_chk`: false when the control reads form the same clip / disc clusters (or
    /// have polyA on both sides).
    fn pass_cns_chk(&self, m_cluster: &SiteMap<TeiInfo>, m_clip_cluster: &SiteMap<ClipChk>, l_case: &[&str], n_polya_cutoff: i64) -> Result<bool> {
        let ins_chrm = fld(l_case, 0)?;
        let ins_pos = py_int(fld(l_case, 1)?)?;
        let s_case_lclip_clst = fld(l_case, 19)?;
        let s_case_rclip_clst = fld(l_case, 20)?;
        let s_case_ldisc_clst = fld(l_case, 21)?;
        let s_case_rdisc_clst = fld(l_case, 22)?;
        if let Some(c) = m_clip_cluster.get(ins_chrm).and_then(|m| m.get(&ins_pos)) {
            let s_ctrl_lclip = format!("{}:{}", c.lc_start, c.lc_end);
            let s_ctrl_rclip = format!("{}:{}", c.rc_start, c.rc_end);
            if self.is_clip_cluster_consist(s_case_lclip_clst, s_case_rclip_clst, &s_ctrl_lclip, &s_ctrl_rclip)? {
                return Ok(false);
            }
        }
        let rec = match m_cluster.get(ins_chrm).and_then(|m| m.get(&ins_pos)) {
            Some(r) => r,
            None => {
                println!("{}:{} doesn't form clip and disc cluster!", ins_chrm, ins_pos);
                return Ok(true);
            }
        };
        let (nlpolya, nrpolya) = (rec.4, rec.5);
        if nlpolya > n_polya_cutoff && nrpolya > n_polya_cutoff {
            return Ok(false);
        }
        let b_clip_consist = self.is_clip_cluster_consist(s_case_lclip_clst, s_case_rclip_clst, &rec.6, &rec.7)?;
        let b_disc_consist = self.is_disc_cluster_consist(s_case_ldisc_clst, s_case_rdisc_clst, &rec.8, &rec.9)?;
        Ok(!(b_clip_consist || b_disc_consist))
    }

    /// `pass_td_chk`: false when the control has the same transduction source or polyA tails.
    fn pass_td_chk(&self, m_td_cluster: &SiteMap<(String, i64)>, m_td_polya: &SiteMap<[i64; 6]>, l_fields: &[&str], n_polya_cutoff: i64) -> Result<bool> {
        let ins_chrm = fld(l_fields, 0)?;
        let ins_pos = py_int(fld(l_fields, 1)?)?;
        let s_case_src = fld(l_fields, 23)?;
        let s_ctrl_src = match m_td_cluster.get(ins_chrm).and_then(|m| m.get(&ins_pos)) {
            None => return Ok(true),
            Some(r) => r.0.as_str(),
        };
        if s_case_src.contains(s_ctrl_src) || s_ctrl_src == s_case_src || s_ctrl_src.contains(s_case_src) {
            println!("{}:{} transduction is filtered out, as it has the same source as control", ins_chrm, ins_pos);
            return Ok(false);
        }
        if let Some(p) = m_td_polya.get(ins_chrm).and_then(|m| m.get(&ins_pos)) {
            let n_polya = p[0] + p[1];
            let n_polyt = p[2] + p[3];
            if n_polya >= n_polya_cutoff || n_polyt >= n_polya_cutoff {
                println!("{}:{} transduction is filtered out, as there are polyA tails found in control", ins_chrm, ins_pos);
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn split2(s: &str) -> Vec<&str> {
        s.split(':').collect()
    }

    /// `_is_clip_cluster_consist`
    fn is_clip_cluster_consist(&self, s_case_lclip: &str, s_case_rclip: &str, s_ctrl_lclip: &str, s_ctrl_rclip: &str) -> Result<bool> {
        let l_ctrl_lclip = Self::split2(s_ctrl_lclip);
        let l_ctrl_rclip = Self::split2(s_ctrl_rclip);
        let l_case_lclip = Self::split2(s_case_lclip);
        let l_case_rclip = Self::split2(s_case_rclip);
        let v = |l: &[&str], i: usize| -> Result<i64> { int_float(fld(l, i)?) };
        let d = self.clip_consist_dist;
        let mut b_l = false;
        let mut b_r = false;
        if s_ctrl_rclip != "-1:-1" && s_case_rclip != "-1:-1"
            && ((v(&l_case_rclip, 0)? - v(&l_ctrl_rclip, 0)?).abs() < d || (v(&l_case_rclip, 1)? - v(&l_ctrl_rclip, 1)?).abs() < d)
        {
            b_r = true;
        }
        if s_ctrl_lclip != "-1:-1" && s_case_lclip != "-1:-1"
            && ((v(&l_case_lclip, 0)? - v(&l_ctrl_lclip, 0)?).abs() < d || (v(&l_case_lclip, 1)? - v(&l_ctrl_lclip, 1)?).abs() < d)
        {
            b_l = true;
        }
        if s_ctrl_rclip != "-1:-1" && s_case_rclip != "-1:-1"
            && ((v(&l_case_rclip, 0)? - v(&l_ctrl_rclip, 1)?).abs() < d || (v(&l_case_rclip, 1)? - v(&l_ctrl_rclip, 0)?).abs() < d)
        {
            b_r = true;
        }
        if s_ctrl_lclip != "-1:-1" && s_case_lclip != "-1:-1"
            && ((v(&l_case_lclip, 0)? - v(&l_ctrl_lclip, 1)?).abs() < d || (v(&l_case_lclip, 1)? - v(&l_ctrl_lclip, 0)?).abs() < d)
        {
            b_l = true;
        }
        Ok(b_l || b_r)
    }

    /// `_is_disc_cluster_consist` (the last two checks compare the same ends twice, as in the
    /// original)
    fn is_disc_cluster_consist(&self, s_case_ldisc: &str, s_case_rdisc: &str, s_ctrl_ldisc: &str, s_ctrl_rdisc: &str) -> Result<bool> {
        let l_ctrl_ldisc = Self::split2(s_ctrl_ldisc);
        let l_ctrl_rdisc = Self::split2(s_ctrl_rdisc);
        let l_case_ldisc = Self::split2(s_case_ldisc);
        let l_case_rdisc = Self::split2(s_case_rdisc);
        let v = |l: &[&str], i: usize| -> Result<i64> { int_float(fld(l, i)?) };
        let d = self.disc_consist_dist;
        let mut b_l = false;
        let mut b_r = false;
        if s_ctrl_rdisc != "-1:-1" && s_case_rdisc != "-1:-1"
            && ((v(&l_case_rdisc, 0)? - v(&l_ctrl_rdisc, 0)?).abs() < d || (v(&l_case_rdisc, 1)? - v(&l_ctrl_rdisc, 1)?).abs() < d)
        {
            b_r = true;
        }
        if s_ctrl_ldisc != "-1:-1" && s_case_ldisc != "-1:-1"
            && ((v(&l_case_ldisc, 0)? - v(&l_ctrl_ldisc, 0)?).abs() < d || (v(&l_case_ldisc, 1)? - v(&l_ctrl_ldisc, 1)?).abs() < d)
        {
            b_l = true;
        }
        if s_ctrl_rdisc != "-1:-1" && s_case_rdisc != "-1:-1"
            && ((v(&l_case_rdisc, 0)? - v(&l_ctrl_rdisc, 1)?).abs() < d || (v(&l_case_rdisc, 1)? - v(&l_ctrl_rdisc, 1)?).abs() < d)
        {
            b_r = true;
        }
        if s_ctrl_ldisc != "-1:-1" && s_case_ldisc != "-1:-1"
            && ((v(&l_case_ldisc, 1)? - v(&l_ctrl_ldisc, 0)?).abs() < d || (v(&l_case_ldisc, 1)? - v(&l_ctrl_ldisc, 1)?).abs() < d)
        {
            b_l = true;
        }
        Ok(b_l || b_r)
    }

    /// `parse_high_confident_somatic`: keep the high-confidence calls that are in the raw
    /// somatic call set.
    pub fn parse_high_confident_somatic(&self, sf_hc_sites: &str, sf_raw_somatic: &str, sf_out: &str) -> Result<()> {
        let mut m_raw: HashSet<(String, i64)> = HashSet::new();
        for line in utils::read_raw_lines(sf_raw_somatic)? {
            let fields = split_ws(&line);
            m_raw.insert((fld(&fields, 0)?.to_string(), py_int(fld(&fields, 1)?)?));
        }
        let mut fout = utils::create(sf_out)?;
        for line in utils::read_raw_lines(sf_hc_sites)? {
            let fields = split_ws(&line);
            let key = (fld(&fields, 0)?.to_string(), py_int(fld(&fields, 1)?)?);
            if !m_raw.contains(&key) {
                continue;
            }
            fout.write_all(line.as_bytes())?;
        }
        fout.flush()?;
        Ok(())
    }
}

/// A value of the dicts built by `SomaticMEICaller`: the original line, or the int 1 that
/// `filter_by_given_list` stores for chromosomes absent from the database.
#[derive(Clone, Debug)]
pub enum SiteVal {
    Line(String),
    One,
}

/// The parts of `SomaticMEICaller` used by `MosaicCaller.filter_by_blacklist`.
pub struct SomaticMeiCaller;

impl SomaticMeiCaller {
    /// `process_chrm_name`: drop a leading "chr".
    pub fn process_chrm_name(schrm: &str) -> String {
        if schrm.len() > 3 && schrm.starts_with("chr") {
            schrm[3..].to_string()
        } else {
            schrm.to_string()
        }
    }

    /// `load_sites`: chrm (without "chr") -> pos -> line (first occurrence kept).
    pub fn load_sites(sf_ce: &str) -> Result<SiteMap<SiteVal>> {
        let mut m: SiteMap<SiteVal> = IndexMap::new();
        for line in utils::read_raw_lines(sf_ce)? {
            let fields = split_ws(&line);
            let chrm = Self::process_chrm_name(fld(&fields, 0)?);
            let pos = py_int(fld(&fields, 1)?)?;
            m.entry(chrm).or_default().entry(pos).or_insert(SiteVal::Line(line.clone()));
        }
        Ok(m)
    }

    /// `filter_by_given_list`: drop the sites within `islack` of a database site.
    pub fn filter_by_given_list(sf_db: &str, m_meis: &SiteMap<SiteVal>, islack: i64) -> Result<SiteMap<SiteVal>> {
        let m_db = Self::load_sites(sf_db)?;
        let mut m_candidates: SiteMap<SiteVal> = IndexMap::new();
        for (chrm, m) in m_meis.iter() {
            match m_db.get(chrm) {
                None => {
                    let e = m_candidates.entry(chrm.clone()).or_default();
                    for &pos in m.keys() {
                        e.insert(pos, SiteVal::One);
                    }
                }
                Some(db) => {
                    for (&pos, v) in m.iter() {
                        let b_exist = (-islack..islack).any(|i| db.contains_key(&(pos + i)));
                        if !b_exist {
                            m_candidates.entry(chrm.clone()).or_default().insert(pos, v.clone());
                        }
                    }
                }
            }
        }
        Ok(m_candidates)
    }

    /// `export_slcted_candidates` (writing the int 1 raises a TypeError in the original)
    pub fn export_slcted_candidates(m_meis: &SiteMap<SiteVal>, sf_out: &str) -> Result<()> {
        let mut fout = utils::create(sf_out)?;
        for m in m_meis.values() {
            for v in m.values() {
                match v {
                    SiteVal::Line(s) => fout.write_all(s.as_bytes())?,
                    SiteVal::One => {
                        fout.flush()?;
                        return Err(anyhow!("TypeError: write() argument must be str, not int"));
                    }
                }
            }
        }
        fout.flush()?;
        Ok(())
    }
}
