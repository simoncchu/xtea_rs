//! Port of `x_basic_info.py`: coverage, read length and insert size estimation from random sites.

use crate::alignments::BamInfo;
use crate::bam::{AlignmentFile, PyAln};
use crate::chrom::is_decoy_contig_chrms;
use crate::globals as gv;
use crate::num::PyNum;
use crate::rng::with_rng;
use crate::utils::{self, split_ws};
use anyhow::Result;
use indexmap::IndexMap;
use std::io::Write;

/// Result of `_collect_basic_info_one_site`.
#[derive(Clone, Debug)]
pub struct SiteInfo {
    pub chrm: String,
    pub pos: i64,
    pub flcov: f64,
    pub frcov: f64,
    pub dom_rlth: i64,
    pub acm_is: i64,
    pub acm_is_squre: i64,
    pub n_pairs: i64,
}

/// (f_ave_cov, rlth, mean_is, std_var) as stored per BAM.
#[derive(Clone, Debug)]
pub struct SampleInfo {
    pub cov: PyNum,
    pub rlth: i64,
    pub mean_is: PyNum,
    pub std_var: PyNum,
}

/// The tuple returned by `get_cov_is_rlth`.
#[derive(Clone, Copy, Debug)]
pub struct BasicRcd {
    pub cov: f64,
    pub rlth: f64,
    pub mean_is: f64,
    pub std_var: f64,
    /// true when the original returned the integer tuple (0, 0, 0, 0)
    pub zero: bool,
}

pub struct XBasicInfo {
    pub working_folder: String,
    pub n_jobs: usize,
    pub sf_reference: String,
}

impl XBasicInfo {
    pub fn new(working_folder: &str, n_jobs: usize, sf_ref: &str) -> Self {
        let mut wf = if working_folder.is_empty() { "./".to_string() } else { working_folder.to_string() };
        if !wf.ends_with('/') {
            wf.push('/');
        }
        XBasicInfo { working_folder: wf, n_jobs, sf_reference: sf_ref.to_string() }
    }

    /// `get_cov_is_rlth`
    pub fn get_cov_is_rlth(&self, sf_bam_list: &str, sf_ref: &str, search_win: i64, b_force: bool) -> Result<BasicRcd> {
        let sf_basic_info = format!("{}{}", self.working_folder, gv::BASIC_INFO_FILE);
        let mut m_info: IndexMap<String, SampleInfo> = IndexMap::new();
        let mut b_gnrt = true;
        if utils::is_file(&sf_basic_info) && !b_force {
            let m_tmp = self.load_basic_info_from_file(&sf_basic_info)?;
            if !m_tmp.is_empty() {
                m_info = m_tmp;
                b_gnrt = false;
            }
        }
        if b_gnrt {
            m_info = self.collect_dump_basic_info_samples(sf_bam_list, sf_ref, search_win)?;
        }
        let mut f_cov = PyNum::Int(0);
        let mut rlth: i64 = 0;
        let mut mean_is = PyNum::Int(0);
        let mut std_var = PyNum::Int(0);
        for (_, v) in m_info.iter() {
            f_cov = crate::num::add(f_cov, v.cov);
            rlth += v.rlth;
            mean_is = crate::num::add(mean_is, v.mean_is);
            std_var = crate::num::add(std_var, v.std_var);
        }
        let n = m_info.len();
        if n == 0 {
            return Ok(BasicRcd { cov: 0.0, rlth: 0.0, mean_is: 0.0, std_var: 0.0, zero: true });
        }
        Ok(BasicRcd {
            cov: f_cov.f(),
            rlth: rlth as f64 / n as f64,
            mean_is: mean_is.f() / n as f64,
            std_var: std_var.f() / n as f64,
            zero: false,
        })
    }

    pub fn collect_dump_basic_info_samples(&self, sf_bam_list: &str, sf_ref: &str, search_win: i64) -> Result<IndexMap<String, SampleInfo>> {
        let m = self.collect_basic_info_samples(sf_bam_list, sf_ref, search_win)?;
        let sf_out = format!("{}{}", self.working_folder, gv::BASIC_INFO_FILE);
        self.dump_basic_info_to_file(&m, &sf_out)?;
        Ok(m)
    }

    pub fn load_basic_info_samples(&self) -> Result<Option<IndexMap<String, SampleInfo>>> {
        let sf = format!("{}{}", self.working_folder, gv::BASIC_INFO_FILE);
        if !utils::is_file(&sf) {
            return Ok(None);
        }
        Ok(Some(self.load_basic_info_from_file(&sf)?))
    }

    /// `load_get_cov` -> (sum of coverage, average coverage); (0, 0) if no file.
    pub fn load_get_cov(&self) -> Result<(f64, f64)> {
        let m = match self.load_basic_info_samples()? {
            None => return Ok((0.0, 0.0)),
            Some(m) => m,
        };
        let mut f_cov = 0.0;
        for (_, v) in m.iter() {
            f_cov += v.cov.f();
        }
        let n = m.len();
        Ok((f_cov, f_cov / n as f64))
    }

    /// `calc_cov_from_bam` (median coverage of random sites).
    pub fn calc_cov_from_bam(&self, sf_bam: &str, sf_ref: &str, search_win: i64) -> Result<f64> {
        let l_sites = self.random_slct_site(sf_bam, sf_ref, gv::N_RANDOM_SITES)?;
        let infos = self.collect_basic_info_of_sites2(&l_sites, sf_bam, search_win, search_win);
        let mut l_cov = Vec::new();
        for r in infos.iter().flatten() {
            let f_cov = (r.flcov + r.frcov) / 2.0;
            if f_cov >= gv::MIN_COV_RANDOM_SITE {
                l_cov.push(f_cov);
            }
        }
        l_cov.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if l_cov.is_empty() {
            return Ok(0.0);
        }
        Ok(l_cov[l_cov.len() / 2])
    }

    /// `collect_basic_info_samples`
    pub fn collect_basic_info_samples(&self, sf_bam_list: &str, sf_ref: &str, search_win: i64) -> Result<IndexMap<String, SampleInfo>> {
        let mut m_sample_info: IndexMap<String, SampleInfo> = IndexMap::new();
        for line in utils::read_lines(sf_bam_list)? {
            let fields = split_ws(&line);
            if fields.is_empty() {
                continue;
            }
            let sf_bam = fields[0].to_string();
            let l_sites = self.random_slct_site(&sf_bam, sf_ref, gv::N_RANDOM_SITES)?;
            let infos = self.collect_basic_info_of_sites2(&l_sites, &sf_bam, search_win, search_win);
            let mut acm_is: i64 = 0;
            let mut acm_is_squre: i64 = 0;
            let mut n_pairs: i64 = 0;
            let mut m_rlth: IndexMap<i64, i64> = IndexMap::new();
            let mut l_cov: Vec<f64> = Vec::new();
            for r in infos.into_iter() {
                let r = r?;
                acm_is += r.acm_is;
                acm_is_squre += r.acm_is_squre;
                n_pairs += r.n_pairs;
                let f_cov = (r.flcov + r.frcov) / 2.0;
                if f_cov >= gv::MIN_COV_RANDOM_SITE {
                    l_cov.push(f_cov);
                }
                *m_rlth.entry(r.dom_rlth).or_insert(0) += 1;
            }
            l_cov.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let f_ave_cov = if !l_cov.is_empty() { PyNum::Float(l_cov[l_cov.len() / 2]) } else { PyNum::Int(0) };
            let rlth = calc_read_length(&m_rlth);
            let (mean_is, std_var) = calc_insert_size(acm_is, acm_is_squre, n_pairs);
            m_sample_info.insert(sf_bam, SampleInfo { cov: f_ave_cov, rlth, mean_is, std_var });
        }
        Ok(m_sample_info)
    }

    /// `collect_basic_info_one_sample` -> (cov, rlth, mean_is, std)
    pub fn collect_basic_info_one_sample(&self, sf_bam: &str, sf_ref: &str, search_win: i64) -> Result<SampleInfo> {
        let l_sites = self.random_slct_site(sf_bam, sf_ref, gv::N_RANDOM_SITES)?;
        let infos = self.collect_basic_info_of_sites2(&l_sites, sf_bam, search_win, search_win);
        let mut acm_is = 0;
        let mut acm_is_squre = 0;
        let mut n_pairs = 0;
        let mut m_rlth: IndexMap<i64, i64> = IndexMap::new();
        let mut l_cov = Vec::new();
        for r in infos.into_iter() {
            let r = r?;
            acm_is += r.acm_is;
            acm_is_squre += r.acm_is_squre;
            n_pairs += r.n_pairs;
            let f_cov = (r.flcov + r.frcov) / 2.0;
            if f_cov >= gv::MIN_COV_RANDOM_SITE {
                l_cov.push(f_cov);
            }
            *m_rlth.entry(r.dom_rlth).or_insert(0) += 1;
        }
        l_cov.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let cov = if !l_cov.is_empty() { PyNum::Float(l_cov[l_cov.len() / 2]) } else { PyNum::Int(0) };
        let rlth = calc_read_length(&m_rlth);
        let (mean_is, std_var) = calc_insert_size(acm_is, acm_is_squre, n_pairs);
        Ok(SampleInfo { cov, rlth, mean_is, std_var })
    }

    /// `dump_basic_info_to_file`
    pub fn dump_basic_info_to_file(&self, m: &IndexMap<String, SampleInfo>, sf_out: &str) -> Result<()> {
        let mut f = utils::create(sf_out)?;
        for (sbam, rcd) in m.iter() {
            write!(f, "{}", sbam)?;
            write!(f, "\t{}", rcd.cov.py_str())?;
            write!(f, "\t{}", rcd.rlth)?;
            write!(f, "\t{}", rcd.mean_is.py_str())?;
            write!(f, "\t{}", rcd.std_var.py_str())?;
            writeln!(f)?;
        }
        f.flush()?;
        Ok(())
    }

    /// `load_basic_info_from_file`
    pub fn load_basic_info_from_file(&self, sf: &str) -> Result<IndexMap<String, SampleInfo>> {
        let mut m = IndexMap::new();
        for line in utils::read_lines(sf)? {
            let fields = split_ws(&line);
            if fields.len() < 5 {
                continue;
            }
            let cov = utils::pfloat(fields[1])?;
            let rlth = utils::pint(fields[2])?;
            let mean_is = utils::pfloat(fields[3])?;
            let std_var = utils::pfloat(fields[4])?;
            m.insert(
                fields[0].to_string(),
                SampleInfo { cov: PyNum::Float(cov), rlth, mean_is: PyNum::Float(mean_is), std_var: PyNum::Float(std_var) },
            );
        }
        Ok(m)
    }

    /// `collect_basic_info_of_sites2`: run `_collect_basic_info_one_site` for all sites in
    /// parallel, results in order.
    pub fn collect_basic_info_of_sites2(&self, l_sites: &[(String, i64)], sf_bam: &str, search_win: i64, focal_win: i64) -> Vec<Result<SiteInfo>> {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(self.n_jobs.max(1)).build().unwrap();
        pool.install(|| {
            utils::par_map(l_sites.len(), |i| {
                let (chrm, pos) = &l_sites[i];
                collect_basic_info_one_site(chrm, *pos, search_win, focal_win, sf_bam, &self.sf_reference)
            })
        })
    }

    /// `_random_slct_site`
    pub fn random_slct_site(&self, sf_bam: &str, sf_ref: &str, n_sites: i64) -> Result<Vec<(String, i64)>> {
        let bam_info = BamInfo::new(sf_bam, sf_ref);
        let m_chrm_info = bam_info.get_all_chrom_name_length()?;
        let mut l_chrm = Vec::new();
        let mut l_lenth = Vec::new();
        for (chrm, len) in m_chrm_info.iter() {
            if is_decoy_contig_chrms(chrm) {
                continue;
            }
            l_chrm.push(chrm.clone());
            l_lenth.push(*len);
        }
        let n_chrm = l_chrm.len() as i64;
        if n_chrm == 0 {
            return Err(anyhow::anyhow!("empty range for randrange() (no usable chromosome in {})", sf_bam));
        }
        let mut l_sites = Vec::new();
        with_rng(|rng| {
            for _ in 0..n_sites {
                let i_tmp = rng.randint(0, n_chrm - 1) as usize;
                let len = l_lenth[i_tmp];
                if len < 2000000 {
                    continue;
                }
                let pos = rng.randint(1000000, len - 1000000);
                l_sites.push((l_chrm[i_tmp].clone(), pos));
            }
        });
        Ok(l_sites)
    }
}

/// `_get_map_end_pos`: start + lengths of all ops except S, H, I.
pub fn get_map_end_pos(cigar: &[(u32, u32)], start_pos: i64) -> i64 {
    let mut end = start_pos;
    for &(t, l) in cigar {
        if t == 4 || t == 5 || t == 1 {
            continue;
        }
        end += l as i64;
    }
    end
}

/// `_calc_read_length`: dominant read length (first maximum in insertion order).
pub fn calc_read_length(m: &IndexMap<i64, i64>) -> i64 {
    let mut dom = 0;
    let mut cnt = 0;
    for (k, v) in m.iter() {
        if *v > cnt {
            cnt = *v;
            dom = *k;
        }
    }
    dom
}

/// `_calc_insert_size`
pub fn calc_insert_size(acm_is: i64, acm_is_squre: i64, n_pairs: i64) -> (PyNum, PyNum) {
    if n_pairs <= 0 {
        return (PyNum::Int(0), PyNum::Int(0));
    }
    let mean_is = acm_is as f64 / n_pairs as f64;
    let e_is_squre = acm_is_squre as f64 / n_pairs as f64;
    let e_mean_squre = mean_is.powi(2);
    let v = e_is_squre - e_mean_squre;
    let std_var = v.powf(0.5);
    (PyNum::Float(mean_is), PyNum::Float(std_var))
}

/// `_collect_basic_info_one_site`
pub fn collect_basic_info_one_site(chrm: &str, insertion_pos: i64, search_win: i64, focal_win: i64, sf_bam: &str, sf_ref: &str) -> Result<SiteInfo> {
    let mut start_pos = insertion_pos - search_win;
    if start_pos <= 0 {
        start_pos = 1;
    }
    let end_pos = insertion_pos + search_win;
    let bam_info = BamInfo::new(sf_bam, sf_ref);
    let b_with_chr = bam_info.is_chrm_contain_chr()?;
    let chrm_in_bam = crate::annotation::process_chrm_name(chrm, b_with_chr);
    let mut samfile = AlignmentFile::open(sf_bam, sf_ref)?;
    let lfocal_start = insertion_pos - focal_win;
    let lfocal_end = insertion_pos;
    let rfocal_start = insertion_pos + 1;
    let rfocal_end = insertion_pos + focal_win;
    let mut lsum: i64 = 0;
    let mut rsum: i64 = 0;
    let mut l_insert_size: Vec<i64> = Vec::new();
    let mut m_read_lenth: IndexMap<i64, i64> = IndexMap::new();
    let n_refs = samfile.references().len() as i32;
    let names: Vec<String> = samfile.references().to_vec();
    samfile.fetch(&chrm_in_bam, start_pos, end_pos, |a| {
        if a.is_sec() || a.is_supp() || a.is_dup() || a.is_unmap() {
            return;
        }
        let map_pos = a.reference_start();
        let cigar = a.cigar_tuples();
        let tmp_end = get_map_end_pos(&cigar, map_pos);
        // overlap of [map_pos, tmp_end) with [lfocal_start, lfocal_end] and [rfocal_start, rfocal_end]
        let lo = map_pos.max(lfocal_start);
        let hi = (tmp_end - 1).min(lfocal_end);
        if hi >= lo {
            lsum += hi - lo + 1;
        }
        let lo = map_pos.max(rfocal_start);
        let hi = (tmp_end - 1).min(rfocal_end);
        if hi >= lo {
            rsum += hi - lo + 1;
        }
        if a.is_read1() {
            let ntid = a.next_reference_id();
            if ntid < -1 || ntid >= n_refs {
                return;
            }
            let (mate_chrm, mate_pos) = if ntid >= 0 { (names[ntid as usize].as_str(), a.next_reference_start()) } else { ("*", 0) };
            let is_tmp = if chrm != mate_chrm || mate_chrm == "*" {
                0
            } else {
                let dist = (map_pos - mate_pos).abs();
                if dist > gv::MAX_NORMAL_INSERT_SIZE {
                    0
                } else {
                    dist
                }
            };
            if is_tmp != 0 {
                l_insert_size.push(is_tmp);
            }
        }
        if cigar.len() > 1 && (cigar[0].0 == 5 || cigar[cigar.len() - 1].0 == 5) {
            return;
        }
        let rlen = a.query_length();
        if rlen != 0 {
            *m_read_lenth.entry(rlen).or_insert(0) += 1;
        }
    })?;
    let flcov = lsum as f64 / focal_win as f64;
    let frcov = rsum as f64 / focal_win as f64;
    let dom = calc_read_length(&m_read_lenth);
    let mut acm_sq: i64 = 0;
    let mut acm: i64 = 0;
    for is in &l_insert_size {
        let adj = is + dom;
        acm_sq += adj * adj;
        acm += adj;
    }
    Ok(SiteInfo {
        chrm: chrm.to_string(),
        pos: insertion_pos,
        flcov,
        frcov,
        dom_rlth: dom,
        acm_is: acm,
        acm_is_squre: acm_sq,
        n_pairs: l_insert_size.len() as i64,
    })
}
