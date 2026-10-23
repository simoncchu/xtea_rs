//! Port of `x_coverage.py` (`ReadDepth`): local read depth on the two sides of candidate sites.

use crate::alignments::BamInfo;
use crate::annotation::process_chrm_name;
use crate::bam::{AlignmentFile, PyAln, Record, BAM_CHARD_CLIP, BAM_CINS, BAM_CSOFT_CLIP};
use crate::intermediate_sites::SiteMap;
use crate::num::{add, PyNum};
use crate::utils::{self, split_ws};
use anyhow::Result;
use indexmap::IndexMap;
use std::io::Write;

pub struct ReadDepth {
    pub working_folder: String,
    pub n_jobs: usize,
    pub sf_reference: String,
}

/// `_get_map_end_pos`: all operations but S, H and I advance the position.
#[inline]
pub fn get_map_end_pos(l_cigar: &[(u32, u32)], start_pos: i64) -> i64 {
    let mut end_pos = start_pos;
    for &(t, l) in l_cigar {
        if t == BAM_CSOFT_CLIP || t == BAM_CHARD_CLIP || t == BAM_CINS {
            continue;
        }
        end_pos += l as i64;
    }
    end_pos
}

/// Number of integer positions `i` with `a <= i < b` and `lo <= i <= hi`
/// (what the per-base loops of the original accumulate).
#[inline]
fn overlap(a: i64, b: i64, lo: i64, hi: i64) -> i64 {
    let s = a.max(lo);
    let e = (b - 1).min(hi);
    if e >= s {
        e - s + 1
    } else {
        0
    }
}

impl ReadDepth {
    pub fn new(working_folder: &str, n_jobs: usize, sf_ref: &str) -> Self {
        ReadDepth { working_folder: working_folder.to_string(), n_jobs, sf_reference: sf_ref.to_string() }
    }

    fn pool(&self) -> rayon::ThreadPool {
        rayon::ThreadPoolBuilder::new().num_threads(self.n_jobs.max(1)).build().unwrap()
    }

    fn bams(sf_bam_list: &str) -> Result<Vec<String>> {
        let mut v = Vec::new();
        for line in utils::read_lines(sf_bam_list)? {
            let f = split_ws(&line);
            if let Some(b) = f.first() {
                v.push(b.to_string());
            }
        }
        Ok(v)
    }

    /// `_calc_depth_one_site` -> (flcov, frcov)
    pub fn calc_depth_one_site(&self, chrm: &str, insertion_pos: i64, search_win: i64, focal_win: i64, sf_bam: &str) -> Result<(f64, f64)> {
        let mut start_pos = insertion_pos - search_win;
        if start_pos <= 0 {
            start_pos = 1;
        }
        let end_pos = insertion_pos + search_win;
        let bam_info = BamInfo::new(sf_bam, &self.sf_reference);
        let b_with_chr = bam_info.is_chrm_contain_chr()?;
        let chrm_in_bam = process_chrm_name(chrm, b_with_chr);
        let mut samfile = AlignmentFile::open(sf_bam, &self.sf_reference)?;
        let (mut nl, mut nr) = (0i64, 0i64);
        samfile.fetch(&chrm_in_bam, start_pos, end_pos, |a: &Record| {
            if a.is_dup() || a.is_unmap() {
                return;
            }
            let map_pos = a.reference_start();
            let tmp_end_pos = get_map_end_pos(&a.cigar_tuples(), map_pos);
            nl += overlap(map_pos, tmp_end_pos, insertion_pos - focal_win, insertion_pos);
            nr += overlap(map_pos, tmp_end_pos, insertion_pos + 1, insertion_pos + focal_win);
        })?;
        Ok((nl as f64 / focal_win as f64, nr as f64 / focal_win as f64))
    }

    /// `_calc_depth_one_site2` -> (flcov, frcov, flcov2, frcov2); integer zeros when the site is
    /// within `search_win` of the chromosome start (as in the original).
    pub fn calc_depth_one_site2(&self, chrm: &str, insertion_pos: i64, search_win: i64, focal_win: i64, focal_win2: i64, sf_bam: &str) -> Result<[PyNum; 4]> {
        let start_pos = insertion_pos - search_win;
        if start_pos <= 0 {
            return Ok([PyNum::Int(0); 4]);
        }
        let end_pos = insertion_pos + search_win;
        let bam_info = BamInfo::new(sf_bam, &self.sf_reference);
        let b_with_chr = bam_info.is_chrm_contain_chr()?;
        let chrm_in_bam = process_chrm_name(chrm, b_with_chr);
        let mut samfile = AlignmentFile::open(sf_bam, &self.sf_reference)?;
        let mut n = [0i64; 4];
        samfile.fetch(&chrm_in_bam, start_pos, end_pos, |a: &Record| {
            if a.is_dup() || a.is_unmap() {
                return;
            }
            let map_pos = a.reference_start();
            let tmp_end_pos = get_map_end_pos(&a.cigar_tuples(), map_pos);
            n[0] += overlap(map_pos, tmp_end_pos, insertion_pos - focal_win, insertion_pos);
            n[1] += overlap(map_pos, tmp_end_pos, insertion_pos + 1, insertion_pos + focal_win);
            n[2] += overlap(map_pos, tmp_end_pos, insertion_pos - focal_win2, insertion_pos);
            n[3] += overlap(map_pos, tmp_end_pos, insertion_pos + 1, insertion_pos + focal_win2);
        })?;
        Ok([
            PyNum::Float(n[0] as f64 / focal_win as f64),
            PyNum::Float(n[1] as f64 / focal_win as f64),
            PyNum::Float(n[2] as f64 / focal_win2 as f64),
            PyNum::Float(n[3] as f64 / focal_win2 as f64),
        ])
    }

    /// `calc_coverage_of_sites`: summed (over BAMs) (flcov, frcov) for the sites of `sf_sites`.
    pub fn calc_coverage_of_sites(&self, sf_sites: &str, sf_bam_list: &str, search_win: i64, focal_win: i64) -> Result<SiteMap<[PyNum; 2]>> {
        let mut m_site_cov: SiteMap<[PyNum; 2]> = IndexMap::new();
        let mut l_sites = Vec::new();
        for line in utils::read_lines(sf_sites)? {
            let f = split_ws(&line);
            if f.len() < 2 {
                continue;
            }
            l_sites.push((f[0].to_string(), utils::pint(f[1])?));
        }
        // NOTE: the original reads the BAM list with `line.rstrip()` (whole line)
        for line in utils::read_lines(sf_bam_list)? {
            let sf_bam = utils::rstrip(&line).to_string();
            let res: Vec<Result<(f64, f64)>> = self
                .pool()
                .install(|| utils::par_map(l_sites.len(), |i| self.calc_depth_one_site(&l_sites[i].0, l_sites[i].1, search_win, focal_win, &sf_bam)));
            for (i, r) in res.into_iter().enumerate() {
                let (fl, fr) = r?;
                let e = m_site_cov.entry(l_sites[i].0.clone()).or_default().entry(l_sites[i].1).or_insert([PyNum::Int(0); 2]);
                e[0] = add(e[0], PyNum::Float(fl));
                e[1] = add(e[1], PyNum::Float(fr));
            }
        }
        Ok(m_site_cov)
    }

    /// `calc_coverage_of_two_regions`
    pub fn calc_coverage_of_two_regions(&self, l_sites: &[(String, i64)], sf_bam_list: &str, search_win: i64, focal_win: i64, focal_win2: i64) -> Result<SiteMap<[PyNum; 4]>> {
        let mut m_site_cov: SiteMap<[PyNum; 4]> = IndexMap::new();
        for sf_bam in Self::bams(sf_bam_list)? {
            let res: Vec<Result<[PyNum; 4]>> = self.pool().install(|| {
                utils::par_map(l_sites.len(), |i| self.calc_depth_one_site2(&l_sites[i].0, l_sites[i].1, search_win, focal_win, focal_win2, &sf_bam))
            });
            for (i, r) in res.into_iter().enumerate() {
                let v = r?;
                let e = m_site_cov.entry(l_sites[i].0.clone()).or_default().entry(l_sites[i].1).or_insert([PyNum::Int(0); 4]);
                for k in 0..4 {
                    e[k] = add(e[k], v[k]);
                }
            }
        }
        Ok(m_site_cov)
    }

    /// `dump_coverage_info`
    pub fn dump_coverage_info(m_site_cov: &SiteMap<[PyNum; 4]>, sf_out: &str) -> Result<()> {
        let mut fout = utils::create(sf_out)?;
        for (ins_chrm, m) in m_site_cov.iter() {
            for (ins_pos, c) in m.iter() {
                writeln!(fout, "{}\t{}\t{}\t{}\t{}\t{}", ins_chrm, ins_pos, c[0], c[1], c[2], c[3])?;
            }
        }
        fout.flush()?;
        Ok(())
    }

    /// `calc_coverage_of_sites2`: coverage of the given sites in one BAM, skipping those below
    /// `MIN_COV_RANDOM_SITE`.
    pub fn calc_coverage_of_sites2(&self, l_sites: &[(String, i64)], sf_bam: &str, search_win: i64, focal_win: i64) -> Result<Vec<f64>> {
        let res: Vec<Result<(f64, f64)>> = self
            .pool()
            .install(|| utils::par_map(l_sites.len(), |i| self.calc_depth_one_site(&l_sites[i].0, l_sites[i].1, search_win, focal_win, sf_bam)));
        let mut v = Vec::new();
        for r in res {
            let (fl, fr) = r?;
            let f_cov = (fl + fr) / 2.0;
            if f_cov < crate::globals::MIN_COV_RANDOM_SITE {
                continue;
            }
            v.push(f_cov);
        }
        Ok(v)
    }
}
