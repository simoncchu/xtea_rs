//! IGV batch-script generation for screenshots of called sites (port of `x_igv.py`).

use crate::py::fmt::{parse_py_int, py_str_repr};
use crate::utils::{self, split_ws};
use anyhow::{anyhow, Result};
use indexmap::IndexMap;
use std::io::Write;

#[derive(Default)]
pub struct XIgv;

fn pint(s: &str) -> Result<i64> {
    parse_py_int(s).ok_or_else(|| anyhow!("invalid literal for int() with base 10: '{}'", s))
}

impl XIgv {
    pub fn new() -> Self {
        XIgv
    }

    fn load_sites(sf_sites: &str, i_extnd: i64, with_pos_in_id: bool) -> Result<IndexMap<String, Vec<(String, String)>>> {
        let mut m_sites: IndexMap<String, Vec<(String, String)>> = IndexMap::new();
        for line in utils::read_raw_lines(sf_sites)? {
            let fields = split_ws(&line);
            if fields.len() < 3 {
                println!("({}, {})", py_str_repr("Wrong site "), py_str_repr(&line));
                continue;
            }
            let s_id = fields[0];
            let ins_chrm = fields[1];
            let ins_pos = pint(fields[2])?;
            let s_region = format!("{}:{}-{}", ins_chrm, ins_pos - i_extnd, ins_pos + i_extnd);
            let s_region_id = if with_pos_in_id {
                format!("{}_{}_{}_{}", ins_chrm, ins_pos, ins_pos - i_extnd, ins_pos + i_extnd)
            } else {
                format!("{}_{}_{}", ins_chrm, ins_pos - i_extnd, ins_pos + i_extnd)
            };
            m_sites.entry(s_id.to_string()).or_default().push((s_region, s_region_id));
        }
        Ok(m_sites)
    }

    /// `prepare_igv_scripts_multi_bams`: each line of `sf_bams` is "sample-id bam1 [bam2 ...]".
    pub fn prepare_igv_scripts_multi_bams(&self, sf_sites: &str, sf_bams: &str, s_wfolder: &str, i_extnd: i64, sf_gnm: &str, sf_out: &str) -> Result<()> {
        let mut m_bams: IndexMap<String, Vec<String>> = IndexMap::new();
        for line in utils::read_raw_lines(sf_bams)? {
            let fields = split_ws(&line);
            if fields.len() < 2 {
                println!("({}, {})", py_str_repr("Wrong bam file list "), py_str_repr(&line));
                continue;
            }
            m_bams.insert(fields[0].to_string(), fields[1..].iter().map(|s| s.to_string()).collect());
        }
        let m_sites = Self::load_sites(sf_sites, i_extnd, true)?;
        let mut fout = utils::create(sf_out)?;
        for (s_id, l_sites) in m_sites.iter() {
            let l_sf_bam = match m_bams.get(s_id) {
                None => {
                    println!("({}, {})", py_str_repr(s_id), py_str_repr("not in provided bam list!"));
                    continue;
                }
                Some(b) => b,
            };
            self.new_bams_cmd_to_file(&mut fout, l_sf_bam, sf_gnm, s_wfolder)?;
            for (s_region, s_region_id) in l_sites.iter() {
                let s_screenshot_id = format!("{}.{}.png", s_id, s_region_id);
                self.new_site_cmd_to_file(&mut fout, s_region, &s_screenshot_id)?;
            }
        }
        fout.flush()?;
        Ok(())
    }

    /// `gnrt_sites_single_sample`: prefix every site with the (last) sample id of the BAM
    /// list; returns the new sites file.
    pub fn gnrt_sites_single_sample(&self, sf_sites: &str, sf_bams: &str) -> Result<String> {
        let mut s_id = "null".to_string();
        for line in utils::read_raw_lines(sf_bams)? {
            let fields = split_ws(&line);
            s_id = fields.first().ok_or_else(|| anyhow!("list index out of range"))?.to_string();
        }
        let sf_new_sites = format!("{}.igv_sites", sf_sites);
        let lines = utils::read_raw_lines(sf_sites)?;
        let mut w = utils::create(&sf_new_sites)?;
        for line in lines.iter() {
            write!(w, "{}\t{}", s_id, line)?;
        }
        w.flush()?;
        Ok(sf_new_sites)
    }

    fn new_bams_cmd_to_file<W: Write>(&self, fout: &mut W, l_sf_bam: &[String], sf_gnm: &str, s_folder: &str) -> Result<()> {
        writeln!(fout, "new")?;
        writeln!(fout, "genome {}", sf_gnm)?;
        for sf_bam in l_sf_bam {
            writeln!(fout, "load {}", sf_bam)?;
        }
        writeln!(fout, "snapshotDirectory {}", s_folder)?;
        Ok(())
    }

    fn new_bam_cmd_to_file<W: Write>(&self, fout: &mut W, sf_bam: &str, sf_gnm: &str, s_folder: &str) -> Result<()> {
        writeln!(fout, "new")?;
        writeln!(fout, "genome {}", sf_gnm)?;
        writeln!(fout, "load {}", sf_bam)?;
        writeln!(fout, "snapshotDirectory {}", s_folder)?;
        Ok(())
    }

    fn new_site_cmd_to_file<W: Write>(&self, fout: &mut W, s_region: &str, s_name: &str) -> Result<()> {
        writeln!(fout, "goto {}", s_region)?;
        writeln!(fout, "sort position")?;
        writeln!(fout, "collapse")?;
        writeln!(fout, "snapshot {}", s_name)?;
        Ok(())
    }

    /// `prepare_igv_scripts`: one BAM per sample ("sample-id bam" per line).
    pub fn prepare_igv_scripts(&self, sf_sites: &str, sf_bams: &str, s_wfolder: &str, i_extnd: i64, sf_gnm: &str, sf_out: &str) -> Result<()> {
        let mut m_bams: IndexMap<String, String> = IndexMap::new();
        for line in utils::read_raw_lines(sf_bams)? {
            let fields = split_ws(&line);
            if fields.len() < 2 {
                println!("({}, {})", py_str_repr("Wrong bam file list "), py_str_repr(&line));
                continue;
            }
            m_bams.insert(fields[0].to_string(), fields[1].to_string());
        }
        let m_sites = Self::load_sites(sf_sites, i_extnd, false)?;
        let mut fout = utils::create(sf_out)?;
        for (s_id, l_sites) in m_sites.iter() {
            let sf_bam = match m_bams.get(s_id) {
                None => {
                    println!("({}, {})", py_str_repr(s_id), py_str_repr("not in provided bam list!"));
                    continue;
                }
                Some(b) => b,
            };
            self.new_bam_cmd_to_file(&mut fout, sf_bam, sf_gnm, s_wfolder)?;
            for (s_region, s_region_id) in l_sites.iter() {
                let s_screenshot_id = format!("{}.{}.png", s_id, s_region_id);
                self.new_site_cmd_to_file(&mut fout, s_region, &s_screenshot_id)?;
            }
        }
        fout.flush()?;
        Ok(())
    }
}
