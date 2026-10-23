//! Gene annotation of the called insertions from a GFF3 file (port of
//! `x_gene_annotation.py`).

use crate::globals as gv;
use crate::intervaltree::IntervalTree;
use crate::utils::{self, split_ws};
use anyhow::{anyhow, Result};
use indexmap::IndexMap;
use std::collections::HashMap;
use std::io::Write;

/// (end_pos, b_rc, gene_id, gene_type, gene_name) of a gene, keyed by its (extended) start.
#[derive(Clone, Debug)]
struct GeneRcd {
    end_pos: i64,
    #[allow(dead_code)]
    b_rc: bool,
    gene_id: String,
    gene_type: String,
    gene_name: String,
}

pub struct Gff3 {
    pub sf_gff3: String,
    m_gene_annotation: IndexMap<String, IndexMap<i64, Vec<GeneRcd>>>,
    m_region_info: HashMap<String, IndexMap<String, (i64, i64)>>,
    pub b_with_chr: bool,
    m_interval_tree: HashMap<String, IntervalTree>,
}

fn py_int(s: &str) -> Result<i64> {
    crate::py::fmt::parse_py_int(s).ok_or_else(|| anyhow!("invalid literal for int() with base 10: '{}'", s))
}

impl Gff3 {
    pub fn new(sf_gff3: &str) -> Self {
        Gff3 {
            sf_gff3: sf_gff3.to_string(),
            m_gene_annotation: IndexMap::new(),
            m_region_info: HashMap::new(),
            b_with_chr: true,
            m_interval_tree: HashMap::new(),
        }
    }

    pub fn set_with_chr(&mut self, b: bool) {
        self.b_with_chr = b;
    }

    /// `_process_chrm_name2`: drop a leading "chr".
    fn process_chrm_name2(chrm: &str) -> String {
        if chrm.len() > 3 && chrm.starts_with("chr") {
            chrm[3..].to_string()
        } else {
            chrm.to_string()
        }
    }

    /// `load_gene_annotation_with_extnd`
    pub fn load_gene_annotation_with_extnd(&mut self, iextnd: i64) -> Result<()> {
        for line in utils::read_raw_lines(&self.sf_gff3)? {
            if utils::rstrip(&line).is_empty() {
                continue;
            }
            let fields = split_ws(&line);
            if fields.len() < 9 {
                continue;
            }
            if fields[0].starts_with('#') {
                continue;
            }
            let chrm = Self::process_chrm_name2(fields[0]);
            let ori_start_pos = py_int(fields[3])?;
            let ori_end_pos = py_int(fields[4])?;
            let mut start_pos = ori_start_pos - iextnd;
            let end_pos = ori_end_pos + iextnd;
            let region_type = fields[2];
            if region_type == "transcript" {
                continue;
            }
            let b_rc = fields[6] == "-";
            let mut gene_id = String::new();
            let mut gene_type = String::new();
            let mut region_id = String::new();
            let mut gene_name = String::new();
            for tmp_field in fields[8].split(';') {
                let tt: Vec<&str> = tmp_field.split('=').collect();
                let val = || -> Result<String> {
                    tt.get(1).map(|s| s.to_string()).ok_or_else(|| anyhow!("list index out of range"))
                };
                if tt[0] == "ID" {
                    region_id = val()?;
                }
                if tt[0] == "gene_id" {
                    gene_id = val()?;
                }
                if tt[0] == "gene_type" {
                    gene_type = val()?;
                }
                if tt[0] == "gene_name" {
                    gene_name = val()?;
                }
            }
            if region_type == "gene" {
                let m_chrm = self.m_gene_annotation.entry(chrm.clone()).or_default();
                if m_chrm.contains_key(&start_pos) {
                    // genes sharing a start position: shift the key a little
                    start_pos -= 1;
                    while m_chrm.contains_key(&start_pos) {
                        start_pos -= 1;
                    }
                }
                m_chrm.entry(start_pos).or_default().push(GeneRcd {
                    end_pos,
                    b_rc,
                    gene_id: gene_id.clone(),
                    gene_type: gene_type.clone(),
                    gene_name: gene_name.clone(),
                });
                let mut regions = IndexMap::new();
                regions.insert(format!("{}:{}", gv::UP_STREAM_REGION, gene_id), (start_pos, ori_start_pos));
                regions.insert(format!("{}:{}", gv::DOWN_STREAM_REGION, gene_id), (ori_end_pos, end_pos));
                self.m_region_info.insert(gene_id, regions);
            } else {
                self.m_region_info.entry(gene_id).or_default().insert(region_id, (ori_start_pos, ori_end_pos));
            }
        }
        Ok(())
    }

    /// `index_gene_annotation_interval_tree`
    pub fn index_gene_annotation_interval_tree(&mut self) {
        for (chrm, m) in self.m_gene_annotation.iter() {
            let mut tree = IntervalTree::new();
            for (start_pos, l) in m.iter() {
                tree.addi(*start_pos, l[0].end_pos);
            }
            self.m_interval_tree.insert(chrm.clone(), tree);
        }
    }

    /// `query_by_position`: start positions of the genes overlapping `pos`.
    pub fn query_by_position(&self, chrm1: &str, pos: i64) -> Vec<i64> {
        let chrm = Self::process_chrm_name2(chrm1);
        match self.m_interval_tree.get(&chrm) {
            None => Vec::new(),
            Some(t) => t.at(pos).into_iter().map(|(b, _)| b).collect(),
        }
    }

    /// `get_gene_id_type`: (gene_id, gene_type, region description) of each hit gene.
    pub fn get_gene_id_type(&self, chrm1: &str, pos: i64) -> Vec<(String, String, String)> {
        let chrm = Self::process_chrm_name2(chrm1);
        let mut l_rslt = Vec::new();
        for gene_istart in self.query_by_position(&chrm, pos) {
            let g = &self.m_gene_annotation[&chrm][&gene_istart][0];
            let mut s_region_type = "intron".to_string();
            if let Some(regions) = self.m_region_info.get(&g.gene_id) {
                for (region_type, &(istart, iend)) in regions.iter() {
                    if pos >= istart && pos <= iend {
                        s_region_type = format!("{}:{}", region_type, g.gene_name);
                    }
                }
            }
            if s_region_type == "intron" {
                s_region_type = format!("intron:{}:{}", g.gene_id, g.gene_name);
            }
            l_rslt.push((g.gene_id.clone(), g.gene_type.clone(), s_region_type));
        }
        l_rslt
    }

    /// `annotate_results`: append the gene information to each record.
    pub fn annotate_results(&self, sf_ori_rslt: &str, sf_out: &str) -> Result<()> {
        let lines = utils::read_raw_lines(sf_ori_rslt)?;
        let mut w = utils::create(sf_out)?;
        for line in lines.iter() {
            if line.starts_with('#') {
                writeln!(w, "{}\tgene_info", utils::rstrip(line))?;
                continue;
            }
            let fields = split_ws(line);
            let chrm = *fields.first().ok_or_else(|| anyhow!("list index out of range"))?;
            let pos = py_int(fields.get(1).ok_or_else(|| anyhow!("list index out of range"))?)?;
            let l_hits = self.get_gene_id_type(chrm, pos);
            let mut s_hit_genes = String::new();
            for rcd in l_hits.iter() {
                if !s_hit_genes.is_empty() {
                    s_hit_genes.push('/');
                }
                s_hit_genes.push_str(&rcd.2);
            }
            if s_hit_genes.is_empty() {
                s_hit_genes = gv::NON_GENE.to_string();
            }
            writeln!(w, "{}\t{}", utils::rstrip(line), s_hit_genes)?;
        }
        w.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gene_regions() {
        let dir = std::env::temp_dir().join(format!("xtea_gff_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let sf = dir.join("a.gff3");
        std::fs::write(
            &sf,
            "##gff-version 3\n\
chr1\tHAVANA\tgene\t10000\t20000\t.\t+\t.\tID=G1;gene_id=G1;gene_type=protein_coding;gene_name=AAA\n\
chr1\tHAVANA\ttranscript\t10000\t20000\t.\t+\t.\tID=T1;gene_id=G1\n\
chr1\tHAVANA\texon\t10000\t10500\t.\t+\t.\tID=exon:T1:1;gene_id=G1\n\
chr1\tHAVANA\tgene\t10000\t12000\t.\t-\t.\tID=G2;gene_id=G2;gene_type=lncRNA;gene_name=BBB\n",
        )
        .unwrap();
        let mut g = Gff3::new(sf.to_str().unwrap());
        g.load_gene_annotation_with_extnd(1500).unwrap();
        g.index_gene_annotation_interval_tree();
        // G2 shares the start of G1 and is shifted by one
        let hits = g.get_gene_id_type("chr1", 10200);
        let names: Vec<String> = hits.iter().map(|h| h.2.clone()).collect();
        assert!(names.contains(&"exon:T1:1:AAA".to_string()));
        assert!(names.contains(&"intron:G2:BBB".to_string()));
        let up = g.get_gene_id_type("1", 9000);
        assert!(up.iter().any(|h| h.2 == "up_stream:G1:AAA"));
        assert!(g.get_gene_id_type("chr1", 30000).is_empty());
        std::fs::remove_dir_all(&dir).ok();
    }
}
