//! Port of `bwa_align.py` (re-alignment of clipped/discordant reads with bwa / samtools).

use crate::bam::{AlignmentFile, AlignmentWriter, PyAln};
use crate::fasta::read_fastx;
use crate::globals as gv;
use crate::utils;
use anyhow::Result;
use std::io::Write;

pub struct BwAlign {
    pub bwa_path: String,
    pub bwa_realign_cutoff: i64,
    pub n_jobs: i64,
    pub bwa_seed_freq: i64,
    pub bwa_seed_medium_freq: i64,
}

impl BwAlign {
    pub fn new(bwa_path: &str, bwa_realign_cutoff: i64, n_jobs: i64) -> Self {
        BwAlign { bwa_path: bwa_path.to_string(), bwa_realign_cutoff, n_jobs, bwa_seed_freq: 10, bwa_seed_medium_freq: 70 }
    }

    /// `realign_clipped_reads`
    pub fn realign_clipped_reads(&self, sf_ref: &str, sf_reads: &str, sf_out_sam: &str) {
        let cmd = format!(
            "{} mem -t {} -T {} -k {} -o {} {} {}",
            self.bwa_path, self.n_jobs, self.bwa_realign_cutoff, self.bwa_realign_cutoff, sf_out_sam, sf_ref, sf_reads
        );
        utils::run_cmd_small_output(&cmd);
    }

    /// `realign_clipped_read_with_polyA`
    pub fn realign_clipped_read_with_polya(&self, sf_ref: &str, sf_reads: &str, sf_out_sam: &str) -> Result<()> {
        let sf_reads1 = format!("{}.non_polyA.fq", sf_reads);
        let sf_reads2 = format!("{}.polyA.fq", sf_reads);
        {
            let mut fout_non = utils::create(&sf_reads1)?;
            let mut fout_polya = utils::create(&sf_reads2)?;
            for entry in read_fastx(sf_reads)? {
                let len_seq = entry.sequence.len() as i64;
                if len_seq < gv::BWA_REALIGN_CUTOFF {
                    writeln!(fout_polya, "{}", entry.to_pysam_string())?;
                } else {
                    writeln!(fout_non, "{}", entry.to_pysam_string())?;
                }
            }
            fout_non.flush()?;
            fout_polya.flush()?;
        }
        let sf_sam1 = format!("{}.non_polyA.sam", sf_out_sam);
        self.realign_clipped_reads(sf_ref, &sf_reads1, &sf_sam1);
        let sf_sam2 = format!("{}.polyA.sam", sf_out_sam);
        self.realign_clipped_polya(sf_ref, &sf_reads2, &sf_sam2);
        {
            let mut fout = utils::create(sf_out_sam)?;
            for line in utils::read_raw_lines(&sf_sam1)? {
                fout.write_all(line.as_bytes())?;
            }
            for line in utils::read_raw_lines(&sf_sam2)? {
                if line.starts_with('@') {
                    continue;
                }
                fout.write_all(line.as_bytes())?;
            }
            fout.flush()?;
        }
        for f in [&sf_reads1, &sf_reads2, &sf_sam1, &sf_sam2] {
            utils::remove_file(f);
        }
        Ok(())
    }

    pub fn index_ref_file(&self, sf_fa: &str) {
        let cmd = format!("{} index {}", self.bwa_path, sf_fa);
        utils::run_cmd_small_output(&cmd);
    }

    /// `realign_clipped_polyA`
    pub fn realign_clipped_polya(&self, sf_ref: &str, sf_reads: &str, sf_out_sam: &str) {
        let cmd = format!(
            "{} mem -t {} -T {} -k {} -o {} -c {} {} {}",
            self.bwa_path,
            self.n_jobs,
            gv::MINIMUM_POLYA_CLIP,
            gv::MINIMUM_POLYA_CLIP,
            sf_out_sam,
            self.bwa_seed_medium_freq,
            sf_ref,
            sf_reads
        );
        utils::run_cmd_small_output(&cmd);
    }

    /// `realign_clipped_reads_low_mem`
    pub fn realign_clipped_reads_low_mem(&self, sf_ref: &str, sf_reads: &str, sf_out_sam: &str) {
        let mut n_cores = self.n_jobs - 1;
        if n_cores <= 0 {
            n_cores = 1;
        }
        let cmd = format!(
            "{} mem -t {} -T {} -k {} -c {} -o {} {} {}",
            self.bwa_path, n_cores, self.bwa_realign_cutoff, self.bwa_realign_cutoff, self.bwa_seed_freq, sf_out_sam, sf_ref, sf_reads
        );
        utils::run_cmd_to_file(&cmd, &format!("{}.std_out", sf_out_sam));
    }

    /// `realign_clipped_reads_low_mem_add_head`
    pub fn realign_clipped_reads_low_mem_add_head(&self, sf_ref: &str, sf_reads: &str, sf_head: &str, sf_out_sam: &str) {
        let mut n_cores = self.n_jobs - 1;
        if n_cores <= 0 {
            n_cores = 1;
        }
        let cmd = format!(
            "{} mem -t {} -T {} -k {} -c {} -D 0.9 -h 2 -H {} {} {} | {} view -hS -F 4 -o {} -",
            self.bwa_path,
            n_cores,
            self.bwa_realign_cutoff,
            self.bwa_realign_cutoff,
            self.bwa_seed_freq,
            sf_head,
            sf_ref,
            sf_reads,
            gv::SAMTOOLS_PATH,
            sf_out_sam
        );
        println!("Run command: {}", cmd);
        utils::run_cmd_to_file(&cmd, &format!("{}.std_out", sf_out_sam));
    }

    /// `align_exon_to_contig`
    pub fn align_exon_to_contig(&self, sf_ref: &str, sf_reads: &str, sf_out_bam: &str) {
        let cmd = format!(
            "{} mem -t {} {} {} | samtools view -hSb - | {} sort -o {} - && samtools index {}",
            self.bwa_path, self.n_jobs, sf_ref, sf_reads, gv::SAMTOOLS_PATH, sf_out_bam, sf_out_bam
        );
        println!("Run command: {}", cmd);
        utils::run_cmd_small_output(&cmd);
    }

    /// `realign_disc_reads`
    pub fn realign_disc_reads(&self, sf_ref: &str, sf_reads: &str, sf_out_sam: &str) {
        let cmd = format!("{} mem -t {} -o {} {} {}", self.bwa_path, self.n_jobs, sf_out_sam, sf_ref, sf_reads);
        utils::run_cmd_small_output(&cmd);
    }

    /// `realign_disc_reads_low_mem`
    pub fn realign_disc_reads_low_mem(&self, sf_ref: &str, sf_reads: &str, sf_out_sam: &str) {
        let cmd = format!("{} mem -t {} -c {} -o {} {} {}", self.bwa_path, self.n_jobs, self.bwa_seed_freq, sf_out_sam, sf_ref, sf_reads);
        utils::run_cmd_small_output(&cmd);
    }

    /// `realign_reads_to_bam`
    pub fn realign_reads_to_bam(&self, samtools: &str, sf_ref: &str, sf_reads: &str, sf_out_bam: &str) {
        let cmd = format!("{} mem -t {} {} {} | {} view -hSb - | {} sort - ", self.bwa_path, self.n_jobs, sf_ref, sf_reads, samtools, samtools);
        println!("Run command: {}", cmd);
        utils::run_cmd_to_file(&cmd, sf_out_bam);
        let cmd2 = format!("{} index {}", samtools, sf_out_bam);
        utils::run_cmd_small_output(&cmd2);
    }

    /// `is_full_map`
    pub fn is_full_map(cigar: &[(u32, u32)], max_clip_len: i64) -> bool {
        if cigar.is_empty() {
            return false;
        }
        let first = cigar[0];
        let last = cigar[cigar.len() - 1];
        if cigar.len() == 1 && first.0 == 0 {
            true
        } else if (first.0 == 4 || first.0 == 5) && (first.1 as i64) <= max_clip_len && last.0 == 0 {
            true
        } else {
            (last.0 == 4 || last.0 == 5) && (last.1 as i64) <= max_clip_len && first.0 == 0
        }
    }

    /// `get_fully_mapped_algnmts`
    pub fn get_fully_mapped_algnmts(&self, sf_sam: &str, sf_reference: &str, max_clip_len: i64, sf_fully_sam: &str, sf_unmap_fa: &str, sf_polya_fa: &str) -> Result<()> {
        let mut bamfile = AlignmentFile::open(sf_sam, sf_reference)?;
        let header = bamfile.header_view();
        let mut out = AlignmentWriter::create(sf_fully_sam, &header, true)?;
        let mut fout_fa = utils::create(sf_unmap_fa)?;
        let mut fout_polya = utils::create(sf_polya_fa)?;
        let mut err: Option<anyhow::Error> = None;
        bamfile.fetch_all(|a| {
            if err.is_some() {
                return;
            }
            if a.is_dup() || a.is_supp() || a.is_sec() {
                return;
            }
            let cig = a.cigar_tuples();
            if Self::is_full_map(&cig, max_clip_len) {
                if let Err(e) = out.write(a) {
                    err = Some(e);
                }
            } else {
                let seq = a.query_sequence().unwrap_or_default();
                let len_seq = seq.len() as i64;
                let r = if len_seq >= gv::MINIMUM_POLYA_CLIP && len_seq < gv::BWA_REALIGN_CUTOFF {
                    write!(fout_polya, ">{}\n{}\n", a.query_name(), seq)
                } else {
                    write!(fout_fa, ">{}\n{}\n", a.query_name(), seq)
                };
                if let Err(e) = r {
                    err = Some(e.into());
                }
            }
        })?;
        if let Some(e) = err {
            return Err(e);
        }
        drop(out);
        fout_fa.flush()?;
        fout_polya.flush()?;
        Ok(())
    }

    /// `append_to_sam`
    pub fn append_to_sam(&self, sf_in_sam: &str, sf_out: &str) -> Result<()> {
        let mut fout = utils::append(sf_out)?;
        for line in utils::read_raw_lines(sf_in_sam)? {
            if !line.starts_with('@') {
                fout.write_all(line.as_bytes())?;
            }
        }
        fout.flush()?;
        Ok(())
    }

    /// `_gnrt_SQ_from_fa`
    pub fn gnrt_sq_from_fa(&self, sf_fa: &str, sf_sq: &str) -> Result<()> {
        let mut fout = utils::create(sf_sq)?;
        for e in read_fastx(sf_fa)? {
            write!(fout, "@SQ\tSN:{}\tLN:{}\n", e.name, e.sequence.len())?;
        }
        fout.flush()?;
        Ok(())
    }

    /// `two_stage_realign`
    pub fn two_stage_realign(&self, sf_ref_cns: &str, sf_ref2: &str, sf_reads: &str, sf_out: &str) -> Result<()> {
        let sf_sam_cns = format!("{}_cns.sam", sf_out);
        let mut n_cores = self.n_jobs - 1;
        if n_cores <= 0 {
            n_cores = 1;
        }
        let cmd = format!(
            "{} mem -t {} -T {} -k {} -o {} {} {}",
            self.bwa_path, n_cores, self.bwa_realign_cutoff, self.bwa_realign_cutoff, sf_sam_cns, sf_ref_cns, sf_reads
        );
        utils::run_cmd_small_output(&cmd);

        let sf_fully_sam = format!("{}.sam", sf_out);
        let sf_unmap_fa = format!("{}_unmapped_phase1.fa", sf_out);
        let max_clip_len = 2;
        let sf_polya_fa = format!("{}_partial_polyA.fa", sf_out);
        self.get_fully_mapped_algnmts(&sf_sam_cns, sf_ref_cns, max_clip_len, &sf_fully_sam, &sf_unmap_fa, &sf_polya_fa)?;

        let sf_polya_sam = format!("{}_partial_polyA.sam", sf_out);
        self.realign_clipped_polya(sf_ref_cns, &sf_polya_fa, &sf_polya_sam);

        let sf_sq = format!("{}_sq.txt", sf_out);
        self.gnrt_sq_from_fa(sf_ref2, &sf_sq)?;
        {
            let mut fout = utils::append(&sf_sq)?;
            let sam = AlignmentFile::open(&sf_sam_cns, "")?;
            for (sn, ln) in sam.sq() {
                write!(fout, "@SQ\tSN:{}\tLN:{}\n", sn, ln)?;
            }
            fout.flush()?;
        }
        self.realign_clipped_reads_low_mem_add_head(sf_ref2, &sf_unmap_fa, &sf_sq, sf_out);
        self.append_to_sam(&sf_fully_sam, sf_out)?;
        self.append_to_sam(&sf_polya_sam, sf_out)?;

        for f in [&sf_sam_cns, &sf_fully_sam, &sf_unmap_fa, &sf_polya_fa, &sf_polya_sam] {
            utils::remove_file(f);
        }
        Ok(())
    }

    /// `calmd_from_algnmt`
    pub fn calmd_from_algnmt(&self, sf_sam: &str, sf_ref: &str, sf_calmd_sam: &str) {
        let cmd = format!("{} calmd -e {} {}", gv::SAMTOOLS_PATH, sf_sam, sf_ref);
        utils::run_cmd_to_file(&cmd, sf_calmd_sam);
    }
}
