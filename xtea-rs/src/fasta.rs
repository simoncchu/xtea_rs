//! pysam.FastaFile / pysam.FastxFile equivalents.

use anyhow::{anyhow, Result};
use rust_htslib::htslib;
use std::ffi::{CStr, CString};
use std::io::{BufRead, BufReader};

/// Indexed FASTA reader (htslib faidx), mirroring `pysam.FastaFile`.
pub struct FastaFile {
    inner: *mut htslib::faidx_t,
    names: Vec<String>,
}

unsafe impl Send for FastaFile {}

impl FastaFile {
    pub fn open(path: &str) -> Result<Self> {
        let cpath = CString::new(path)?;
        let inner = unsafe { htslib::fai_load(cpath.as_ptr()) };
        if inner.is_null() {
            return Err(anyhow!("could not open FASTA file {}", path));
        }
        let n = unsafe { htslib::faidx_nseq(inner) };
        let mut names = Vec::with_capacity(n.max(0) as usize);
        for i in 0..n {
            let p = unsafe { htslib::faidx_iseq(inner, i) };
            let s = unsafe { CStr::from_ptr(p) }.to_string_lossy().into_owned();
            names.push(s);
        }
        Ok(FastaFile { inner, names })
    }

    /// `FastaFile.references`
    pub fn references(&self) -> &[String] {
        &self.names
    }

    pub fn has_reference(&self, name: &str) -> bool {
        self.names.iter().any(|n| n == name)
    }

    /// Length of a sequence (`get_reference_length`), or None when absent.
    pub fn length(&self, name: &str) -> Option<i64> {
        let c = CString::new(name).ok()?;
        let has = unsafe { htslib::faidx_has_seq(self.inner, c.as_ptr()) };
        if has == 0 {
            return None;
        }
        Some(unsafe { htslib::faidx_seq_len64(self.inner, c.as_ptr()) })
    }

    /// `FastaFile.fetch(reference, start, end)` (0-based, half-open).  Follows pysam's argument
    /// checks: a negative start or start > end is an error, like pysam's ValueError.
    pub fn fetch(&self, name: &str, start: i64, end: i64) -> Result<String> {
        if start > end {
            return Err(anyhow!("invalid coordinates: start ({}) > stop ({})", start, end));
        }
        if start < 0 {
            return Err(anyhow!("start out of range ({})", start));
        }
        if start == end {
            return Ok(String::new());
        }
        let c = CString::new(name)?;
        let mut len: htslib::hts_pos_t = 0;
        let ptr = unsafe { htslib::faidx_fetch_seq64(self.inner, c.as_ptr(), start, end - 1, &mut len) };
        if ptr.is_null() {
            return Err(anyhow!("failure when retrieving sequence on '{}'", name));
        }
        let s = unsafe {
            let bytes = std::slice::from_raw_parts(ptr as *const u8, len.max(0) as usize);
            let s = String::from_utf8_lossy(bytes).into_owned();
            libc::free(ptr as *mut libc::c_void);
            s
        };
        Ok(s)
    }

    /// Like `fetch` but silently returns an empty string on errors.
    pub fn fetch_or_empty(&self, name: &str, start: i64, end: i64) -> String {
        self.fetch(name, start, end).unwrap_or_default()
    }
}

impl Drop for FastaFile {
    fn drop(&mut self) {
        unsafe { htslib::fai_destroy(self.inner) };
    }
}

/// One record of `pysam.FastxFile`.
#[derive(Clone, Debug, Default)]
pub struct FastxRecord {
    pub name: String,
    pub comment: Option<String>,
    pub sequence: String,
    pub quality: Option<String>,
}

impl FastxRecord {
    /// `str(entry)` as produced by pysam (FASTA or FASTQ text without trailing newline).
    pub fn to_pysam_string(&self) -> String {
        let mut head = self.name.clone();
        if let Some(c) = &self.comment {
            head.push(' ');
            head.push_str(c);
        }
        match &self.quality {
            Some(q) => format!("@{}\n{}\n+\n{}", head, self.sequence, q),
            None => format!(">{}\n{}", head, self.sequence),
        }
    }
}

/// Read all records of a FASTA/FASTQ file (plain or gzip compressed) like `pysam.FastxFile`.
pub fn read_fastx(path: &str) -> Result<Vec<FastxRecord>> {
    let f = std::fs::File::open(path).map_err(|e| anyhow!("cannot open {}: {}", path, e))?;
    let mut magic = [0u8; 2];
    let reader: Box<dyn BufRead> = {
        use std::io::Read;
        let mut f2 = std::fs::File::open(path)?;
        let n = f2.read(&mut magic)?;
        if n == 2 && magic == [0x1f, 0x8b] {
            Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(f)))
        } else {
            Box::new(BufReader::new(f))
        }
    };
    let mut out = Vec::new();
    let mut lines = reader.lines();
    let mut pending: Option<String> = None;
    loop {
        let header = match pending.take() {
            Some(h) => h,
            None => match lines.next() {
                Some(l) => l?,
                None => break,
            },
        };
        let header = header.trim_end_matches('\r').to_string();
        if header.is_empty() {
            continue;
        }
        let first = header.as_bytes()[0];
        let rest = &header[1..];
        let (name, comment) = match rest.split_once(|c: char| c == ' ' || c == '\t') {
            Some((n, c)) => (n.to_string(), Some(c.to_string())),
            None => (rest.to_string(), None),
        };
        if first == b'>' {
            let mut seq = String::new();
            for l in lines.by_ref() {
                let l = l?;
                let l = l.trim_end_matches('\r');
                if l.starts_with('>') {
                    pending = Some(l.to_string());
                    break;
                }
                seq.push_str(l);
            }
            out.push(FastxRecord { name, comment, sequence: seq, quality: None });
        } else if first == b'@' {
            let mut seq = String::new();
            let mut found_plus = false;
            for l in lines.by_ref() {
                let l = l?;
                let l = l.trim_end_matches('\r');
                if l.starts_with('+') {
                    found_plus = true;
                    break;
                }
                seq.push_str(l);
            }
            let mut qual = String::new();
            if found_plus {
                while qual.len() < seq.len() {
                    match lines.next() {
                        Some(l) => {
                            let l = l?;
                            qual.push_str(l.trim_end_matches('\r'));
                        }
                        None => break,
                    }
                }
            }
            out.push(FastxRecord { name, comment, sequence: seq, quality: Some(qual) });
        } else {
            // skip unexpected lines
            continue;
        }
    }
    Ok(out)
}

/// Reverse complement used throughout xTea (`gnrt_reverse_complementary` with the standard
/// chr_map: A/a->T, T/t->A, C/c->G, G/g->C, N/n->N, others -> N).
pub fn reverse_complement_xtea(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| match b {
            b'A' | b'a' => 'T',
            b'T' | b't' => 'A',
            b'C' | b'c' => 'G',
            b'G' | b'g' => 'C',
            _ => 'N',
        })
        .collect()
}
