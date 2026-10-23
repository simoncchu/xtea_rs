//! A thin pysam-like layer over rust-htslib.
//!
//! The Python implementation accesses alignments through `pysam.AlignmentFile` and
//! `pysam.AlignedSegment`; the helpers here reproduce the semantics of the attributes xTea
//! uses (e.g. `cigar` op codes, `reference_end`, `query_sequence` being `None` for '*', the
//! argument checks of `fetch`).

use anyhow::{anyhow, Result};
use rust_htslib::bam::{self, record::Aux, Read};
use std::collections::HashMap;
use std::path::Path;

pub use rust_htslib::bam::Record;

/// pysam CIGAR operation codes.
pub const BAM_CMATCH: u32 = 0;
pub const BAM_CINS: u32 = 1;
pub const BAM_CDEL: u32 = 2;
pub const BAM_CREF_SKIP: u32 = 3;
pub const BAM_CSOFT_CLIP: u32 = 4;
pub const BAM_CHARD_CLIP: u32 = 5;
pub const BAM_CPAD: u32 = 6;
pub const BAM_CEQUAL: u32 = 7;
pub const BAM_CDIFF: u32 = 8;

/// pysam's MAX_POS
pub const MAX_POS: i64 = (1i64 << 31) - 1;

enum Inner {
    Indexed(bam::IndexedReader),
    Sequential(bam::Reader),
}

/// Equivalent of `pysam.AlignmentFile` opened for reading.
pub struct AlignmentFile {
    inner: Inner,
    names: Vec<String>,
    lengths: Vec<i64>,
    tids: HashMap<String, i32>,
    pub path: String,
}

fn index_exists(path: &str) -> bool {
    let cands = [
        format!("{}.bai", path),
        format!("{}.csi", path),
        format!("{}.crai", path),
        path.trim_end_matches(".bam").to_string() + ".bai",
        path.trim_end_matches(".cram").to_string() + ".crai",
    ];
    cands.iter().any(|c| Path::new(c).exists())
}

impl AlignmentFile {
    /// Open a SAM/BAM/CRAM file.  An index is used when available (needed by `fetch` with a
    /// region); `reference` is used for CRAM decoding (pysam's `reference_filename`).
    pub fn open(path: &str, reference: &str) -> Result<Self> {
        if !Path::new(path).exists() {
            return Err(anyhow!("file `{}` not found", path));
        }
        let lower = path.to_ascii_lowercase();
        let try_index = !(lower.ends_with(".sam")) && index_exists(path);
        let inner = if try_index {
            match bam::IndexedReader::from_path(path) {
                Ok(mut r) => {
                    if !reference.is_empty() && Path::new(reference).exists() {
                        let _ = r.set_reference(reference);
                    }
                    Inner::Indexed(r)
                }
                Err(_) => {
                    let mut r = bam::Reader::from_path(path).map_err(|e| anyhow!("cannot open {}: {}", path, e))?;
                    if !reference.is_empty() && Path::new(reference).exists() {
                        let _ = r.set_reference(reference);
                    }
                    Inner::Sequential(r)
                }
            }
        } else {
            let mut r = bam::Reader::from_path(path).map_err(|e| anyhow!("cannot open {}: {}", path, e))?;
            if !reference.is_empty() && Path::new(reference).exists() {
                let _ = r.set_reference(reference);
            }
            Inner::Sequential(r)
        };
        let header = match &inner {
            Inner::Indexed(r) => r.header().clone(),
            Inner::Sequential(r) => r.header().clone(),
        };
        let mut names = Vec::new();
        let mut lengths = Vec::new();
        let mut tids = HashMap::new();
        for (i, n) in header.target_names().iter().enumerate() {
            let s = String::from_utf8_lossy(n).into_owned();
            lengths.push(header.target_len(i as u32).unwrap_or(0) as i64);
            tids.insert(s.clone(), i as i32);
            names.push(s);
        }
        Ok(AlignmentFile { inner, names, lengths, tids, path: path.to_string() })
    }

    /// `samfile.references`
    pub fn references(&self) -> &[String] {
        &self.names
    }

    /// `samfile.lengths`
    pub fn lengths(&self) -> &[i64] {
        &self.lengths
    }

    /// `samfile.get_tid(name)` (-1 when absent).
    pub fn get_tid(&self, name: &str) -> i32 {
        *self.tids.get(name).unwrap_or(&-1)
    }

    /// Name for a tid (None for -1 / out of range), like `reference_name`.
    pub fn reference_name(&self, tid: i32) -> Option<&str> {
        if tid < 0 {
            return None;
        }
        self.names.get(tid as usize).map(|s| s.as_str())
    }

    /// Map of chromosome id to name with -1 -> "*" (xTea's `_get_chrm_id_name`).
    pub fn chrm_id_name(&self) -> HashMap<i32, String> {
        let mut m = HashMap::new();
        for (i, n) in self.names.iter().enumerate() {
            m.insert(i as i32, n.clone());
        }
        m.insert(-1, "*".to_string());
        m
    }

    pub fn header_view(&self) -> bam::HeaderView {
        match &self.inner {
            Inner::Indexed(r) => r.header().clone(),
            Inner::Sequential(r) => r.header().clone(),
        }
    }

    /// `samfile.header.to_dict()['SQ']` style access: (name, length) pairs in header order.
    pub fn sq(&self) -> Vec<(String, i64)> {
        self.names.iter().cloned().zip(self.lengths.iter().cloned()).collect()
    }

    /// `samfile.header.to_dict()['SQ']` with all tags: one (tag, value) list per `@SQ` line.
    /// Like pysam, falls back to the binary reference list when the header text has no `@SQ`
    /// lines; `None` when there are no references at all (`'SQ' not in header`).
    pub fn header_sq_records(&self) -> Option<Vec<Vec<(String, String)>>> {
        let hv = self.header_view();
        let text = String::from_utf8_lossy(hv.as_bytes()).into_owned();
        let mut out = Vec::new();
        for line in text.split('\n') {
            let line = line.trim_matches(|c| c == ' ' || c == '\0');
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.first() != Some(&"@SQ") {
                continue;
            }
            let mut rec: Vec<(String, String)> = Vec::new();
            for f in &fields[1..] {
                if let Some((k, v)) = f.split_once(':') {
                    let v = if k == "LN" {
                        v.trim().parse::<i64>().map(|x| x.to_string()).unwrap_or_else(|_| v.to_string())
                    } else {
                        v.to_string()
                    };
                    match rec.iter().position(|(ek, _)| ek == k) {
                        Some(p) => rec[p].1 = v,
                        None => rec.push((k.to_string(), v)),
                    }
                }
            }
            out.push(rec);
        }
        if out.is_empty() {
            for (n, l) in self.names.iter().zip(self.lengths.iter()) {
                out.push(vec![("LN".to_string(), l.to_string()), ("SN".to_string(), n.clone())]);
            }
        }
        if out.is_empty() {
            None
        } else {
            Some(out)
        }
    }

    /// `fetch(contig, start, stop)` calling `f` for each record (pysam semantics).
    pub fn fetch<F: FnMut(&Record)>(&mut self, contig: &str, start: i64, stop: i64, mut f: F) -> Result<()> {
        if start > stop {
            return Err(anyhow!("invalid coordinates: start ({}) > stop ({})", start, stop));
        }
        if start < 0 || start >= MAX_POS {
            return Err(anyhow!("start out of range ({})", start));
        }
        let stop = stop.min(MAX_POS);
        let tid = self.get_tid(contig);
        if tid < 0 {
            return Err(anyhow!("invalid contig `{}`", contig));
        }
        match &mut self.inner {
            Inner::Indexed(r) => {
                r.fetch((tid, start, stop)).map_err(|e| anyhow!("fetch failed: {}", e))?;
                let mut rec = Record::new();
                while let Some(res) = r.read(&mut rec) {
                    res.map_err(|e| anyhow!("error reading {}: {}", self.path, e))?;
                    f(&rec);
                }
                Ok(())
            }
            Inner::Sequential(_) => Err(anyhow!("fetch called on bamfile without index")),
        }
    }

    /// `fetch(contig)` - whole contig.
    pub fn fetch_contig<F: FnMut(&Record)>(&mut self, contig: &str, f: F) -> Result<()> {
        self.fetch(contig, 0, MAX_POS, f)
    }

    /// `fetch()` without region on an un-indexed file / `fetch(until_eof=True)`: all records in
    /// file order.
    pub fn fetch_all<F: FnMut(&Record)>(&mut self, mut f: F) -> Result<()> {
        match &mut self.inner {
            Inner::Sequential(r) => {
                let mut rec = Record::new();
                while let Some(res) = r.read(&mut rec) {
                    res.map_err(|e| anyhow!("error reading {}: {}", self.path, e))?;
                    f(&rec);
                }
                Ok(())
            }
            Inner::Indexed(_) => {
                // re-open sequentially to stream the whole file in order
                let mut r = bam::Reader::from_path(&self.path).map_err(|e| anyhow!("cannot open {}: {}", self.path, e))?;
                let mut rec = Record::new();
                while let Some(res) = r.read(&mut rec) {
                    res.map_err(|e| anyhow!("error reading {}: {}", self.path, e))?;
                    f(&rec);
                }
                Ok(())
            }
        }
    }

    /// Collect the records of a region into owned values.
    pub fn fetch_vec(&mut self, contig: &str, start: i64, stop: i64) -> Result<Vec<Record>> {
        let mut v = Vec::new();
        self.fetch(contig, start, stop, |r| v.push(r.clone()))?;
        Ok(v)
    }
}

/// Helpers reproducing pysam `AlignedSegment` attributes.
pub trait PyAln {
    fn is_dup(&self) -> bool;
    fn is_supp(&self) -> bool;
    fn is_unmap(&self) -> bool;
    fn mate_is_unmapped(&self) -> bool;
    fn is_sec(&self) -> bool;
    fn is_read1(&self) -> bool;
    fn is_read2(&self) -> bool;
    fn is_rev(&self) -> bool;
    fn mate_is_rev(&self) -> bool;
    fn reference_start(&self) -> i64;
    fn reference_id(&self) -> i32;
    fn next_reference_id(&self) -> i32;
    fn next_reference_start(&self) -> i64;
    fn mapping_quality(&self) -> i64;
    fn template_length(&self) -> i64;
    fn query_name(&self) -> String;
    fn cigar_tuples(&self) -> Vec<(u32, u32)>;
    fn query_sequence(&self) -> Option<String>;
    fn query_qualities(&self) -> Option<Vec<u8>>;
    fn reference_end(&self) -> Option<i64>;
    fn query_length(&self) -> i64;
    fn has_tag(&self, tag: &str) -> bool;
    fn get_tag_string(&self, tag: &str) -> Option<String>;
}

impl PyAln for Record {
    #[inline]
    fn is_dup(&self) -> bool {
        self.flags() & 0x400 != 0
    }
    #[inline]
    fn is_supp(&self) -> bool {
        self.flags() & 0x800 != 0
    }
    #[inline]
    fn is_unmap(&self) -> bool {
        self.flags() & 0x4 != 0
    }
    #[inline]
    fn mate_is_unmapped(&self) -> bool {
        self.flags() & 0x8 != 0
    }
    #[inline]
    fn is_sec(&self) -> bool {
        self.flags() & 0x100 != 0
    }
    #[inline]
    fn is_read1(&self) -> bool {
        self.flags() & 0x40 != 0
    }
    #[inline]
    fn is_read2(&self) -> bool {
        self.flags() & 0x80 != 0
    }
    #[inline]
    fn is_rev(&self) -> bool {
        self.flags() & 0x10 != 0
    }
    #[inline]
    fn mate_is_rev(&self) -> bool {
        self.flags() & 0x20 != 0
    }
    #[inline]
    fn reference_start(&self) -> i64 {
        self.pos()
    }
    #[inline]
    fn reference_id(&self) -> i32 {
        self.tid()
    }
    #[inline]
    fn next_reference_id(&self) -> i32 {
        self.mtid()
    }
    #[inline]
    fn next_reference_start(&self) -> i64 {
        self.mpos()
    }
    #[inline]
    fn mapping_quality(&self) -> i64 {
        self.mapq() as i64
    }
    #[inline]
    fn template_length(&self) -> i64 {
        self.insert_size()
    }
    fn query_name(&self) -> String {
        String::from_utf8_lossy(self.qname()).into_owned()
    }
    fn cigar_tuples(&self) -> Vec<(u32, u32)> {
        self.raw_cigar().iter().map(|&c| (c & 0xf, c >> 4)).collect()
    }
    fn query_sequence(&self) -> Option<String> {
        if self.seq_len() == 0 {
            return None;
        }
        Some(String::from_utf8_lossy(&self.seq().as_bytes()).into_owned())
    }
    fn query_qualities(&self) -> Option<Vec<u8>> {
        let q = self.qual();
        if q.is_empty() || q[0] == 0xff {
            return None;
        }
        Some(q.to_vec())
    }
    fn reference_end(&self) -> Option<i64> {
        if self.is_unmap() || self.raw_cigar().is_empty() {
            return None;
        }
        let mut end = self.pos();
        for &c in self.raw_cigar() {
            let op = c & 0xf;
            let len = (c >> 4) as i64;
            match op {
                BAM_CMATCH | BAM_CDEL | BAM_CREF_SKIP | BAM_CEQUAL | BAM_CDIFF => end += len,
                _ => {}
            }
        }
        Some(end)
    }
    #[inline]
    fn query_length(&self) -> i64 {
        self.seq_len() as i64
    }
    fn has_tag(&self, tag: &str) -> bool {
        self.aux(tag.as_bytes()).is_ok()
    }
    fn get_tag_string(&self, tag: &str) -> Option<String> {
        match self.aux(tag.as_bytes()) {
            Ok(a) => Some(aux_to_py_string(&a)),
            Err(_) => None,
        }
    }
}

/// `str(get_tag(...))` for the scalar tag types.
pub fn aux_to_py_string(a: &Aux) -> String {
    match a {
        Aux::Char(c) => (*c as char).to_string(),
        Aux::I8(v) => v.to_string(),
        Aux::U8(v) => v.to_string(),
        Aux::I16(v) => v.to_string(),
        Aux::U16(v) => v.to_string(),
        Aux::I32(v) => v.to_string(),
        Aux::U32(v) => v.to_string(),
        Aux::Float(v) => crate::py::fmt::float_repr(*v as f64),
        Aux::Double(v) => crate::py::fmt::float_repr(*v),
        Aux::String(s) => s.to_string(),
        Aux::HexByteArray(s) => s.to_string(),
        _ => String::new(),
    }
}

/// `pysam.AlignmentFile(sf, 'w'|'wb', template=...)` writer.
pub struct AlignmentWriter {
    w: bam::Writer,
}

impl AlignmentWriter {
    pub fn create(path: &str, header: &bam::HeaderView, sam: bool) -> Result<Self> {
        let h = bam::Header::from_template(header);
        let fmt = if sam { bam::Format::Sam } else { bam::Format::Bam };
        let w = bam::Writer::from_path(path, &h, fmt).map_err(|e| anyhow!("cannot create {}: {}", path, e))?;
        Ok(AlignmentWriter { w })
    }
    pub fn write(&mut self, r: &Record) -> Result<()> {
        self.w.write(r).map_err(|e| anyhow!("write failed: {}", e))
    }
}
