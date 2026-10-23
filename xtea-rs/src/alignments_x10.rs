//! 10X (linked-read) helpers of `x_alignments.py` (XBamInfo).

use crate::alignments::change_barcode_format;
use crate::bam::{AlignmentFile, PyAln};
use anyhow::Result;
use indexmap::IndexSet;

/// `XBamInfo.parse_barcodes_for_one_site` (insertion-ordered set of barcodes)
pub fn parse_barcodes_for_one_site(bamfile: &mut AlignmentFile, chrm: &str, pos: i64, i_extend: i64) -> Result<IndexSet<String>> {
    let mut set = IndexSet::new();
    bamfile.fetch(chrm, pos - i_extend, pos + i_extend, |a| {
        if a.is_dup() || a.is_supp() || a.is_unmap() || a.is_sec() {
            return;
        }
        if let Some(bx) = a.get_tag_string("BX") {
            set.insert(change_barcode_format(&bx));
        }
    })?;
    Ok(set)
}
