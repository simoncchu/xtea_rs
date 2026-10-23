//! Port of `XReference` from `x_reference.py` (the `XChromosome` part lives in `chrom.rs`).

use indexmap::IndexMap;

/// `XReference.break_ref_to_bins`: chrom -> [(bin_start, bin_end)].
pub fn break_ref_to_bins(m_chrm_length: &IndexMap<String, i64>, bin_size: i64) -> IndexMap<String, Vec<(i64, i64)>> {
    let mut m_bins: IndexMap<String, Vec<(i64, i64)>> = IndexMap::new();
    for (chrm, &chrm_lenth) in m_chrm_length.iter() {
        let v = m_bins.entry(chrm.clone()).or_default();
        let mut tmp_pos = 0i64;
        while tmp_pos < chrm_lenth {
            let mut block_end = tmp_pos + bin_size;
            if block_end > chrm_lenth {
                block_end = chrm_lenth;
            }
            v.push((tmp_pos, block_end));
            tmp_pos = block_end;
        }
    }
    m_bins
}
