//! Chromosome name helpers: `XChromosome.is_decoy_contig_chrms` (x_reference.py) and
//! `XIntemediateSites.is_decoy_contig_chrms` (x_intermediate_sites.py), which differ.

use crate::globals as gv;

/// `XChromosome().is_decoy_contig_chrms(chrm)`
pub fn is_decoy_contig_chrms(chrm: &str) -> bool {
    if chrm.split('_').count() > 1 {
        return true;
    } else if chrm == "hs37d5" {
        return true;
    }
    if chrm.len() > 3 {
        let p = |n: usize| -> &str { if chrm.len() >= n { &chrm[..n] } else { chrm } };
        if p(3) == "HLA" || p(3) == "HPV" || p(3) == "HIV" || p(3) == "CMV" || p(3) == "MCV" || p(2) == "SV" || p(4) == "KSHV" || p(5) == "decoy" || p(6) == "random" {
            return true;
        }
    }
    if chrm.split('-').count() > 1 {
        return true;
    }
    if chrm.contains("hs38") || chrm.contains("hs37") {
        return true;
    }
    if !gv::mit_switch_on() && (chrm == "MT" || chrm == "chrMT" || chrm == "chrM") {
        return true;
    }
    chrm.split('.').count() > 1
}

/// `XIntemediateSites().is_decoy_contig_chrms(chrm)` (a slightly different rule set).
pub fn is_decoy_contig_chrms_sites(chrm: &str) -> bool {
    if chrm.split('_').count() > 1 {
        return true;
    } else if chrm == "hs37d5" {
        return true;
    }
    if !gv::mit_switch_on() && (chrm == "MT" || chrm == "chrMT" || chrm == "chrM") {
        return true;
    }
    chrm.split('.').count() > 1
}

/// `TELocator._is_decoy_contig_chrms`
pub fn is_decoy_contig_chrms_locator(chrm: &str) -> bool {
    chrm.split('_').count() > 1 || chrm == "hs37d5"
}
