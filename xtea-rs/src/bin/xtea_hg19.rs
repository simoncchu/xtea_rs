//! `xtea_hg19`: Rust port of `gnrt_pipeline_local.py` (`bin/xtea_hg19`), pipeline generator for hg19.
fn main() {
    std::process::exit(xtea::pipeline::local_hg19::main());
}
