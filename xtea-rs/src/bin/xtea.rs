//! `xtea`: Rust port of `gnrt_pipeline_local_v38.py` (`bin/xtea`), pipeline generator for hg38.
fn main() {
    std::process::exit(xtea::pipeline::local_v38::main());
}
