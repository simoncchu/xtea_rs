//! `xtea_germline_cloud`: Rust port of `gnrt_pipeline_germline_cloud_v38.py`, the germline
//! wrapper of the cloud workflow (hg38).
fn main() {
    std::process::exit(xtea::pipeline::germline_cloud::main());
}
