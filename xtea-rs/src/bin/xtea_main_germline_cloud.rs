//! `xtea_main_germline_cloud`: Rust port of `x_TEA_main_germline_cloud.py` (the entry point
//! used by the cloud germline workflow).

fn main() {
    let opts: xtea::cli::CloudOpts = xtea::py::optparse::parse_opts("xtea_main_germline_cloud");
    if let Err(e) = xtea::cli::run_cloud(opts) {
        eprintln!("Error: {:#}", e);
        std::process::exit(1);
    }
}
