//! `xtea_cmp`: Rust port of `x_cmp.py`.
use xtea::tools::cmp::{run, Opts};

fn main() {
    let opts: Opts = xtea::py::optparse::parse_opts("xtea_cmp");
    if let Err(e) = run(opts) {
        eprintln!("Error: {:#}", e);
        std::process::exit(1);
    }
}
