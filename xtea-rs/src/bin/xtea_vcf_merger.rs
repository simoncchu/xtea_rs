//! `xtea_vcf_merger`: Rust port of `x_vcf_merger.py`.
use xtea::tools::vcf_merger::{run, Opts};

fn main() {
    let opts: Opts = xtea::tools::parse_optparse_args();
    if let Err(e) = run(opts) {
        eprintln!("Error: {:#}", e);
        std::process::exit(1);
    }
}
