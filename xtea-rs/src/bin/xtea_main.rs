//! `xtea_main`: Rust port of `x_TEA_main.py`.
use clap::Parser;

fn main() {
    let opts = xtea::cli::Opts::parse();
    if let Err(e) = xtea::cli::run(opts) {
        eprintln!("Error: {:#}", e);
        std::process::exit(1);
    }
}
