//! `xtea_main`: Rust port of `x_TEA_main.py`.

fn main() {
    // optparse semantics: unique-prefix long options, values taken verbatim, last repeat wins,
    // optparse-style help
    let opts: xtea::cli::Opts = xtea::py::optparse::parse_opts("xtea_main");
    if let Err(e) = xtea::cli::run(opts) {
        eprintln!("Error: {:#}", e);
        std::process::exit(1);
    }
}
