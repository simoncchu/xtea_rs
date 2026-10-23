//! `xtea_population_af`: Rust port of `x_population_AF.py` (`xtea_population_af <vcf> <out_csv>`;
//! like the script, the arguments are taken positionally from argv without option parsing).

fn main() {
    let args: Vec<std::ffi::OsString> = std::env::args_os().collect();
    if let Err(e) = xtea::tools::population_af::main_with_args(&args) {
        eprintln!("Error: {:#}", e);
        std::process::exit(1);
    }
}
