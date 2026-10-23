//! End-to-end parity of the stand-alone tools (`xtea_vcf_merger`, `xtea_cmp`,
//! `xtea_population_af`) with the original Python scripts.  For every case of
//! `tests/data/tools/cases.tsv`, `<case>/input` is copied to a temporary directory, the tool
//! runs there with the listed arguments, and its stdout, exit status and every file it creates
//! or changes must equal what CPython 3.8 produced (`<case>/expected`, regenerated with
//! `tests/data/tools/gen_tools_fixtures.py`; large files are compared by FNV-1a hash).
//! The tools sort through the shell like the scripts, so GNU `sort` must be available.

use std::collections::BTreeMap;
use std::path::Path;
use std::process::Command;

fn walk(dir: &Path, base: &Path, out: &mut BTreeMap<String, Vec<u8>>) {
    for e in std::fs::read_dir(dir).unwrap() {
        let p = e.unwrap().path();
        if p.is_dir() {
            walk(&p, base, out);
        } else {
            out.insert(p.strip_prefix(base).unwrap().to_string_lossy().into_owned(), std::fs::read(&p).unwrap());
        }
    }
}

fn snapshot(dir: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut m = BTreeMap::new();
    if dir.exists() {
        walk(dir, dir, &mut m);
    }
    m
}

fn copy_tree(src: &Path, dst: &Path) {
    std::fs::create_dir_all(dst).unwrap();
    for e in std::fs::read_dir(src).unwrap() {
        let p = e.unwrap().path();
        let t = dst.join(p.file_name().unwrap());
        if p.is_dir() {
            copy_tree(&p, &t);
        } else {
            std::fs::copy(&p, &t).unwrap();
        }
    }
}

fn fnv1a64(data: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in data {
        h = (h ^ b as u64).wrapping_mul(0x100000001b3);
    }
    format!("{:016x}\n", h)
}

#[test]
fn tools_match_python_outputs() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/tools");
    let cases = std::fs::read_to_string(root.join("cases.tsv")).unwrap();
    let mut n_cases = 0;
    for line in cases.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let (name, tool, args) = (f[0], f[1], f[2]);
        let exe = match tool {
            "merger" => env!("CARGO_BIN_EXE_xtea_vcf_merger"),
            "cmp" => env!("CARGO_BIN_EXE_xtea_cmp"),
            "af" => env!("CARGO_BIN_EXE_xtea_population_af"),
            _ => panic!("unknown tool {}", tool),
        };
        let run = std::env::temp_dir().join(format!("xtea_tools_parity_{}_{}", std::process::id(), name));
        let _ = std::fs::remove_dir_all(&run);
        copy_tree(&root.join(name).join("input"), &run);
        let before = snapshot(&run);
        let out = Command::new(exe).args(args.split(' ')).current_dir(&run).output().unwrap();
        let after = snapshot(&run);

        let exp_dir = root.join(name).join("expected");
        assert_eq!(
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&std::fs::read(exp_dir.join("stdout.txt")).unwrap()),
            "{}: stdout",
            name
        );
        let exp_code: i32 = std::fs::read_to_string(exp_dir.join("exit_code.txt")).unwrap().trim().parse().unwrap();
        assert_eq!(out.status.success(), exp_code == 0, "{}: exit status", name);

        let changed: BTreeMap<String, Vec<u8>> = after.into_iter().filter(|(k, v)| before.get(k) != Some(v)).collect();
        let expected = snapshot(&exp_dir.join("files"));
        let mut exp_names: Vec<&str> = expected.keys().map(|k| k.strip_suffix(".fnv1a64").unwrap_or(k)).collect();
        exp_names.sort();
        let got_names: Vec<&str> = changed.keys().map(|k| k.as_str()).collect();
        assert_eq!(got_names, exp_names, "{}: created/changed files", name);
        for (k, v) in &expected {
            match k.strip_suffix(".fnv1a64") {
                Some(orig) => assert_eq!(fnv1a64(&changed[orig]).as_bytes(), v.as_slice(), "{}: {} differs", name, orig),
                None => assert!(changed[k] == *v, "{}: {} differs", name, k),
            }
        }
        std::fs::remove_dir_all(&run).unwrap();
        n_cases += 1;
    }
    assert_eq!(n_cases, 6);
}
