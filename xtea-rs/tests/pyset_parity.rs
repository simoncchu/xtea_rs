use xtea::py::set::PySet;

fn parse_i(s: &str) -> i64 { s.parse().unwrap() }
fn parse_t(s: &str) -> (i64, i64) { let (a, b) = s.split_once(':').unwrap(); (a.parse().unwrap(), b.parse().unwrap()) }

#[test]
fn pyset_matches_cpython_operation_log() {
    let text = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/set_ops_py38.txt")).unwrap();
    let mut tuple_mode = false;
    let mut ia: [PySet<i64>; 2] = [PySet::new(), PySet::new()];
    let mut ta: [PySet<(i64, i64)>; 2] = [PySet::new(), PySet::new()];
    let idx = |n: &str| if n == "a" { 0 } else { 1 };
    let mut checks = 0;
    for line in text.lines() {
        let f: Vec<&str> = line.split(' ').collect();
        match f[0] {
            "NEW" => { tuple_mode = f[1] == "T"; ia = [PySet::new(), PySet::new()]; ta = [PySet::new(), PySet::new()]; }
            "ADD" => { if tuple_mode { ta[idx(f[1])].add(parse_t(f[2])) } else { ia[idx(f[1])].add(parse_i(f[2])) } }
            "DISCARD" => { if tuple_mode { ta[idx(f[1])].discard(&parse_t(f[2])); } else { ia[idx(f[1])].discard(&parse_i(f[2])); } }
            "CLEAR" => { if tuple_mode { ta[idx(f[1])].clear() } else { ia[idx(f[1])].clear() } }
            "COPY" => { let (d, s) = (idx(f[1]), idx(f[2])); if tuple_mode { ta[d] = PySet::copy_of(&ta[s]) } else { ia[d] = PySet::copy_of(&ia[s]) } }
            "MERGE" => { let (d, s) = (idx(f[1]), idx(f[2])); if tuple_mode { let o = ta[s].clone(); ta[d].merge(&o) } else { let o = ia[s].clone(); ia[d].merge(&o) } }
            "DIFF" => { let (d, s) = (idx(f[1]), idx(f[2])); if tuple_mode { let o = ta[s].clone(); ta[d].difference_update(&o) } else { let o = ia[s].clone(); ia[d].difference_update(&o) } }
            "CHECK" => {
                checks += 1;
                let exp = if f.len() > 2 { f[2] } else { "" };
                let got = if tuple_mode {
                    ta[idx(f[1])].iter().map(|(a, b)| format!("{}:{}", a, b)).collect::<Vec<_>>().join(",")
                } else {
                    ia[idx(f[1])].iter().map(|a| a.to_string()).collect::<Vec<_>>().join(",")
                };
                assert_eq!(got, exp, "mismatch at line: {}", line);
            }
            _ => panic!("bad op {}", line),
        }
    }
    assert!(checks > 1000);
}
