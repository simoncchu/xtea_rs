use xtea::intervaltree::IntervalTree;

#[test]
fn intervaltree_matches_python_vendored_library() {
    let text = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/itree_py38.txt")).unwrap();
    let mut t = IntervalTree::new();
    let mut nq = 0;
    for line in text.lines() {
        let f: Vec<&str> = line.split(' ').collect();
        match f[0] {
            "TREE" => t = IntervalTree::new(),
            "ADD" => t.addi(f[1].parse().unwrap(), f[2].parse().unwrap()),
            "Q" => {
                nq += 1;
                let p: i64 = f[1].parse().unwrap();
                let exp = if f.len() > 2 { f[2] } else { "" };
                let got = t.at(p).iter().map(|(b, e)| format!("{}:{}", b, e)).collect::<Vec<_>>().join(",");
                assert_eq!(got, exp, "query {}", line);
            }
            _ => panic!(),
        }
    }
    assert!(nq > 1000);
}
