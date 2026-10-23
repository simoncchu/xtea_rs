#!/usr/bin/env python3
"""Regenerate the fixtures of tests/tools_parity.rs with the original Python 3.8 scripts.

usage: gen_tools_fixtures.py <python3.8 interpreter> <xtea python source dir> <sim1 data dir>

<sim1 data dir> holds SIM1/{L1,Alu,SVA}/candidate_disc_filtered_cns.txt and truth.tsv (outputs
of the xTea simulation validation run).  For every case, <case>/input is the working directory
the tool runs in and <case>/expected holds the files the Python run created or changed, its
stdout (stdout.txt) and exit status (exit_code.txt).  Files larger than 64 kB are stored as
their FNV-1a 64 hash (<file>.fnv1a64).  `sort` must be GNU coreutils sort.
"""
import os
import random
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PY, PYSRC, SIMDATA = sys.argv[1], sys.argv[2], sys.argv[3]
SCRIPTS = {"merger": "x_vcf_merger.py", "cmp": "x_cmp.py", "af": "x_population_AF.py"}
MAX_PLAIN = 64 * 1024


def write(path, text):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "wb") as f:
        f.write(text.encode() if isinstance(text, str) else text)


def fnv1a64(data):
    h = 0xCBF29CE484222325
    for b in data:
        h = ((h ^ b) * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return "%016x\n" % h


HEAD = ["##fileformat=VCFv4.2", "##reference=genome.fa", "##contig=<ID=chr1,length=6000000>",
        "##contig=<ID=chr2,length=4000000>", "##contig=<ID=chr10,length=3000000>", "##contig=<ID=chr1_random,length=5>",
        "##contig=<ID=chrUn_x,length=5>", "##contig=<ID=HLA-A,length=5>", "##source=xTea_contig",
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">']


def rec(chrm, pos, rep, svlen, end, sub, td, gt, ref_rep="not_in_copy", gene="not_gene_region"):
    info = "SVTYPE=INS:ME:%s;SVLEN=%d;END=%d;TSD=+AC;TSDLEN=2;SUBTYPE=%s;TD_SRC=%s;STRAND=+;AF=0.5;REF_REP=%s;GENE_INFO=%s" % (
        rep, svlen, end, sub, td, ref_rep, gene)
    return "\t".join([chrm, str(pos), ".", "A", "<INS:ME:%s>" % rep, ".", "PASS", info, "GT", gt])


def vcf(sample, records, sep="\n"):
    return sep.join(HEAD + ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample] + records) + sep


def build_merge_sample(d):
    write(os.path.join(d, "sample_info.txt"), "sA male S1 EUR\nsB female S2 EUR\nsC male S3 AFR\nsA male S1b EUR\n")
    write(os.path.join(d, "rslt_list.txt"), "sA res/sA\nsB res/sB/\nsC res/sC\nsD res/sD\n")
    a = "res/sA/"
    write(os.path.join(d, a + "L1/sA_LINE1.vcf"), vcf("sA", [
        rec("chr1", 1000, "LINE1", 900, 1015, "one_side", "not_transduction", "0/1"),
        rec("chr1", 5000, "LINE1", 900, 5015, "two_side_tprt_both", "not_transduction", "1/1"),
        rec("chr1", 9000, "LINE1", 900, 9015, "two_side", "chr2:5~sibling", "0/1"),
        rec("chr10", 777, "LINE1", 900, 790, "one_side", "not_transduction", "0/1"),
        rec("chr2", 50, "LINE1", 900, 60, "one_side", "not_transduction", "0/1"), "chrZ\t1"]))
    write(os.path.join(d, a + "SVA/sA_SVA.vcf"), vcf("sA", [
        rec("chr1", 1030, "SVA", 1500, 1040, "two_side_tprt_both", "not_transduction", "0/1"),  # replaces 1000
        rec("chr1", 5040, "SVA", 1500, 5050, "two_side_tprt_both", "not_transduction", "0/1"),  # both kept
        rec("chr1", 9049, "SVA", 1500, 9060, "one_side", "not_transduction", "0/1"),  # both kept
        rec("chr10", 777, "SVA", 1500, 790, "two_side_tprt_both", "not_transduction", "1/1"),  # same position
        rec("chr2", 100, "SVA", 1500, 110, "two_side_tprt_both", "not_transduction", "0/1")], sep="\r\n"))
    write(os.path.join(d, a + "Alu/sA_ALU.vcf"), vcf("sA", [
        rec("chr1", 1035, "ALU", 280, 1050, "two_side_tprt_both", "not_transduction", "0/1"),
        rec("chr1", 990, "ALU", 280, 1000, "two_side_tprt_both", "not_transduction", "0/1"),  # hits the blanked 1000
        rec("chr10", 700, "ALU", 280, 710, "one_side", "not_transduction", "0/1")]))
    write(os.path.join(d, a + "HERV/sA_HERV.vcf"), vcf("sA", [
        rec("chr10", 727, "HERV-K", 5000, 740, "one_side", "not_transduction", "0/1"),
        rec("chr2", 99, "HERV-K", 5000, 110, "one_side", "not_transduction", "0/1")]))
    write(os.path.join(d, "res/sB/SVA/sB_SVA.vcf"), vcf("sB", [rec("chr2", 300, "SVA", 1500, 310, "one_side", "x", "0/1")]))
    write(os.path.join(d, "res/sC/Alu/sC_ALU.vcf"), vcf("sC", [rec("chr1", 30, "ALU", 280, 45, "one_side", "x", "0/1")]))
    os.makedirs(os.path.join(d, "out"))


def build_merge_population(d):
    rng = random.Random(7)
    samples = [("NA1", "EUR"), ("NA2", "EUR"), ("hg3", "AFR"), ("HG4", "AFR"), ("HG5", "EAS"), ("S6", "AMR")]
    write(os.path.join(d, "sample_info.txt"), "".join("%s male B%d %s\n" % (s, i, p) for i, (s, p) in enumerate(samples)) +
          "HG9 female B9 SAS\n")
    write(os.path.join(d, "rslt_list.txt"), "".join("%s v/%s.vcf\n" % (s, s) for s, _ in samples) + "UNK v/UNK.vcf\n")
    base = [("chr1", rng.randint(1000, 5900000)) for _ in range(10)] + [("chr2", rng.randint(1000, 3900000)) for _ in range(6)]
    base += [("chr10", 300000 + 9 * i) for i in range(12)]  # one dense cluster (CPython set order matters)
    base += [("chr1", 100000), ("chr1", 100050), ("chr1", 100101), ("chr1", 100400)]
    tds = ["not_transduction"] * 4 + ["chr1:10-20~0~orphan", "chr2:7~sibling", ""]
    gts = ["0/1", "0/1", "1/1", "0|1", "1|1", "1/0", "./."]
    for s, p in samples + [("UNK", "EUR")]:
        recs = []
        for (c, pos) in base:
            if rng.random() < 0.45:
                continue
            q = pos + rng.choice([0, 0, 1, -3, 12, 35])
            recs.append((c, q, rec(c, q, rng.choice(["LINE1", "ALU", "SVA"]), rng.choice([281, 6019]), q + rng.choice([14, 15]),
                                   "two_side", rng.choice(tds), rng.choice(gts), rng.choice(["r1", "r2"]),
                                   rng.choice(["g1", "g2"]))))
        recs.sort()
        write(os.path.join(d, "v", s + ".vcf"), vcf(s, [r[2] for r in recs]))
    os.makedirs(os.path.join(d, "out"))


def truth_noheader():
    with open(os.path.join(SIMDATA, "truth.tsv")) as f:
        return "".join("\t".join(l.split("\t")[:3]) + "\n" for l in f if not l.startswith("#"))


def build_cmp_single(d):
    shutil.copyfile(os.path.join(SIMDATA, "SIM1/L1/candidate_disc_filtered_cns.txt"), os.path.join(d, "calls.txt"))
    write(os.path.join(d, "bcmk.txt"), truth_noheader())
    write(os.path.join(d, "exclude.txt"), "#excluded\nchr1\t149000\nchr2\t1572600\n")
    os.makedirs(os.path.join(d, "out"))


def build_cmp_group(d):
    lines = []
    for t, rep in (("L1", "LINE1"), ("Alu", "ALU"), ("SVA", "SVA")):
        with open(os.path.join(SIMDATA, "SIM1", t, "candidate_disc_filtered_cns.txt")) as f:
            write(os.path.join(d, t + ".txt"), "".join("\t".join(l.split("\t")[:6]) + "\n" for l in f))
        lines.append("%s.txt %s xTea_30x" % (t, rep))
    lines.append("L1.txt LINE1 xTea_10x")
    write(os.path.join(d, "list.txt"), "\n".join(lines) + "\n")
    write(os.path.join(d, "bcmk.txt"), truth_noheader())
    os.makedirs(os.path.join(d, "out"))


def build_cmp_unicode(d):
    write(os.path.join(d, "calls.txt"), "chr1 149129 a\nchr1\x1c٤٥٧٢٧٥\x1cb\n"
          "chr1　946109　c\nchr1\t946109\tdup\nchr\t10\tx\nchrchr1\t20\ty\nchr2\t1572575\nchr2\t3583035\tlast")
    write(os.path.join(d, "bcmk.txt"), "chr1\t149129\nchr1\t457293\r\nchr1\t946125\rchr\t12\nchr2\t1572564\nchrchr1\t25\n")
    os.makedirs(os.path.join(d, "out"))


def build_population_af(d):
    shutil.copyfile(os.path.join(HERE, "merge_population", "expected", "files", "out", "pop.vcf"), os.path.join(d, "pop.vcf"))
    with open(os.path.join(d, "pop.vcf"), "a") as f:
        for i, (p, v) in enumerate([("EUR", "0.5"), ("CAS", "1e-05"), ("OCN", "0.10000000000000002"), ("WAS", "nan"),
                                    ("SAS", "０.５"), ("AFR", "1_0.5e-1")]):
            f.write("chr3\t%d\t.\tN\t<X>\t.\tPASS\tSVTYPE=X;%s_AF=%s;Z=1\tGT\t0/1\n" % (i, p, v))
        f.write("chr3\t99\t.\tN\t<X>\t.\tPASS\tEUR_AF=0.9;AFR_AF=0.8;NOEQ\tGT\t1/1\n")
    os.makedirs(os.path.join(d, "out"))


CASES = [
    ("merge_sample", "merger", "-S -i sample_info.txt -x rslt_list.txt -y 15 -p out", build_merge_sample),
    ("merge_population", "merger", "-P -i sample_info.txt -x rslt_list.txt -w 50 -p out -o out/pop.vcf", build_merge_population),
    ("cmp_single", "cmp", "-i calls.txt --bcmk bcmk.txt --exclude exclude.txt --extnd 150 -o out/res", build_cmp_single),
    ("cmp_group", "cmp", "-L --list list.txt --bcmk bcmk.txt -o out/grp", build_cmp_group),
    ("cmp_unicode", "cmp", "-i calls.txt --bcmk bcmk.txt --extnd 20 -o out/u", build_cmp_unicode),
    ("population_af", "af", "pop.vcf out/af.csv", build_population_af),
]


def snapshot(d):
    out = {}
    for root, _, files in os.walk(d):
        for f in files:
            p = os.path.join(root, f)
            with open(p, "rb") as fh:
                out[os.path.relpath(p, d)] = fh.read()
    return out


def main():
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    with open(os.path.join(HERE, "cases.tsv"), "w") as fc:
        for name, tool, args, build in CASES:
            case = os.path.join(HERE, name)
            shutil.rmtree(case, ignore_errors=True)
            os.makedirs(os.path.join(case, "input"))
            build(os.path.join(case, "input"))
            for root, dirs, files in os.walk(os.path.join(case, "input")):
                if not dirs and not files:
                    write(os.path.join(root, ".keep"), "")  # git does not keep empty directories
            run = os.path.join(HERE, "_run")
            shutil.rmtree(run, ignore_errors=True)
            shutil.copytree(os.path.join(case, "input"), run)
            before = snapshot(run)
            p = subprocess.run([PY, os.path.join(PYSRC, SCRIPTS[tool])] + args.split(), cwd=run, env=env,
                               stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            after = snapshot(run)
            exp = os.path.join(case, "expected")
            write(os.path.join(exp, "stdout.txt"), p.stdout)
            write(os.path.join(exp, "exit_code.txt"), "%d\n" % p.returncode)
            for rel, data in sorted(after.items()):
                if before.get(rel) == data:
                    continue
                if len(data) > MAX_PLAIN:
                    write(os.path.join(exp, "files", rel + ".fnv1a64"), fnv1a64(data))
                else:
                    write(os.path.join(exp, "files", rel), data)
            shutil.rmtree(run)
            fc.write("%s\t%s\t%s\n" % (name, tool, args))
            print(name, "exit", p.returncode, "changed files", len([k for k in after if before.get(k) != after[k]]))


if __name__ == "__main__":
    main()
