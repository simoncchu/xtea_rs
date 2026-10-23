#!/usr/bin/env python3
"""Build a synthetic test set for xTea (Python vs Rust parity tests).

Creates, under the output folder:
  ref/genome.fa              synthetic reference with embedded L1/Alu/SVA copies
  rep_lib_annotation/        a repeat library with the same layout as the real one
                             (RepeatMasker .out files, consensus copied from the real library)
  gene.gff3                  synthetic GENCODE-like gene annotation
  truth.tsv                  simulated non-reference insertions
  reads_1.fq / reads_2.fq    simulated 2x150bp reads of a diploid sample

Usage: simulate.py <real_rep_lib_dir> <out_dir> [--cov 30] [--seed 1]
"""
import argparse
import os
import random
import shutil
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs = {}
    name = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def rand_seq(rng, n, gc=0.41):
    a = (1 - gc) / 2
    c = gc / 2
    return "".join(rng.choices("ACGT", weights=[a, c, c, a], k=n))


def mutate(rng, seq, rate):
    out = []
    for ch in seq:
        r = rng.random()
        if r < rate * 0.9:
            out.append(rng.choice([b for b in "ACGT" if b != ch]))
        elif r < rate * 0.95:
            continue  # deletion
        elif r < rate:
            out.append(ch)
            out.append(rng.choice("ACGT"))  # insertion
        else:
            out.append(ch)
    return "".join(out)


class RefBuilder:
    """Builds a chromosome as background + embedded TE copies, tracking annotation."""

    def __init__(self, rng, name, length):
        self.rng = rng
        self.name = name
        self.length = length
        self.pieces = []  # (start_in_background, te_record) inserted at sorted positions
        self.ann = []

    def build(self, copies):
        # copies: list of dicts with seq, subfamily, family, div, strand, cns_start, cns_end, cns_len
        rng = self.rng
        n = len(copies)
        bg_len = self.length - sum(len(c["seq"]) for c in copies)
        pts = sorted(rng.sample(range(50000, bg_len - 50000), n))
        bg = rand_seq(rng, bg_len)
        out = []
        last = 0
        cur = 0
        for p, c in zip(pts, copies):
            out.append(bg[last:p])
            cur += p - last
            start = cur
            out.append(c["seq"])
            cur += len(c["seq"])
            end = cur
            last = p
            self.ann.append(dict(c, chrom=self.name, start=start, end=end))
        out.append(bg[last:])
        return "".join(out)


def make_copy(rng, cns_name, cns, subfam, family, div, frag_len=None, strand=None):
    n = len(cns)
    if frag_len is None or frag_len >= n:
        cs, ce = 0, n
    else:
        # 5' truncated: keep the 3' end
        cs, ce = n - frag_len, n
    seq = mutate(rng, cns[cs:ce], div / 100.0)
    if strand is None:
        strand = rng.choice("+C")
    if strand == "C":
        seq = rc(seq)
    return dict(seq=seq, subfamily=subfam, family=family, div=div, strand=strand, cns_start=cs + 1, cns_end=ce, cns_len=n, cns_name=cns_name)


def rmsk_line(a, idx):
    sw = int(len(a["seq"]) * (1 - a["div"] / 100.0) * 9)
    left_q = "(%d)" % 1000
    if a["strand"] == "+":
        cols = [sw, "%.1f" % a["div"], "0.5", "0.5", a["chrom"], a["start"], a["end"], left_q, "+", a["subfamily"], a["family"],
                a["cns_start"], a["cns_end"], "(%d)" % (a["cns_len"] - a["cns_end"]), idx]
    else:
        cols = [sw, "%.1f" % a["div"], "0.5", "0.5", a["chrom"], a["start"], a["end"], left_q, "C", a["subfamily"], a["family"],
                "(%d)" % (a["cns_len"] - a["cns_end"]), a["cns_end"], a["cns_start"], idx]
    return " ".join(str(c) for c in cols) + "\n"


def write_fasta(path, seqs):
    with open(path, "w") as f:
        for name, s in seqs:
            f.write(">%s\n" % name)
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("real_lib")
    ap.add_argument("out")
    ap.add_argument("--cov", type=float, default=30.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--read_len", type=int, default=150)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    out = args.out
    os.makedirs(out + "/ref", exist_ok=True)

    cns_l1 = read_fasta(args.real_lib + "/consensus/LINE1.fa")["LINE1"]
    alu = read_fasta(args.real_lib + "/consensus/ALU.fa")
    sva = read_fasta(args.real_lib + "/consensus/SVA.fa")

    def l1_frag():
        div = rng.uniform(6, 22)
        sub = "L1PA7" if div < 10 else ("L1MA8" if div < 16 else "L1ME3")
        return make_copy(rng, "LINE1", cns_l1, sub, "LINE/L1", round(div, 1), rng.randint(600, 3200))

    def l1_full():
        return make_copy(rng, "LINE1", cns_l1, "L1HS", "LINE/L1", round(rng.uniform(0.2, 1.2), 1))

    def alu_copy():
        r = rng.random()
        if r < 0.25:
            name, div = rng.choice(["AluY", "AluYa5", "AluYb8"]), rng.uniform(1.5, 5)
        elif r < 0.75:
            name, div = rng.choice(["AluSx", "AluSz", "AluSp"]), rng.uniform(7, 12)
        else:
            name, div = rng.choice(["AluJb", "AluJo"]), rng.uniform(14, 20)
        return make_copy(rng, name, alu[name], name, "SINE/Alu", round(div, 1))

    def sva_copy(full):
        name = rng.choice(["SVA_E", "SVA_F", "SVA_D"])
        s = sva[name]
        fl = None if full else rng.randint(600, len(s) - 100)
        return make_copy(rng, name, s, name, "Retroposon/SVA", round(rng.uniform(0.5, 4) if full else rng.uniform(3, 9), 1), fl)

    layout = {"chr1": (6000000, 35, 6, 260, 5, 3), "chr2": (4000000, 25, 4, 180, 4, 2)}
    chroms = []
    annots = []
    for cname, (clen, n_frag, n_full, n_alu, n_sva_frag, n_sva_full) in layout.items():
        copies = [l1_frag() for _ in range(n_frag)] + [l1_full() for _ in range(n_full)] + [alu_copy() for _ in range(n_alu)]
        copies += [sva_copy(False) for _ in range(n_sva_frag)] + [sva_copy(True) for _ in range(n_sva_full)]
        rng.shuffle(copies)
        rb = RefBuilder(rng, cname, clen)
        seq = rb.build(copies)
        chroms.append((cname, seq))
        annots.extend(rb.ann)
    chroms.append(("chrM", rand_seq(rng, 16569)))
    chroms.append(("chr1_KI270706v1_random", rand_seq(rng, 60000)))
    ref = dict(chroms)
    write_fasta(out + "/ref/genome.fa", chroms)

    # ------------------------------------------------------------------ repeat library
    lib = out + "/rep_lib_annotation"
    for d in ["LINE/hg38", "Alu/hg38", "SVA/hg38", "HERV/hg38", "blacklist/hg38"]:
        os.makedirs(lib + "/" + d, exist_ok=True)
    if os.path.exists(lib + "/consensus"):
        shutil.rmtree(lib + "/consensus")
    shutil.copytree(args.real_lib + "/consensus", lib + "/consensus")
    annots.sort(key=lambda a: (a["chrom"], a["start"]))
    l1 = [a for a in annots if a["family"] == "LINE/L1"]
    al = [a for a in annots if a["family"] == "SINE/Alu"]
    sv = [a for a in annots if a["family"] == "Retroposon/SVA"]
    with open(lib + "/LINE/hg38/hg38_L1.fa.out", "w") as f:
        for i, a in enumerate(l1):
            f.write(rmsk_line(a, i + 1))
    with open(lib + "/LINE/hg38/hg38_L1_larger_500_with_all_L1HS.out", "w") as f:
        for i, a in enumerate(l1):
            if a["end"] - a["start"] > 500 or a["subfamily"] == "L1HS":
                f.write(rmsk_line(a, i + 1))
    with open(lib + "/LINE/hg38/L1HS_larger_5K.out", "w") as f:
        for i, a in enumerate(l1):
            if a["subfamily"] == "L1HS" and a["end"] - a["start"] > 5000:
                f.write(rmsk_line(a, i + 1))
    with open(lib + "/Alu/hg38/hg38_Alu.out", "w") as f:
        for i, a in enumerate(al):
            f.write(rmsk_line(a, i + 1))
    with open(lib + "/SVA/hg38/hg38_SVA.out", "w") as f:
        for i, a in enumerate(sv):
            f.write(rmsk_line(a, i + 1))
    with open(lib + "/SVA/hg38/SVA_full_length.out", "w") as f:
        for i, a in enumerate(sv):
            if a["end"] - a["start"] > 1200:
                f.write(rmsk_line(a, i + 1))
    with open(lib + "/blacklist/hg38/centromere.bed", "w") as f:
        f.write("chr1\t3000000\t3020000\tcentromere\n")
        f.write("chr2\t2000000\t2015000\tcentromere\n")

    # ------------------------------------------------------------------ gene annotation
    with open(out + "/gene.gff3", "w") as f:
        f.write("##gff-version 3\n#description: synthetic\n")
        gid = 0
        for cname, clen in [("chr1", 6000000), ("chr2", 4000000)]:
            pos = 100000
            while pos < clen - 200000:
                gid += 1
                glen = rng.randint(20000, 120000)
                strand = rng.choice("+-")
                g = "ENSG%011d.1" % gid
                t = "ENST%011d.1" % gid
                gname = "GENE%d" % gid
                gtype = rng.choice(["protein_coding", "lncRNA", "protein_coding"])
                f.write("%s\tHAVANA\tgene\t%d\t%d\t.\t%s\t.\tID=%s;gene_id=%s;gene_type=%s;gene_name=%s;level=2\n" % (cname, pos, pos + glen, strand, g, g, gtype, gname))
                f.write("%s\tHAVANA\ttranscript\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s;gene_id=%s;transcript_id=%s;gene_type=%s;gene_name=%s\n" % (cname, pos, pos + glen, strand, t, g, g, t, gtype, gname))
                nex = rng.randint(2, 8)
                ex_pos = sorted(rng.sample(range(pos, pos + glen - 300), nex))
                for k, ep in enumerate(ex_pos):
                    f.write("%s\tHAVANA\texon\t%d\t%d\t.\t%s\t.\tID=exon:%s:%d;Parent=%s;gene_id=%s;transcript_id=%s;gene_type=%s;gene_name=%s;exon_number=%d\n" % (cname, ep, ep + rng.randint(80, 290), strand, t, k + 1, t, g, t, gtype, gname, k + 1))
                    if gtype == "protein_coding" and k == 0:
                        f.write("%s\tHAVANA\tUTR\t%d\t%d\t.\t%s\t.\tID=UTR5:%s;Parent=%s;gene_id=%s;transcript_id=%s;gene_type=%s;gene_name=%s\n" % (cname, ep, ep + 50, strand, t, t, g, t, gtype, gname))
                pos += glen + rng.randint(20000, 150000)

    # ------------------------------------------------------------------ sample haplotypes
    full_l1 = [a for a in l1 if a["subfamily"] == "L1HS"]
    events = []
    used = []

    def pick_pos(chrom):
        clen = len(ref[chrom])
        while True:
            p = rng.randint(80000, clen - 80000)
            if all(not (c == chrom and abs(p - q) < 6000) for c, q in used):
                used.append((chrom, p))
                return p

    def pick_pos_in_repeat(chrom):
        cands = [a for a in annots if a["chrom"] == chrom and a["family"] in ("LINE/L1", "SINE/Alu")]
        while True:
            a = rng.choice(cands)
            p = rng.randint(a["start"] + 20, a["end"] - 20)
            if all(not (c == chrom and abs(p - q) < 6000) for c, q in used):
                used.append((chrom, p))
                return p

    def transduced_flank(src, n):
        # 3' downstream flank of the source copy (in the copy's orientation)
        s = ref[src["chrom"]]
        if src["strand"] == "+":
            return s[src["end"]:src["end"] + n]
        return rc(s[src["start"] - n:src["start"]])

    def add(kind, chrom, pos, body, sub, strand, tsd_len, polya_len, zyg, note):
        events.append(dict(kind=kind, chrom=chrom, pos=pos, body=body, sub=sub, strand=strand, tsd=tsd_len, polya=polya_len, zyg=zyg, note=note))

    zygs = ["het1", "het2", "hom"]
    for i in range(14):
        chrom = "chr1" if i % 3 else "chr2"
        pos = pick_pos_in_repeat(chrom) if i % 5 == 4 else pick_pos(chrom)
        r = rng.random()
        if r < 0.25:
            body = cns_l1  # full length
            note = "full"
        else:
            ln = rng.randint(500, 3500)
            body = cns_l1[-ln:]
            note = "trunc%d" % ln
        if i % 4 == 1:
            # 5' inversion via twin priming
            inv = rng.randint(200, 600)
            body = rc(body[:inv]) + body[inv + rng.randint(0, 30):]
            note += "_5inv"
        body = mutate(rng, body, 0.004)
        if i in (2, 6, 9, 12):
            src = rng.choice(full_l1)
            fl = transduced_flank(src, rng.randint(150, 700))
            body = body + "A" * rng.randint(8, 15) + fl
            note += "_td:%s:%d" % (src["chrom"], src["start"])
        add("L1", chrom, pos, body, "L1HS", rng.choice("+-"), rng.randint(7, 18), rng.randint(15, 45), rng.choice(zygs), note)
    for i in range(2):
        chrom = "chr1" if i == 0 else "chr2"
        src = rng.choice(full_l1)
        fl = transduced_flank(src, rng.randint(400, 1200))
        add("L1_orphan", chrom, pick_pos(chrom), fl, "orphan", rng.choice("+-"), rng.randint(8, 16), rng.randint(20, 40), "hom", "orphan_td:%s:%d" % (src["chrom"], src["start"]))
    for i in range(16):
        chrom = "chr1" if i % 2 else "chr2"
        pos = pick_pos_in_repeat(chrom) if i % 6 == 5 else pick_pos(chrom)
        name = rng.choice(["AluYa5", "AluYb8", "AluY", "AluYc1"])
        body = mutate(rng, alu[name], 0.006)
        add("Alu", chrom, pos, body, name, rng.choice("+-"), rng.randint(8, 18), rng.randint(12, 35), rng.choice(zygs), "")
    for i in range(7):
        chrom = "chr1" if i % 2 else "chr2"
        name = rng.choice(["SVA_E", "SVA_F"])
        s = sva[name]
        if i % 3 == 0:
            body = s[-rng.randint(700, len(s) - 50):]
            note = "trunc"
        else:
            body = s
            note = "full"
        body = mutate(rng, body, 0.004)
        add("SVA", chrom, pick_pos(chrom), body, name, rng.choice("+-"), rng.randint(8, 18), rng.randint(15, 40), rng.choice(zygs), note)

    def apply(chrom, hap):
        s = ref[chrom]
        evs = sorted([e for e in events if e["chrom"] == chrom and (e["zyg"] == "hom" or e["zyg"] == "het%d" % hap)], key=lambda e: e["pos"])
        out_s = []
        last = 0
        for e in evs:
            p = e["pos"]
            tsd = s[p - e["tsd"]:p]
            ins = e["body"] + "A" * e["polya"]
            if e["strand"] == "-":
                ins = rc(ins)
            out_s.append(s[last:p])
            out_s.append(ins)
            out_s.append(tsd)
            last = p
        out_s.append(s[last:])
        return "".join(out_s)

    haps = []
    for h in (1, 2):
        hs = []
        for cname, _ in chroms:
            hs.append((cname, apply(cname, h)))
        haps.append(hs)

    with open(out + "/truth.tsv", "w") as f:
        f.write("#chrom\tpos\tkind\tsubfamily\tstrand\tins_len\ttsd_len\tpolyA\tzygosity\tnote\n")
        for e in sorted(events, key=lambda e: (e["chrom"], e["pos"])):
            f.write("%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\n" % (e["chrom"], e["pos"], e["kind"], e["sub"], e["strand"], len(e["body"]), e["tsd"], e["polya"], e["zyg"], e["note"]))

    # ------------------------------------------------------------------ reads
    rl = args.read_len
    qual_chars = "?@ABCDEFGHI"
    n_pair = 0
    with open(out + "/reads_1.fq", "w") as f1, open(out + "/reads_2.fq", "w") as f2:
        for hi, hs in enumerate(haps):
            for cname, s in hs:
                L = len(s)
                npairs = int(L * args.cov / 2.0 / (2 * rl))
                for _ in range(npairs):
                    isz = int(rng.gauss(450, 60))
                    isz = max(rl + 20, min(isz, 1000))
                    st = rng.randint(0, L - isz)
                    frag = s[st:st + isz]
                    if "N" in frag:
                        continue
                    if rng.random() < 0.5:
                        frag = rc(frag)
                    r1 = frag[:rl]
                    r2 = rc(frag[-rl:])
                    reads = []
                    for r in (r1, r2):
                        q = []
                        seq = []
                        for k, ch in enumerate(r):
                            if rng.random() < 0.002 + (0.004 if k > rl - 20 else 0.0):
                                ch = rng.choice([b for b in "ACGT" if b != ch])
                                q.append(rng.choice("+,-./0"))
                            else:
                                q.append(rng.choice(qual_chars))
                            seq.append(ch)
                        reads.append(("".join(seq), "".join(q)))
                    n_pair += 1
                    name = "sim%d_%s_%d" % (hi + 1, cname, n_pair)
                    f1.write("@%s/1\n%s\n+\n%s\n" % (name, reads[0][0], reads[0][1]))
                    f2.write("@%s/2\n%s\n+\n%s\n" % (name, reads[1][0], reads[1][1]))
    print("pairs:", n_pair, file=sys.stderr)


if __name__ == "__main__":
    main()
