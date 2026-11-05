#!/usr/bin/env python3
"""Synthetic Hi-C test data for preprocess-hic with known answers.

Usage: make_test_data.py OUTDIR

Writes OUTDIR/genome.fa (+ chrom.sizes): two random chromosomes, with a 3 kb segment of chrA copied into chrB
(reads there get MAPQ 0), and two "runs" of 2x50 bp read pairs, OUTDIR/fastq/run{A,B}_{1,2}.fastq.gz, built
from exact genome substrings (no sequencing errors, no ligation junction inside a read). Each read name is
<n>_<category>; expected.json has the number of read pairs per category. Categories (DpnII+HinfI fragments):
  cis_valid, trans, dangling_end (same fragment, inward), self_circle (same fragment, outward),
  religation (adjacent fragments, inward), dup (exact copy of a cis_valid pair; half of them in the other run),
  multi (one read in the duplicated segment), unmapped (both reads random), half_mapped (one read random).
"""
import gzip
import json
import os
import random
import sys

L = 50
SIZES = {"chrA": 400_000, "chrB": 300_000}
REPEAT = ("chrA", 100_000, 103_000, "chrB", 50_000)  # chrA[100000:103000] copied to chrB[50000:53000]
N = {"cis_valid": 300, "trans": 150, "dangling_end": 60, "self_circle": 40, "religation": 30,
     "dup": 50, "multi": 40, "unmapped": 30, "half_mapped": 20}
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def cuts(seq):
    """DpnII (^GATC) and HinfI (G^ANTC) cut positions, written independently of digest_genome.py."""
    c = set()
    for i in range(len(seq) - 4):
        if seq[i:i + 4] == "GATC":
            c.add(i)
        if seq[i] == "G" and seq[i + 1] == "A" and seq[i + 3:i + 5] == "TC":
            c.add(i + 1)
    return sorted(p for p in c if 0 < p < len(seq))


def main(out):
    rng = random.Random(20261105)
    os.makedirs(os.path.join(out, "fastq"), exist_ok=True)
    genome = {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in SIZES.items()}
    c1, s, e, c2, s2 = REPEAT
    genome[c2] = genome[c2][:s2] + genome[c1][s:e] + genome[c2][s2 + e - s:]
    with open(os.path.join(out, "genome.fa"), "w") as f:
        for c, sq in genome.items():
            f.write(f">{c} synthetic\n")
            for i in range(0, len(sq), 60):
                f.write(sq[i:i + 60] + "\n")
    with open(os.path.join(out, "chrom.sizes"), "w") as f:
        for c, n in SIZES.items():
            f.write(f"{c}\t{n}\n")

    def in_repeat(c, a, b):
        return (c == REPEAT[0] and a < REPEAT[2] and b > REPEAT[1]) or \
               (c == REPEAT[3] and a < REPEAT[3 + 1] + REPEAT[2] - REPEAT[1] and b > REPEAT[4])

    frags = {}
    for c, sq in genome.items():
        b = [0] + cuts(sq) + [len(sq)]
        frags[c] = [(i, b[i], b[i + 1]) for i in range(len(b) - 1) if not in_repeat(c, b[i] - 100, b[i + 1] + 100)]
    chroms = list(genome)

    def read(c, pos5, strand):  # pos5: 0-based 5' base
        if strand == "+":
            return genome[c][pos5:pos5 + L]
        return revcomp(genome[c][pos5 - L + 1:pos5 + 1])

    def plus_in(f):   # 5' of a + read whose span is inside fragment f
        return rng.randint(f[1] + 2, f[2] - L - 2)

    def minus_in(f):  # 5' of a - read whose span is inside fragment f
        return rng.randint(f[1] + L + 2, f[2] - 3)

    usable = {c: [f for f in frags[c] if f[2] - f[1] >= L + 10] for c in chroms}
    long_ = [(c, f) for c in chroms for f in frags[c] if f[2] - f[1] >= 300]
    rng.shuffle(long_)
    pairs = []  # (category, (c, pos5, strand), (c, pos5, strand)) for side "left" and "right"

    def side_in(c, f, strand):
        return (c, plus_in(f) if strand == "+" else minus_in(f), strand)

    for _ in range(N["cis_valid"]):
        c = rng.choice(chroms)
        while True:
            fa, fb = rng.sample(usable[c], 2)
            if abs(fa[0] - fb[0]) >= 2:
                break
        pairs.append(("cis_valid", side_in(c, fa, rng.choice("+-")), side_in(c, fb, rng.choice("+-"))))
    for _ in range(N["trans"]):
        pairs.append(("trans", side_in("chrA", rng.choice(usable["chrA"]), rng.choice("+-")),
                      side_in("chrB", rng.choice(usable["chrB"]), rng.choice("+-"))))
    for _ in range(N["dangling_end"]):
        c, f = long_.pop()
        p = rng.randint(f[1] + 2, f[1] + 40)
        q = rng.randint(p + 150, f[2] - 3)
        pairs.append(("dangling_end", (c, p, "+"), (c, q, "-")))
    for _ in range(N["self_circle"]):
        c, f = long_.pop()
        q = rng.randint(f[1] + L + 2, f[1] + L + 40)
        p = rng.randint(q + 100, f[2] - L - 2)
        pairs.append(("self_circle", (c, q, "-"), (c, p, "+")))
    by_index = {c: {f[0]: f for f in usable[c]} for c in chroms}
    n = 0
    while n < N["religation"]:
        c = rng.choice(chroms)
        f = rng.choice(usable[c])
        g = by_index[c].get(f[0] + 1)
        if g is None:
            continue
        pairs.append(("religation", side_in(c, f, "+"), side_in(c, g, "-")))
        n += 1
    valid = [p for p in pairs if p[0] == "cis_valid"]
    for p in rng.sample(valid, N["dup"]):
        pairs.append(("dup",) + p[1:])
    for _ in range(N["multi"]):
        pos = rng.randint(REPEAT[1] + 200, REPEAT[2] - 200)
        pairs.append(("multi", ("chrA", pos, "+"), side_in("chrA", rng.choice(usable["chrA"]), rng.choice("+-"))))
    rand = lambda: "".join(rng.choice("ACGT") for _ in range(L))  # noqa: E731
    for _ in range(N["unmapped"]):
        pairs.append(("unmapped", None, None))
    for _ in range(N["half_mapped"]):
        pairs.append(("half_mapped", side_in("chrB", rng.choice(usable["chrB"]), "+"), None))

    rng.shuffle(pairs)
    runs = {"runA": [], "runB": []}
    for i, (cat, a, b) in enumerate(pairs):
        ra = read(*a) if a else rand()
        rb = read(*b) if b else rand()
        if rng.random() < 0.5:
            ra, rb = rb, ra
        runs["runA" if i % 5 < 3 else "runB"].append((f"{i}_{cat}", ra, rb))
    for run, recs in runs.items():
        for mate in (1, 2):
            with gzip.open(os.path.join(out, "fastq", f"{run}_{mate}.fastq.gz"), "wt") as f:
                for name, r1, r2 in recs:
                    sq = r1 if mate == 1 else r2
                    f.write(f"@{name}\n{sq}\n+\n{'I' * len(sq)}\n")
    nfr = sum(len(cuts(sq)) + 1 for sq in genome.values())
    with open(os.path.join(out, "expected.json"), "w") as f:
        json.dump({"read_pairs": len(pairs), "categories": N, "fragments": nfr}, f, indent=2)


if __name__ == "__main__":
    main(sys.argv[1])
