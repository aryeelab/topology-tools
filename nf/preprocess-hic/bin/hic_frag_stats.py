#!/usr/bin/env python3
"""Restriction-fragment QC of Hi-C pairs annotated by `pairtools restrict`.

Usage:
  hic_frag_stats.py PAIRS [PAIRS ...]   > counts.tsv   (plain, .gz or - for stdin)
  hic_frag_stats.py --merge COUNTS [COUNTS ...] > counts.tsv

Classifies each pair whose two sides are mapped (chrom != '!'), with HiC-Pro's categories:
  trans                     different chromosomes
  dangling_end              same fragment, inward (strand1 +, strand2 -): unligated DNA
  self_circle               same fragment, outward (-, +): a fragment ligated to itself
  same_fragment_same_strand same fragment, ++ or --
  religation                adjacent fragments, inward: re-ligated at the site that was cut
  adjacent_other            adjacent fragments, other orientations
  cis_valid                 every other cis pair
Pairs must be upper-triangular (pos1 <= pos2 within a chromosome), as pairtools writes them; fragment
numbers (rfrag1, rfrag2) count from the chromosome start. Output: one "key<TAB>count" line per category plus
total (all pairs read), mapped (two mapped sides) and cis. Columns are found by name in the #columns: header.
"""
import gzip
import sys

KEYS = ["total", "mapped", "trans", "cis", "dangling_end", "self_circle", "same_fragment_same_strand",
        "religation", "adjacent_other", "cis_valid"]
NEEDED = ["chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2", "rfrag1", "rfrag2"]


def open_text(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def count_pairs(f, n=None):
    n = n if n is not None else dict.fromkeys(KEYS, 0)
    idx = None
    for line in f:
        if line.startswith("#"):
            if line.startswith("#columns:"):
                cols = line.split()[1:]
                missing = [c for c in NEEDED if c not in cols]
                if missing:
                    sys.exit(f"hic_frag_stats.py: columns {missing} not in the pairs header (run pairtools restrict first)")
                idx = [cols.index(c) for c in NEEDED]
            continue
        if idx is None:
            sys.exit("hic_frag_stats.py: data line before a #columns: header")
        x = line.rstrip("\n").split("\t")
        c1, p1, c2, p2, s1, s2, f1, f2 = (x[i] for i in idx)
        n["total"] += 1
        if c1 == "!" or c2 == "!":
            continue
        n["mapped"] += 1
        if c1 != c2:
            n["trans"] += 1
            continue
        n["cis"] += 1
        if int(p1) > int(p2):
            sys.exit(f"hic_frag_stats.py: pair not upper-triangular: {line.strip()}")
        d = int(f2) - int(f1)
        orient = s1 + s2
        if d == 0:
            if orient == "+-":
                n["dangling_end"] += 1
            elif orient == "-+":
                n["self_circle"] += 1
            else:
                n["same_fragment_same_strand"] += 1
        elif d == 1:
            if orient == "+-":
                n["religation"] += 1
            else:
                n["adjacent_other"] += 1
        else:
            n["cis_valid"] += 1
    return n


def merge(paths):
    n = dict.fromkeys(KEYS, 0)
    for p in paths:
        with open(p) as f:
            for line in f:
                if line.strip():
                    k, v = line.split("\t")
                    if k not in n:
                        sys.exit(f"hic_frag_stats.py: unknown key {k!r} in {p}")
                    n[k] += int(v)
    return n


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        sys.exit(__doc__)
    if argv[0] == "--merge":
        n = merge(argv[1:])
    else:
        n = dict.fromkeys(KEYS, 0)
        for p in argv:
            f = open_text(p)
            count_pairs(f, n)
            if f is not sys.stdin:
                f.close()
    for k in KEYS:
        print(f"{k}\t{n[k]}")


if __name__ == "__main__":
    main()
