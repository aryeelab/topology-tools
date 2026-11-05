#!/usr/bin/env python3
"""In-silico digest of a genome with one or more restriction enzymes, as a BED of fragments.

Usage: digest_genome.py --fasta GENOME.fa --chrom-sizes CHROM.sizes --enzymes DpnII,HinfI [-o frags.bed]

The BED (chrom, start, end; 0-based, half-open) is the input of `pairtools restrict -f`. Chromosomes are
written in chrom-sizes order, and only those listed there; a chromosome in the sizes file but not in the
FASTA, or with a different length, is an error. Each fragment runs from one cut to the next; the first
starts at 0 and the last ends at the chromosome length. Sites are matched on the forward strand only
(all supported enzymes are palindromic, N included) and case-insensitively; N in the genome never matches.
Several enzymes: the union of their cut positions (Arima Hi-C+ = DpnII,HinfI).
"""
import argparse
import re
import sys

# name: (recognition site, cut offset on the forward strand from the start of the site)
ENZYMES = {
    "DpnII": ("GATC", 0),     # ^GATC
    "MboI": ("GATC", 0),      # ^GATC
    "Sau3AI": ("GATC", 0),    # ^GATC
    "HinfI": ("GANTC", 1),    # G^ANTC
    "DdeI": ("CTNAG", 1),     # C^TNAG
    "MseI": ("TTAA", 1),      # T^TAA
    "CviJI": ("RGCY", 2),     # RG^CY
    "NlaIII": ("CATG", 4),    # CATG^
    "HindIII": ("AAGCTT", 1),  # A^AGCTT
    "NcoI": ("CCATGG", 1),    # C^CATGG
    "BglII": ("AGATCT", 1),   # A^GATCT
}
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "N": "[ACGT]", "R": "[AG]", "Y": "[CT]",
         "W": "[AT]", "S": "[CG]", "K": "[GT]", "M": "[AC]"}


def site_regex(site: str) -> re.Pattern:
    # Lookahead so that overlapping sites are all found.
    return re.compile("(?=" + "".join(IUPAC[b] for b in site) + ")")


def read_fasta(path):
    """Yield (name, sequence) with the name cut at the first whitespace. Plain or gzipped FASTA."""
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = gzip.open
    name, chunks = None, []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def cut_positions(seq: str, patterns) -> list:
    seq = seq.upper()
    cuts = set()
    for rx, offset in patterns:
        for m in rx.finditer(seq):
            p = m.start() + offset
            if 0 < p < len(seq):
                cuts.add(p)
    return sorted(cuts)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--enzymes", required=True, help="comma-separated, e.g. DpnII,HinfI; known: " + ",".join(ENZYMES))
    ap.add_argument("-o", "--output", default="-")
    a = ap.parse_args(argv)

    names = [e.strip() for e in a.enzymes.split(",") if e.strip()]
    unknown = [e for e in names if e not in ENZYMES]
    if unknown or not names:
        sys.exit(f"unknown enzyme(s) {unknown}; known: {','.join(ENZYMES)}")
    patterns = [(site_regex(ENZYMES[e][0]), ENZYMES[e][1]) for e in names]

    sizes = []
    with open(a.chrom_sizes) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                c, n = line.split()[:2]
                sizes.append((c, int(n)))
    wanted = dict(sizes)
    frags = {}
    for name, seq in read_fasta(a.fasta):
        if name not in wanted:
            continue
        if len(seq) != wanted[name]:
            sys.exit(f"{name}: length {len(seq)} in FASTA but {wanted[name]} in {a.chrom_sizes}")
        frags[name] = cut_positions(seq, patterns)
    missing = [c for c, _ in sizes if c not in frags]
    if missing:
        sys.exit(f"chromosomes in {a.chrom_sizes} but not in {a.fasta}: {missing}")

    out = sys.stdout if a.output == "-" else open(a.output, "w")
    n = 0
    for c, length in sizes:
        bounds = [0] + frags[c] + [length]
        for s, e in zip(bounds[:-1], bounds[1:]):
            out.write(f"{c}\t{s}\t{e}\n")
            n += 1
    if out is not sys.stdout:
        out.close()
    print(f"digest_genome.py: {','.join(names)}: {n} fragments on {len(sizes)} chromosomes", file=sys.stderr)


if __name__ == "__main__":
    main()
