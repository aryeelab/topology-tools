#!/usr/bin/env python3
"""Summarise the QC of one preprocess-hic sample as JSON, and check that the BAM is lossless.

Usage: hic_qc_summary.py --sample-id SID --dedup-stats SID.dedup.stats.txt --restrict-stats SID.restrict_stats.tsv
                         --flagstat SID.flagstat.tsv [--version V] [-o SID.hic_qc.json]

--dedup-stats: `pairtools dedup --output-stats` (pre-dedup totals); --restrict-stats: hic_frag_stats.py
counts of the deduplicated pairs; --flagstat: `samtools flagstat -O tsv` of the BAM.
Exits 1 (after writing the JSON) when the BAM's primary records are not 2 x the read pairs in the stats:
every read of every pair must be in the BAM.
"""
import argparse
import json
import sys


def read_kv(path):
    d = {}
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                k, v = line.rstrip("\n").split("\t")[:2]
                d[k] = v
    return d


def read_flagstat(path):
    d = {}
    with open(path) as f:
        for line in f:
            if line.strip():
                passed, failed, cat = line.rstrip("\n").split("\t")[:3]
                if passed.isdigit() and failed.isdigit():  # skip the "mapped %" style lines
                    d[cat] = int(passed) + int(failed)
    return d


def frac(a, b):
    return round(a / b, 6) if b else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--dedup-stats", required=True)
    ap.add_argument("--restrict-stats", required=True)
    ap.add_argument("--flagstat", required=True)
    ap.add_argument("--version", default=None)
    ap.add_argument("-o", "--output", default="-")
    a = ap.parse_args(argv)

    s = read_kv(a.dedup_stats)
    i = lambda k: int(float(s.get(k, 0)))  # noqa: E731
    r = {k: int(v) for k, v in read_kv(a.restrict_stats).items()}
    fs = read_flagstat(a.flagstat)

    total, mapped, dups, nodups = i("total"), i("total_mapped"), i("total_dups"), i("total_nodups")
    cis, trans = i("cis"), i("trans")
    primary = fs.get("primary")
    lossless = primary == 2 * total
    rm = r.get("mapped", 0)
    qc = {
        "sample_id": a.sample_id,
        "preprocess_hic_version": a.version,
        "read_pairs": total,
        "unmapped_pairs": i("total_unmapped"),
        "single_sided_pairs": i("total_single_sided_mapped"),
        "mapped_pairs": mapped,
        "frac_mapped": frac(mapped, total),
        "duplicate_pairs": dups,
        "frac_duplicates": frac(dups, mapped),
        "nodup_pairs": nodups,
        "frac_nodup_of_read_pairs": frac(nodups, total),
        "pair_types": {k.split("/", 1)[1]: int(float(v)) for k, v in s.items() if k.startswith("pair_types/")},
        "nodup": {
            "cis": cis, "trans": trans, "frac_cis": frac(cis, nodups), "frac_trans": frac(trans, nodups),
            "cis_1kb+": i("cis_1kb+"), "cis_20kb+": i("cis_20kb+"),
            "frac_cis_1kb+": frac(i("cis_1kb+"), nodups), "frac_cis_20kb+": frac(i("cis_20kb+"), nodups),
            "frac_cis_lt_1kb": frac(cis - i("cis_1kb+"), nodups),
        },
        "restriction_fragments": {
            **r,
            **{f"frac_{k}": frac(r.get(k, 0), rm) for k in
               ["dangling_end", "self_circle", "same_fragment_same_strand", "religation", "adjacent_other",
                "cis_valid", "trans"]},
            "frac_valid": frac(r.get("cis_valid", 0) + r.get("trans", 0), rm),
            "note": "deduplicated pairs; fractions of pairs with both sides mapped; QC only, no pairs removed",
        },
        "bam": {"records": fs.get("total (QC-passed reads + QC-failed reads)"), "primary": primary,
                "secondary": fs.get("secondary"), "supplementary": fs.get("supplementary"),
                "primary_mapped": fs.get("primary mapped"),
                "lossless_check": "pass" if lossless else "FAIL",
                "lossless_rule": "primary records == 2 x read_pairs in the dedup stats"},
    }
    if r.get("mapped") is not None and r["mapped"] != nodups:
        qc["restriction_fragments"]["warning"] = f"restrict stats counted {r['mapped']} mapped pairs, dedup stats {nodups} nodups"
    out = sys.stdout if a.output == "-" else open(a.output, "w")
    out.write(json.dumps(qc, indent=2) + "\n")
    if out is not sys.stdout:
        out.close()
    if not lossless:
        sys.exit(f"hic_qc_summary.py: BAM has {primary} primary records, expected 2 x {total} read pairs")


if __name__ == "__main__":
    main()
