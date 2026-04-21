#!/usr/bin/env python3
"""Convergence analysis for parse_pairs_file() sampling.

Determines the smallest sample_size N at which all QC metrics stabilise to
~3-4 significant figures by comparing consecutive sample sizes (N vs 2N).
When |metric(2N) − metric(N)| / metric(N) < 0.01% for every metric, N is
declared stable.

When polars is available the file is loaded ONCE into a DataFrame; all
sample-size evaluations then operate on in-memory subsets (much faster than
re-reading a large file repeatedly).

Usage:
    python tests/convergence_analysis.py --pairs <file.mapped.pairs>

The recommended default sample_size is printed at the end.  Update the
_DEFAULT_SAMPLE_SIZE constant in microc-qc.py accordingly.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO_ROOT / "microc-qc.py"

spec = importlib.util.spec_from_file_location("microc_qc", MODULE_PATH)
microc_qc = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(microc_qc)

SAMPLE_SIZES = [10_000, 50_000, 100_000, 250_000, 500_000, 1_000_000, 2_000_000, 5_000_000]
STABILITY_THRESHOLD = 0.0001   # 0.01% relative change


def frag_dist_l1(d1: dict, d2: dict) -> float:
    """Total-variation distance between two (unnormalised) fragment-length histograms."""
    all_keys = set(d1) | set(d2)
    total1 = sum(d1.values()) or 1
    total2 = sum(d2.values()) or 1
    return sum(abs(d1.get(k, 0) / total1 - d2.get(k, 0) / total2) for k in all_keys) / 2


def relative_change(a: float, b: float) -> float:
    if a == 0:
        return 0.0 if b == 0 else 1.0
    return abs(b - a) / abs(a)


# ---------------------------------------------------------------------------
# Fast in-memory path (polars available)
# ---------------------------------------------------------------------------

def _metrics_from_polars_df(df, total_rows: int, sample_size: int, cis_distance: int) -> dict:
    """Compute QC metrics from a polars DataFrame subsampled to sample_size rows."""
    import polars as pl

    if total_rows > sample_size:
        stride = max(1, total_rows // sample_size)
        sdf = df.gather_every(stride)
    else:
        sdf = df

    scale = total_rows / len(sdf)

    cis_lr = int(round(
        ((sdf["chrom1"] == sdf["chrom2"]) &
         ((sdf["pos2"] - sdf["pos1"]).abs() >= cis_distance)).sum() * scale
    ))

    frag1 = (sdf["pos31"] - sdf["pos51"]).abs() + 1
    frag2 = (sdf["pos32"] - sdf["pos52"]).abs() + 1
    frag_series = pl.concat([frag1.rename("fl"), frag2.rename("fl")])
    frag_vc = frag_series.value_counts().sort("fl")
    raw_fl = dict(zip(frag_vc["fl"].to_list(), frag_vc["count"].to_list()))
    fragment_lengths = {k: int(round(v * scale)) for k, v in raw_fl.items()} if scale > 1 else raw_fl

    rl1 = sdf["read_len1"].filter(sdf["read_len1"] > 0)
    rl2 = sdf["read_len2"].filter(sdf["read_len2"] > 0)
    rl_series = pl.concat([rl1.rename("rl"), rl2.rename("rl")])
    rl_vc = rl_series.value_counts().sort("rl")
    read_lengths: Counter[int] = Counter(
        dict(zip(rl_vc["rl"].to_list(), rl_vc["count"].to_list()))
    )

    return {
        "non_dup_reads": total_rows,
        "cis_long_range_pairs": cis_lr,
        "read_length": microc_qc.infer_consensus_read_length(read_lengths),
        "fragment_length_distribution": dict(sorted(fragment_lengths.items())),
    }


def run_all_polars(pairs_path: Path, cis_distance: int) -> list[tuple[int, dict]]:
    """Load the full file once with polars, then subsample at each N."""
    import polars as pl

    needed = ["chrom1", "pos1", "chrom2", "pos2", "read_len1", "read_len2",
              "pos51", "pos52", "pos31", "pos32"]
    int32_cols = [c for c in needed if c not in ("chrom1", "chrom2")]

    columns, _ = microc_qc._parse_pairs_header(pairs_path)

    t0 = time.perf_counter()
    print("Loading full file into memory...", end=" ", flush=True)
    df = pl.read_csv(
        str(pairs_path),
        separator="\t",
        comment_prefix="#",
        has_header=False,
        new_columns=columns,
        schema_overrides={c: pl.Int32 for c in int32_cols},
    ).select(needed)
    total_rows = len(df)
    elapsed = time.perf_counter() - t0
    print(f"{total_rows:,} rows in {elapsed:.1f}s\n")

    results = []
    for n in SAMPLE_SIZES:
        t0 = time.perf_counter()
        print(f"  sample_size={n:>9,} ...", end=" ", flush=True)
        metrics = _metrics_from_polars_df(df, total_rows, n, cis_distance)
        elapsed = time.perf_counter() - t0
        cis_rate = metrics["cis_long_range_pairs"] / max(total_rows, 1)
        print(f"cis_rate={cis_rate:.4f}  ({elapsed:.2f}s)")
        results.append((n, metrics))
        if total_rows <= n:
            print(f"  (file has {total_rows:,} rows — remaining sample sizes are identical)")
            for remaining in SAMPLE_SIZES[SAMPLE_SIZES.index(n) + 1:]:
                results.append((remaining, metrics))
            break

    return results


# ---------------------------------------------------------------------------
# Fallback: call parse_pairs_file N times (slow for large files)
# ---------------------------------------------------------------------------

def run_all_slow(pairs_path: Path, cis_distance: int) -> list[tuple[int, dict]]:
    results = []
    for n in SAMPLE_SIZES:
        t0 = time.perf_counter()
        print(f"  sample_size={n:>9,} ...", end=" ", flush=True)
        metrics = microc_qc.parse_pairs_file(
            pairs_path, cis_distance=cis_distance, sample_size=n, progress=False
        )
        elapsed = time.perf_counter() - t0
        cis_rate = metrics["cis_long_range_pairs"] / max(metrics["non_dup_reads"], 1)
        print(f"cis_rate={cis_rate:.4f}  non_dup={metrics['non_dup_reads']:,}  ({elapsed:.1f}s)")
        results.append((n, metrics))
    return results


# ---------------------------------------------------------------------------
# Print convergence table and recommendation
# ---------------------------------------------------------------------------

def print_table(results: list[tuple[int, dict]], threshold: float) -> int | None:
    print(f"\n{'N':>12}  {'cis_rate':>10}  {'frag_L1':>10}  {'read_len':>10}  "
          f"{'cis_chg':>10}  {'frag_chg':>10}  {'stable':>8}")
    print("-" * 90)

    recommended: int | None = None

    for i, (n, m) in enumerate(results):
        cis_rate = m["cis_long_range_pairs"] / max(m["non_dup_reads"], 1)
        read_len = m.get("read_length") or 0

        if i == 0:
            cis_chg = frag_chg = float("nan")
            stable_str = "—"
        else:
            prev_n, prev_m = results[i - 1]
            prev_cis_rate = prev_m["cis_long_range_pairs"] / max(prev_m["non_dup_reads"], 1)
            cis_chg = relative_change(prev_cis_rate, cis_rate)
            frag_chg = frag_dist_l1(
                prev_m["fragment_length_distribution"],
                m["fragment_length_distribution"],
            )
            is_stable = cis_chg < threshold and frag_chg < threshold
            stable_str = "YES" if is_stable else "no"
            if is_stable and recommended is None:
                recommended = prev_n

        frag_l1_str = f"{frag_dist_l1(results[i-1][1]['fragment_length_distribution'], m['fragment_length_distribution']):.4%}" if i > 0 else "       —"
        cis_chg_str = f"{cis_chg:.4%}" if i > 0 else "       —"
        frag_chg_str = f"{frag_chg:.4%}" if i > 0 else "       —"

        print(
            f"{n:>12,}  {cis_rate:>10.4f}  {frag_l1_str:>10}  {read_len:>10}  "
            f"{cis_chg_str:>10}  {frag_chg_str:>10}  {stable_str:>8}"
        )

    return recommended


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pairs",
        type=Path,
        required=True,
        help="Path to a .mapped.pairs or .pairs file (plain or gzipped).",
    )
    parser.add_argument(
        "--cis-distance",
        type=int,
        default=10_000,
        help="Minimum same-chromosome distance for cis long-range pairs. Default: 10000.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=STABILITY_THRESHOLD,
        help=f"Relative-change threshold for stability (default: {STABILITY_THRESHOLD}).",
    )
    args = parser.parse_args()

    if not args.pairs.exists():
        sys.exit(f"ERROR: {args.pairs} not found.")

    threshold = args.threshold
    print(f"Convergence analysis on: {args.pairs}")
    print(f"Stability threshold: {threshold:.4%} relative change between consecutive N and 2N\n")

    if microc_qc._POLARS_AVAILABLE:
        results = run_all_polars(args.pairs, args.cis_distance)
    else:
        print("polars not available — falling back to repeated parse_pairs_file calls (slow)\n")
        results = run_all_slow(args.pairs, args.cis_distance)

    recommended = print_table(results, threshold)

    print()
    if recommended is None:
        print(
            "WARNING: metrics did not converge within the tested range.\n"
            "Try running on a larger file or extending SAMPLE_SIZES.\n"
            f"As a conservative default, recommend: {SAMPLE_SIZES[-1]:,}"
        )
    else:
        print(
            f"RECOMMENDATION: metrics stable at N = {recommended:,}\n"
            f"Update _DEFAULT_SAMPLE_SIZE = {recommended:_} in microc-qc.py"
        )


if __name__ == "__main__":
    main()
