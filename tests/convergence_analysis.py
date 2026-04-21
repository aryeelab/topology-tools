#!/usr/bin/env python3
"""Convergence analysis for parse_pairs_file() sampling.

Determines the smallest sample_size N at which all QC metrics stabilise to
~3-4 significant figures by comparing consecutive sample sizes (N vs 2N).
When |metric(2N) − metric(N)| / metric(N) < 0.01% for every metric, N is
declared stable.

Usage:
    python tests/convergence_analysis.py --pairs <file.mapped.pairs>

The recommended default sample_size is printed at the end.  Update the
_DEFAULT_SAMPLE_SIZE constant in microc-qc.py accordingly.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO_ROOT / "microc-qc.py"

spec = importlib.util.spec_from_file_location("microc_qc", MODULE_PATH)
microc_qc = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(microc_qc)

SAMPLE_SIZES = [10_000, 50_000, 100_000, 250_000, 500_000, 1_000_000, 2_000_000, 5_000_000]
STABILITY_THRESHOLD = 0.0001   # 0.01% relative change between consecutive doublings
MIN_DOUBLINGS_STABLE = 1       # consecutive stable doublings required


def frag_dist_l1(d1: dict, d2: dict) -> float:
    """Total-variation distance between two (unnormalised) fragment-length histograms."""
    all_keys = set(d1) | set(d2)
    total1 = sum(d1.values()) or 1
    total2 = sum(d2.values()) or 1
    return sum(abs(d1.get(k, 0) / total1 - d2.get(k, 0) / total2) for k in all_keys) / 2


def run_sample(pairs_path: Path, sample_size: int) -> dict:
    print(f"  sample_size={sample_size:>9,} ...", end=" ", flush=True)
    metrics = microc_qc.parse_pairs_file(pairs_path, sample_size=sample_size, progress=False)
    cis_rate = metrics["cis_long_range_pairs"] / max(metrics["non_dup_reads"], 1)
    print(f"cis_rate={cis_rate:.4f}  non_dup={metrics['non_dup_reads']:,}")
    return metrics


def relative_change(a: float, b: float) -> float:
    """|(b-a)/a| clamped to 1.0 when a==0."""
    if a == 0:
        return 0.0 if b == 0 else 1.0
    return abs(b - a) / abs(a)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pairs",
        type=Path,
        required=True,
        help="Path to a .mapped.pairs or .pairs file (plain or gzipped).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=STABILITY_THRESHOLD,
        help=f"Relative-change threshold for declaring stability (default: {STABILITY_THRESHOLD}).",
    )
    args = parser.parse_args()

    if not args.pairs.exists():
        sys.exit(f"ERROR: {args.pairs} not found.")

    threshold = args.threshold
    print(f"Convergence analysis on: {args.pairs}")
    print(f"Stability threshold: {threshold:.4%} relative change between consecutive N and 2N\n")

    results: list[tuple[int, dict]] = []
    for n in SAMPLE_SIZES:
        metrics = run_sample(args.pairs, n)
        results.append((n, metrics))

    print(f"\n{'N':>12}  {'cis_rate':>10}  {'frag_L1':>10}  {'read_len':>10}  "
          f"{'cis_chg':>10}  {'frag_chg':>10}  {'stable':>8}")
    print("-" * 90)

    recommended: int | None = None

    for i, (n, m) in enumerate(results):
        cis_rate = m["cis_long_range_pairs"] / max(m["non_dup_reads"], 1)
        read_len = m.get("read_length") or 0

        if i == 0:
            cis_chg = float("nan")
            frag_chg = float("nan")
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
                recommended = prev_n  # the smaller of the two stable pair

        frag_l1 = frag_dist_l1(
            {} if i == 0 else results[i - 1][1]["fragment_length_distribution"],
            m["fragment_length_distribution"],
        ) if i > 0 else float("nan")

        cis_chg_str = f"{cis_chg:.4%}" if i > 0 else "       —"
        frag_chg_str = f"{frag_chg:.4%}" if i > 0 else "       —"
        frag_l1_str = f"{frag_l1:.4%}" if i > 0 else "       —"

        print(
            f"{n:>12,}  {cis_rate:>10.4f}  {frag_l1_str:>10}  {read_len:>10}  "
            f"{cis_chg_str:>10}  {frag_chg_str:>10}  {stable_str:>8}"
        )

    print()
    if recommended is None:
        print(
            "WARNING: metrics did not converge within the tested range.\n"
            "Try running with a larger file or extending SAMPLE_SIZES.\n"
            f"As a conservative default, recommend: {SAMPLE_SIZES[-1]:,}"
        )
    else:
        print(
            f"RECOMMENDATION: metrics stable at N = {recommended:,}\n"
            f"Update _DEFAULT_SAMPLE_SIZE = {recommended:_} in microc-qc.py"
        )


if __name__ == "__main__":
    main()
