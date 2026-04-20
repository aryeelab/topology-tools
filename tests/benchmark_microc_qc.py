#!/usr/bin/env python3
"""Benchmark and profile parse_pairs_file().

Usage:
    python tests/benchmark_microc_qc.py [--rows N] [--no-profile]

Creates a synthetic pairs file by tiling small-rcmc.mapped.pairs to ~N data
rows (default 10_000_000), then times parse_pairs_file() and optionally
profiles it with cProfile.
"""
from __future__ import annotations

import argparse
import cProfile
import importlib.util
import io
import pstats
import sys
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO_ROOT / "microc-qc.py"

spec = importlib.util.spec_from_file_location("microc_qc", MODULE_PATH)
microc_qc = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(microc_qc)

SOURCE_PAIRS = REPO_ROOT / "test-output" / "small-rcmc.mapped.pairs"
SOURCE_ROWS = 22_110  # data lines in small-rcmc.mapped.pairs


def build_large_pairs(target_rows: int, dest: Path) -> int:
    """Tile SOURCE_PAIRS until dest has >= target_rows data lines."""
    # Read header and data lines from source
    header_lines: list[str] = []
    data_lines: list[str] = []
    with SOURCE_PAIRS.open() as fh:
        for line in fh:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                data_lines.append(line)

    repeats = max(1, -(-target_rows // SOURCE_ROWS))  # ceiling division
    total_rows = repeats * SOURCE_ROWS

    with dest.open("w") as out:
        out.writelines(header_lines)
        for _ in range(repeats):
            out.writelines(data_lines)

    return total_rows


def run_benchmark(pairs_path: Path, label: str) -> float:
    t0 = time.perf_counter()
    result = microc_qc.parse_pairs_file(pairs_path, progress=False)
    elapsed = time.perf_counter() - t0
    rows = result["non_dup_reads"]
    throughput = rows / elapsed / 1_000_000
    print(f"{label}: {rows:,} rows in {elapsed:.2f}s  ({throughput:.2f}M rows/sec)")
    return elapsed


def run_profile(pairs_path: Path) -> None:
    pr = cProfile.Profile()
    pr.enable()
    microc_qc.parse_pairs_file(pairs_path, progress=False)
    pr.disable()

    stream = io.StringIO()
    ps = pstats.Stats(pr, stream=stream).sort_stats("cumulative")
    ps.print_stats(20)
    print(stream.getvalue())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", type=int, default=10_000_000,
                        help="Approximate number of data rows in synthetic file (default: 10M)")
    parser.add_argument("--no-profile", action="store_true",
                        help="Skip cProfile run")
    args = parser.parse_args()

    if not SOURCE_PAIRS.exists():
        sys.exit(f"ERROR: {SOURCE_PAIRS} not found. Generate test fixtures first.")

    with tempfile.NamedTemporaryFile(suffix=".mapped.pairs", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    try:
        print(f"Building synthetic pairs file with ~{args.rows:,} rows...")
        actual_rows = build_large_pairs(args.rows, tmp_path)
        print(f"  → {actual_rows:,} rows written to {tmp_path}\n")

        print("=== Timing ===")
        run_benchmark(tmp_path, "parse_pairs_file")

        if not args.no_profile:
            print("\n=== cProfile (top 20 by cumulative time) ===")
            run_profile(tmp_path)
    finally:
        tmp_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
