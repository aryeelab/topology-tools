#!/usr/bin/env python3

"""Compute Micro-C QC metrics from pairs, stats, and BAM files."""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, TextIO

try:
    import polars as pl
    _POLARS_AVAILABLE = True
except ImportError:
    _POLARS_AVAILABLE = False

try:
    import pandas as pd
    _PANDAS_AVAILABLE = True
except ImportError:
    _PANDAS_AVAILABLE = False

_PROGRESS_INTERVAL = 500_000
_DEFAULT_SAMPLE_SIZE = 500_000   # empirically determined: see tests/convergence_analysis.py
# On 68M-row chr19 Micro-C data: cis_rate and fragment distribution stable at 250k
# (< 0.5% change when doubling N); 500k chosen as conservative default with margin.


def _progress(msg: str, end: str = "") -> None:
    print(f"\r{msg}", end=end, file=sys.stderr, flush=True)


def open_textfile(path: Path) -> TextIO:
    """Open plain-text or gzipped text files."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")


def infer_consensus_read_length(read_lengths: Counter[int]) -> int | None:
    """Return the modal read length, if any lengths were observed."""
    if not read_lengths:
        return None
    return read_lengths.most_common(1)[0][0]


def _pairs_progress_pct(handle: TextIO, file_size: int) -> str:
    """Return a formatted percent-complete string, or empty string on failure."""
    try:
        # gzip.open("rt") -> TextIOWrapper(GzipFile); plain open -> TextIOWrapper(BufferedReader)
        inner = getattr(handle, "buffer", handle)
        pos = getattr(inner, "fileobj", inner).tell()
        if file_size > 0:
            return f" ({100 * pos // file_size}%)"
    except Exception:
        pass
    return ""


def _count_data_lines(path: Path) -> int:
    """Count non-comment data lines via fast binary byte scan (no per-row parsing)."""
    comment_byte = ord("#")
    newline_byte = ord("\n")
    count = 0
    prev_newline = True  # treat start-of-file as if preceded by a newline
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)  # 1 MB chunks
            if not chunk:
                break
            for b in chunk:
                if prev_newline and b != comment_byte:
                    count += 1
                prev_newline = b == newline_byte
    return count


def _parse_pairs_header(path: Path) -> tuple[list[str], set[str]]:
    """Read comment lines and return (column_names, required_columns_present)."""
    required_columns = {
        "chrom1", "pos1", "chrom2", "pos2", "pair_type",
        "pos51", "pos52", "pos31", "pos32", "read_len1", "read_len2",
    }
    columns: list[str] | None = None
    with open_textfile(path) as fh:
        for line in fh:
            if line.startswith("#columns:"):
                columns = line.strip().split(": ", 1)[1].split()
                break
    if columns is None:
        raise ValueError(f"Pairs file {path} is missing a #columns header")
    missing = required_columns - set(columns)
    if missing:
        raise ValueError(
            f"Pairs file {path} is missing required columns: {', '.join(sorted(missing))}"
        )
    return columns, required_columns


def _parse_pairs_polars(
    path: Path, columns: list[str], cis_distance: int, progress: bool,
    sample_size: int | None = None,
) -> dict:
    """Fastest path: parse pairs file using polars (SIMD-accelerated CSV reader)."""
    needed = ["chrom1", "pos1", "chrom2", "pos2", "read_len1", "read_len2",
              "pos51", "pos52", "pos31", "pos32"]
    int32_cols = [c for c in needed if c not in ("chrom1", "chrom2")]

    if progress:
        _progress("  pairs: reading...")

    df = pl.read_csv(
        str(path),
        separator="\t",
        comment_prefix="#",
        has_header=False,
        new_columns=columns,
        schema_overrides={c: pl.Int32 for c in int32_cols},
    ).select(needed)

    total_rows = len(df)

    if sample_size is not None and total_rows > sample_size:
        stride = max(1, total_rows // sample_size)
        df = df.gather_every(stride)

    scale = total_rows / len(df)

    if progress:
        label = f" (1/{round(scale)}× sample)" if scale > 1.5 else ""
        _progress(f"  pairs: {total_rows:,} reads — aggregating{label}...")

    cis_lr = int(round(
        ((df["chrom1"] == df["chrom2"]) &
         ((df["pos2"] - df["pos1"]).abs() >= cis_distance)).sum() * scale
    ))

    frag1 = (df["pos31"] - df["pos51"]).abs() + 1
    frag2 = (df["pos32"] - df["pos52"]).abs() + 1
    frag_series = pl.concat([frag1.rename("fl"), frag2.rename("fl")])
    frag_vc = frag_series.value_counts().sort("fl")
    raw_fragment_lengths = dict(zip(frag_vc["fl"].to_list(), frag_vc["count"].to_list()))
    if scale > 1.0:
        fragment_lengths = {k: int(round(v * scale)) for k, v in raw_fragment_lengths.items()}
    else:
        fragment_lengths = raw_fragment_lengths

    rl1 = df["read_len1"].filter(df["read_len1"] > 0)
    rl2 = df["read_len2"].filter(df["read_len2"] > 0)
    rl_series = pl.concat([rl1.rename("rl"), rl2.rename("rl")])
    rl_vc = rl_series.value_counts().sort("rl")
    read_lengths: Counter[int] = Counter(
        dict(zip(rl_vc["rl"].to_list(), rl_vc["count"].to_list()))
    )

    if progress:
        _progress(f"  pairs: {total_rows:,} reads (100%)", end="\n")

    return {
        "non_dup_reads": total_rows,
        "cis_long_range_pairs": cis_lr,
        "read_length": infer_consensus_read_length(read_lengths),
        "read_length_distribution": dict(sorted(read_lengths.items())),
        "fragment_length_distribution": dict(sorted(fragment_lengths.items())),
        "fragment_count": 2 * total_rows,
    }


def _parse_pairs_pandas(
    path: Path, columns: list[str], cis_distance: int, progress: bool,
    sample_size: int | None = None,
) -> dict:
    """Fast path: parse pairs file using chunked pandas read_csv."""
    needed = ["chrom1", "pos1", "chrom2", "pos2", "read_len1", "read_len2",
              "pos51", "pos52", "pos31", "pos32"]

    # Pre-count for exact non_dup_reads and stride computation when sampling
    total_rows: int | None = None
    stride = 1
    if sample_size is not None:
        total_rows = _count_data_lines(path)
        stride = max(1, total_rows // sample_size)

    base_row = 0       # global data-row index (before any filtering)
    sampled_rows = 0   # rows actually accumulated for statistics
    cis_long_range_pairs = 0
    fragment_lengths: Counter[int] = Counter()
    read_lengths: Counter[int] = Counter()

    chunks = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=columns,
        usecols=needed,
        chunksize=1_000_000,
        dtype={
            "pos1": "int32", "pos2": "int32",
            "pos51": "int32", "pos52": "int32",
            "pos31": "int32", "pos32": "int32",
            "read_len1": "int16", "read_len2": "int16",
        },
    )

    prev_milestone = 0
    for chunk in chunks:
        chunk_len = len(chunk)

        if stride > 1:
            keep = [i for i in range(chunk_len) if (base_row + i) % stride == 0]
            chunk = chunk.iloc[keep]

        base_row += chunk_len
        sampled_rows += len(chunk)

        rows_seen = total_rows if total_rows is not None else base_row
        if progress:
            milestone = (rows_seen // _PROGRESS_INTERVAL) * _PROGRESS_INTERVAL
            if milestone > prev_milestone:
                _progress(f"  pairs: {rows_seen:,} reads")
                prev_milestone = milestone

        cis = chunk.chrom1 == chunk.chrom2
        cis_long_range_pairs += int(
            (cis & ((chunk.pos2 - chunk.pos1).abs() >= cis_distance)).sum()
        )

        r1 = chunk.read_len1[chunk.read_len1 > 0]
        r2 = chunk.read_len2[chunk.read_len2 > 0]
        for length, count in pd.concat([r1, r2]).value_counts().items():
            read_lengths[int(length)] += int(count)

        frag1 = (chunk.pos31 - chunk.pos51).abs() + 1
        frag2 = (chunk.pos32 - chunk.pos52).abs() + 1
        for length, count in pd.concat([frag1, frag2]).value_counts().items():
            fragment_lengths[int(length)] += int(count)

    non_dup_reads = total_rows if total_rows is not None else base_row
    scale = non_dup_reads / sampled_rows if sampled_rows < non_dup_reads else 1.0

    if progress:
        _progress(f"  pairs: {non_dup_reads:,} reads (100%)", end="\n")

    if scale > 1.0:
        cis_long_range_pairs = int(round(cis_long_range_pairs * scale))
        fragment_lengths = Counter({k: int(round(v * scale)) for k, v in fragment_lengths.items()})

    return {
        "non_dup_reads": non_dup_reads,
        "cis_long_range_pairs": cis_long_range_pairs,
        "read_length": infer_consensus_read_length(read_lengths),
        "read_length_distribution": dict(sorted(read_lengths.items())),
        "fragment_length_distribution": dict(sorted(fragment_lengths.items())),
        "fragment_count": 2 * non_dup_reads,
    }


def _parse_pairs_python(
    path: Path, columns: list[str], cis_distance: int, progress: bool,
    sample_size: int | None = None,
) -> dict:
    """Pure-Python fallback for parse_pairs_file (used when pandas is unavailable)."""
    column_index = {name: idx for idx, name in enumerate(columns)}
    ci = {col: column_index[col] for col in
          ("chrom1", "pos1", "chrom2", "pos2", "read_len1", "read_len2",
           "pos51", "pos52", "pos31", "pos32")}

    # Pre-count for exact total and stride computation when sampling
    total_rows: int | None = None
    stride = 1
    if sample_size is not None:
        total_rows = _count_data_lines(path)
        stride = max(1, total_rows // sample_size)

    data_line = 0  # global data-line counter (always incremented)
    cis_long_range_pairs = 0
    fragment_lengths: Counter[int] = Counter()
    read_lengths: Counter[int] = Counter()
    file_size = path.stat().st_size

    with open_textfile(path) as handle:
        for raw_line in handle:
            if raw_line.startswith("#"):
                continue

            data_line += 1

            if progress and data_line % _PROGRESS_INTERVAL == 0:
                pct = _pairs_progress_pct(handle, file_size)
                _progress(f"  pairs: {data_line:,} reads{pct}")

            # Skip non-sampled rows (stride == 1 means keep all)
            if data_line % stride != 0:
                continue

            fields = raw_line.rstrip("\n").split("\t")

            chrom1 = fields[ci["chrom1"]]
            chrom2 = fields[ci["chrom2"]]
            pos1 = int(fields[ci["pos1"]])
            pos2 = int(fields[ci["pos2"]])
            if chrom1 == chrom2 and abs(pos2 - pos1) >= cis_distance:
                cis_long_range_pairs += 1

            read_len1 = int(fields[ci["read_len1"]])
            read_len2 = int(fields[ci["read_len2"]])
            if read_len1 > 0:
                read_lengths[read_len1] += 1
            if read_len2 > 0:
                read_lengths[read_len2] += 1

            pos51 = int(fields[ci["pos51"]])
            pos52 = int(fields[ci["pos52"]])
            pos31 = int(fields[ci["pos31"]])
            pos32 = int(fields[ci["pos32"]])
            fragment_lengths[abs(pos31 - pos51) + 1] += 1
            fragment_lengths[abs(pos32 - pos52) + 1] += 1

    non_dup_reads = total_rows if total_rows is not None else data_line
    sampled = data_line // stride  # number of rows actually processed
    scale = non_dup_reads / sampled if sampled < non_dup_reads else 1.0

    if progress:
        _progress(f"  pairs: {non_dup_reads:,} reads (100%)", end="\n")

    if scale > 1.0:
        cis_long_range_pairs = int(round(cis_long_range_pairs * scale))
        fragment_lengths = Counter({k: int(round(v * scale)) for k, v in fragment_lengths.items()})

    return {
        "non_dup_reads": non_dup_reads,
        "cis_long_range_pairs": cis_long_range_pairs,
        "read_length": infer_consensus_read_length(read_lengths),
        "read_length_distribution": dict(sorted(read_lengths.items())),
        "fragment_length_distribution": dict(sorted(fragment_lengths.items())),
        "fragment_count": 2 * non_dup_reads,
    }


def parse_pairs_file(
    path: Path,
    cis_distance: int = 10_000,
    progress: bool = True,
    sample_size: int | None = None,
) -> dict:
    """Parse a pairs/pairsam file and compute QC metrics.

    Uses the fastest available library:
      polars  (install: pip install polars)  — ~30× vs pure Python
      pandas  (install: pip install pandas)  — ~4× vs pure Python
      pure Python                            — always available, no extra deps

    When sample_size is set, non_dup_reads is always exact (full file count) while
    rate/distribution metrics are estimated from a stride sample of sample_size reads.
    Set sample_size=0 to disable sampling and read every row.
    """
    if sample_size == 0:
        sample_size = None
    columns, _ = _parse_pairs_header(path)
    if _POLARS_AVAILABLE:
        return _parse_pairs_polars(path, columns, cis_distance, progress, sample_size)
    if _PANDAS_AVAILABLE:
        return _parse_pairs_pandas(path, columns, cis_distance, progress, sample_size)
    return _parse_pairs_python(path, columns, cis_distance, progress, sample_size)


def is_unique_alignment(record, unique_mapq_min: int) -> bool:
    """Determine unique mapping directly from BAM fields."""
    if record.is_unmapped:
        return False
    if record.has_tag("NH"):
        return record.get_tag("NH") == 1
    return record.mapping_quality >= unique_mapq_min


def parse_bam_file(path: Path, unique_mapq_min: int = 20, progress: bool = True) -> dict:
    """Compute read-based metrics from BAM using pysam."""
    try:
        import pysam
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "pysam is required to read BAM files. "
            "Install it with conda or pip, or omit --bam."
        ) from exc

    total_reads = 0
    unique_reads = 0
    read_lengths: Counter[int] = Counter()

    with pysam.AlignmentFile(str(path), "rb") as bam_file:
        for record in bam_file.fetch(until_eof=True):
            if record.is_secondary or record.is_supplementary:
                continue

            total_reads += 1
            if record.query_length:
                read_lengths[record.query_length] += 1
            if is_unique_alignment(record, unique_mapq_min):
                unique_reads += 1
            if progress and total_reads % _PROGRESS_INTERVAL == 0:
                _progress(f"  bam:   {total_reads:,} reads")

    if progress:
        _progress(f"  bam:   {total_reads:,} reads", end="\n")

    return {
        "total_reads": total_reads,
        "unique_reads": unique_reads,
        "read_length": infer_consensus_read_length(read_lengths),
        "read_length_distribution": dict(sorted(read_lengths.items())),
        "unique_mapq_min": unique_mapq_min,
    }


def infer_sample_id(pairs_path: Path, bam_path: Path | None = None) -> str:
    """Infer a readable sample id from the input file names."""
    candidates = [pairs_path]
    if bam_path is not None:
        candidates.append(bam_path)

    for candidate in candidates:
        name = candidate.name
        for suffix in (
            ".mapped.pairsam.gz",
            ".mapped.pairsam",
            ".pairsam.gz",
            ".pairsam",
            ".mapped.pairs.gz",
            ".mapped.pairs",
            ".pairs.gz",
            ".pairs",
            ".bam",
        ):
            if name.endswith(suffix):
                return name[: -len(suffix)]
    return pairs_path.stem


def discover_samples(directory: Path) -> tuple[list[str], list[Path], list[Path]]:
    """Discover matching pairs and BAM files in a directory tree."""
    if not directory.is_dir():
        raise ValueError(f"--dir must point to a directory: {directory}")

    pairs_by_sample: dict[str, Path] = {}
    bam_by_sample: dict[str, Path] = {}

    for path in sorted(directory.rglob("*")):
        if not path.is_file():
            continue

        name = path.name
        if name.endswith(
            (
                ".mapped.pairsam.gz",
                ".mapped.pairsam",
                ".pairsam.gz",
                ".pairsam",
                ".mapped.pairs.gz",
                ".mapped.pairs",
                ".pairs.gz",
                ".pairs",
            )
        ):
            sample_id = infer_sample_id(path)
            if sample_id in pairs_by_sample and pairs_by_sample[sample_id] != path:
                raise ValueError(
                    f"Duplicate pairs files found for sample {sample_id}: "
                    f"{pairs_by_sample[sample_id]} and {path}"
                )
            pairs_by_sample[sample_id] = path
        elif name.endswith(".bam"):
            sample_id = infer_sample_id(path)
            if sample_id in bam_by_sample and bam_by_sample[sample_id] != path:
                raise ValueError(
                    f"Duplicate BAM files found for sample {sample_id}: "
                    f"{bam_by_sample[sample_id]} and {path}"
                )
            bam_by_sample[sample_id] = path

    sample_ids = sorted(set(pairs_by_sample) & set(bam_by_sample))
    if not sample_ids:
        raise ValueError(f"No matching BAM/pairs samples found in {directory}")

    missing_pairs = sorted(set(bam_by_sample) - set(pairs_by_sample))
    missing_bams = sorted(set(pairs_by_sample) - set(bam_by_sample))
    if missing_pairs or missing_bams:
        problems = []
        if missing_pairs:
            problems.append(f"missing pairs for: {', '.join(missing_pairs)}")
        if missing_bams:
            problems.append(f"missing BAMs for: {', '.join(missing_bams)}")
        raise ValueError(f"Incomplete samples in {directory}: {'; '.join(problems)}")

    pairs_paths = [pairs_by_sample[sample_id] for sample_id in sample_ids]
    bam_paths = [bam_by_sample[sample_id] for sample_id in sample_ids]
    return sample_ids, pairs_paths, bam_paths


def build_summary(
    sample_id: str,
    pairs_metrics: dict,
    bam_metrics: dict,
    cis_distance: int,
    pairs_path: Path,
    bam_path: Path | None,
    sample_size: int | None = None,
) -> dict:
    """Assemble the final QC summary document."""
    read_length = pairs_metrics["read_length"]
    read_length_source = "pairs"
    if bam_metrics and bam_metrics.get("read_length") is not None:
        read_length = bam_metrics["read_length"]
        read_length_source = "bam"

    sources: dict = {
        "read_length": read_length_source,
        "total_reads": "bam",
        "unique_reads": (
            f"bam(primary reads; NH==1 when present, else MAPQ>={bam_metrics['unique_mapq_min']})"
        ),
        "non_dup_reads": "pairs(exact line count)",
        "cis_long_range_pairs": f"pairs(same chromosome and distance >= {cis_distance})",
        "fragment_length_distribution": "pairs(abs(pos31-pos51)+1, abs(pos32-pos52)+1)",
    }
    if sample_size:
        sources["sampling"] = (
            f"stride sample of {sample_size:,} from {pairs_metrics['non_dup_reads']:,} reads"
        )

    return {
        "sample_id": sample_id,
        "inputs": {
            "pairs": str(pairs_path),
            "bam": str(bam_path) if bam_path else None,
        },
        "metrics": {
            "read_length": read_length,
            "total_reads": bam_metrics["total_reads"],
            "unique_reads": bam_metrics["unique_reads"],
            "non_dup_reads": pairs_metrics["non_dup_reads"],
            "cis_long_range_pairs": pairs_metrics["cis_long_range_pairs"],
        },
        "fragment_length_distribution": pairs_metrics["fragment_length_distribution"],
        "sources": sources,
    }


def build_batch_summary(
    sample_ids: list[str],
    pairs_paths: list[Path],
    bam_paths: list[Path],
    cis_distance: int,
    unique_mapq_min: int,
    sample_size: int | None = None,
) -> dict:
    """Assemble per-sample summaries for a batch run."""
    n = len(sample_ids)
    print(f"Processing {n} sample{'s' if n != 1 else ''}...", file=sys.stderr)
    samples = []
    for i, (sample_id, pairs_path, bam_path) in enumerate(
        zip(sample_ids, pairs_paths, bam_paths), start=1
    ):
        print(f"[{i}/{n}] {sample_id}", file=sys.stderr)
        pairs_metrics = parse_pairs_file(
            pairs_path, cis_distance=cis_distance, sample_size=sample_size
        )
        bam_metrics = parse_bam_file(bam_path, unique_mapq_min=unique_mapq_min)
        samples.append(
            build_summary(
                sample_id=sample_id,
                pairs_metrics=pairs_metrics,
                bam_metrics=bam_metrics,
                cis_distance=cis_distance,
                pairs_path=pairs_path,
                bam_path=bam_path,
                sample_size=sample_size,
            )
        )
    return {"samples": samples}


def write_sample_summaries(samples: list[dict], output_path: Path | None) -> None:
    """Write one JSON document per sample."""
    if output_path is None:
        output_path = Path.cwd()

    if len(samples) == 1 and output_path.suffix == ".json":
        output_path.write_text(json.dumps(samples[0], indent=2, sort_keys=True) + "\n")
        return

    output_path.mkdir(parents=True, exist_ok=True)
    for sample in samples:
        sample_output = output_path / f"{sample['sample_id']}.qc.json"
        sample_output.write_text(json.dumps(sample, indent=2, sort_keys=True) + "\n")


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute Micro-C QC metrics from parallel lists of pairs and BAM files."
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--pairs",
        nargs="+",
        type=Path,
        help="Pairs or pairsam files. Gzipped input is supported.",
    )
    input_group.add_argument(
        "--dir",
        type=Path,
        help="Directory containing matching BAM and pairs files to discover automatically.",
    )
    parser.add_argument(
        "--bams",
        type=Path,
        nargs="+",
        help="BAM files matching --pairs. Required unless --dir is used.",
    )
    parser.add_argument(
        "--sample-ids",
        type=str,
        nargs="+",
        help="Optional sample ids. Must match the number of BAMs and pairs.",
    )
    parser.add_argument(
        "--unique-mapq-min",
        type=int,
        default=20,
        help=(
            "Fallback MAPQ cutoff for unique BAM mapping when NH tags are absent. "
            "Default: 20."
        ),
    )
    parser.add_argument(
        "--cis-distance",
        type=int,
        default=10_000,
        help="Minimum same-chromosome distance for cis long-range pairs. Default: 10000.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        help=(
            "Optional output path for JSON summaries. "
            "For one sample this may be a .json file or a directory; "
            "for multiple samples this must be a directory. "
            "Defaults to the current working directory."
        ),
    )
    parser.add_argument(
        "--sample-size",
        type=int,
        default=_DEFAULT_SAMPLE_SIZE,
        metavar="N",
        help=(
            "Stride-sample N read pairs from each pairs file for rate and distribution "
            "metrics. non_dup_reads is always exact (full file count). "
            f"Default: {_DEFAULT_SAMPLE_SIZE:,}. Set to 0 to read every row."
        ),
    )
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="JSON indentation for the summary output. Default: 2.",
    )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)

    if args.dir is not None:
        sample_ids, pairs_paths, bam_paths = discover_samples(args.dir)
        if args.sample_ids is not None:
            raise ValueError("--sample-ids cannot be used with --dir")
    else:
        if args.bams is None:
            raise ValueError("--bams is required unless --dir is used")
        if len(args.pairs) != len(args.bams):
            raise ValueError("--pairs and --bams must have the same number of entries")
        if args.sample_ids is not None and len(args.sample_ids) != len(args.pairs):
            raise ValueError("--sample-ids must match the number of pairs/BAM files")

        pairs_paths = args.pairs
        bam_paths = args.bams
        sample_ids = args.sample_ids
        if sample_ids is None:
            sample_ids = [
                infer_sample_id(pairs_path=pairs_path, bam_path=bam_path)
                for pairs_path, bam_path in zip(pairs_paths, bam_paths)
            ]

    summary = build_batch_summary(
        sample_ids=sample_ids,
        pairs_paths=pairs_paths,
        bam_paths=bam_paths,
        cis_distance=args.cis_distance,
        unique_mapq_min=args.unique_mapq_min,
        sample_size=args.sample_size if args.sample_size != 0 else None,
    )

    write_sample_summaries(summary["samples"], args.out)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
