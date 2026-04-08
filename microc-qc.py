#!/usr/bin/env python3

"""Compute Micro-C QC metrics from pairs, stats, and BAM files."""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path
from typing import Iterable, TextIO


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


def parse_pairs_file(path: Path, cis_distance: int = 10_000) -> dict:
    """Stream a pairs/pairsam file and compute QC metrics."""
    required_columns = {
        "chrom1",
        "pos1",
        "chrom2",
        "pos2",
        "pair_type",
        "pos51",
        "pos52",
        "pos31",
        "pos32",
        "read_len1",
        "read_len2",
    }

    columns: list[str] | None = None
    column_index: dict[str, int] | None = None
    non_dup_reads = 0
    cis_long_range_pairs = 0
    fragment_lengths: Counter[int] = Counter()
    read_lengths: Counter[int] = Counter()

    with open_textfile(path) as handle:
        for raw_line in handle:
            if raw_line.startswith("#columns:"):
                columns = raw_line.strip().split(": ", 1)[1].split()
                column_index = {name: idx for idx, name in enumerate(columns)}
                missing = required_columns.difference(column_index)
                if missing:
                    missing_str = ", ".join(sorted(missing))
                    raise ValueError(
                        f"Pairs file {path} is missing required columns: {missing_str}"
                    )
                continue

            if raw_line.startswith("#"):
                continue

            if column_index is None:
                raise ValueError(f"Pairs file {path} is missing a #columns header")

            fields = raw_line.rstrip("\n").split("\t")
            non_dup_reads += 1

            chrom1 = fields[column_index["chrom1"]]
            chrom2 = fields[column_index["chrom2"]]
            pos1 = int(fields[column_index["pos1"]])
            pos2 = int(fields[column_index["pos2"]])
            if chrom1 == chrom2 and abs(pos2 - pos1) >= cis_distance:
                cis_long_range_pairs += 1

            read_len1 = int(fields[column_index["read_len1"]])
            read_len2 = int(fields[column_index["read_len2"]])
            if read_len1 > 0:
                read_lengths[read_len1] += 1
            if read_len2 > 0:
                read_lengths[read_len2] += 1

            pos51 = int(fields[column_index["pos51"]])
            pos52 = int(fields[column_index["pos52"]])
            pos31 = int(fields[column_index["pos31"]])
            pos32 = int(fields[column_index["pos32"]])
            fragment_lengths[abs(pos31 - pos51) + 1] += 1
            fragment_lengths[abs(pos32 - pos52) + 1] += 1

    return {
        "non_dup_reads": non_dup_reads,
        "cis_long_range_pairs": cis_long_range_pairs,
        "read_length": infer_consensus_read_length(read_lengths),
        "read_length_distribution": dict(sorted(read_lengths.items())),
        "fragment_length_distribution": dict(sorted(fragment_lengths.items())),
        "fragment_count": sum(fragment_lengths.values()),
    }


def is_unique_alignment(record, unique_mapq_min: int) -> bool:
    """Determine unique mapping directly from BAM fields."""
    if record.is_unmapped:
        return False
    if record.has_tag("NH"):
        return record.get_tag("NH") == 1
    return record.mapping_quality >= unique_mapq_min


def parse_bam_file(path: Path, unique_mapq_min: int = 20) -> dict:
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
) -> dict:
    """Assemble the final QC summary document."""
    read_length = pairs_metrics["read_length"]
    read_length_source = "pairs"
    if bam_metrics and bam_metrics.get("read_length") is not None:
        read_length = bam_metrics["read_length"]
        read_length_source = "bam"

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
        "sources": {
            "read_length": read_length_source,
            "total_reads": "bam",
            "unique_reads": (
                f"bam(primary reads; NH==1 when present, else MAPQ>={bam_metrics['unique_mapq_min']})"
            ),
            "non_dup_reads": "pairs(lines)",
            "cis_long_range_pairs": f"pairs(same chromosome and distance >= {cis_distance})",
            "fragment_length_distribution": "pairs(abs(pos31-pos51)+1, abs(pos32-pos52)+1)",
        },
    }


def build_batch_summary(
    sample_ids: list[str],
    pairs_paths: list[Path],
    bam_paths: list[Path],
    cis_distance: int,
    unique_mapq_min: int,
) -> dict:
    """Assemble per-sample summaries for a batch run."""
    samples = []
    for sample_id, pairs_path, bam_path in zip(sample_ids, pairs_paths, bam_paths):
        pairs_metrics = parse_pairs_file(pairs_path, cis_distance=cis_distance)
        bam_metrics = parse_bam_file(bam_path, unique_mapq_min=unique_mapq_min)
        samples.append(
            build_summary(
                sample_id=sample_id,
                pairs_metrics=pairs_metrics,
                bam_metrics=bam_metrics,
                cis_distance=cis_distance,
                pairs_path=pairs_path,
                bam_path=bam_path,
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
    )

    write_sample_summaries(summary["samples"], args.out)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
