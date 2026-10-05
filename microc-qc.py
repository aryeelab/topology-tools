#!/usr/bin/env python3

"""Compute Micro-C QC metrics from pairs, stats, and BAM files.

Memory: in the default exact mode every pairs backend (polars, pandas, pure Python; see
--backend) reads the file in fixed-size batches, so peak memory is bounded (about 1-2 GB per sample)
regardless of file size. Runtime scales with file size (pandas: ~9 min per 100 GB of pairs).
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
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

__version__ = "1.0.0"  # microc-qc release; tagged microc-qc/v<version> (see RELEASING.md)

_PROGRESS_INTERVAL = 500_000
_DEFAULT_SAMPLE_SIZE = 0  # 0 = exact (no sampling); streaming aggregations make exact as
# fast and memory-efficient as sampling, so there is no longer a speed/accuracy tradeoff.
_BAM_SAMPLE_READS = 500_000  # primary reads to scan for unique-fraction and read-length estimate
# fragment_length_distribution counts only read ends on primary nuclear chromosomes: numbered
# autosomes and X, with or without a "chr" prefix, including subregion contigs ("chr10:1-2").
# chrM, chrY, unplaced/random/alt scaffolds and decoys are excluded.
_DEFAULT_FRAGLEN_CHROM_PATTERN = r"^(chr)?([0-9]+|X)(:[0-9]+-[0-9]+)?$"
_CHRM_NAMES = ("chrM", "chrMT", "M", "MT")


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
    """Count non-comment data lines using fast chunked binary reads."""
    opener = gzip.open if str(path).endswith(".gz") else open
    count = 0
    buf = b""
    with opener(path, "rb") as fh:
        while True:
            chunk = fh.read(64 << 20)  # 64 MB chunks
            if not chunk:
                break
            lines = (buf + chunk).split(b"\n")
            buf = lines.pop()  # last partial line (no trailing newline yet)
            for line in lines:
                if line and line[0:1] != b"#":
                    count += 1
    if buf and buf[0:1] != b"#":
        count += 1
    return count


def _collect_stride_lines(path: Path, stride: int, progress: bool) -> bytes:
    """Read every stride-th data line; return as concatenated bytes for polars.

    Uses 64 MB binary chunks + C-level split() so Python only loops over the
    resulting list — much faster than line-by-line iteration.  Memory is
    O(sample_size * bytes_per_row), never proportional to the full file size.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    buf = b""
    data_line = 0
    sampled: list[bytes] = []
    with opener(path, "rb") as fh:
        while True:
            chunk = fh.read(64 << 20)
            if not chunk:
                break
            lines = (buf + chunk).split(b"\n")
            buf = lines.pop()
            for line in lines:
                if not line or line[0:1] == b"#":
                    continue
                data_line += 1
                if progress and data_line % _PROGRESS_INTERVAL == 0:
                    _progress(f"  pairs: {data_line:,} reads")
                if data_line % stride == 0:
                    sampled.append(line + b"\n")
    # handle last partial line (file not ending in '\n')
    if buf and buf[0:1] != b"#":
        data_line += 1
        if data_line % stride == 0:
            sampled.append(buf + b"\n")
    return b"".join(sampled)


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


def _polars_metrics_from_df(
    df, total_rows: int, scale: float, cis_distance: int, progress: bool,
    chrom_pattern: str = _DEFAULT_FRAGLEN_CHROM_PATTERN,
) -> dict:
    """Compute QC metrics from an already-loaded (possibly sampled) polars DataFrame."""
    if progress:
        label = f" (1/{round(scale)}× sample)" if scale > 1.5 else ""
        _progress(f"  pairs: {total_rows:,} reads — aggregating{label}...")

    cis_lr = int(round(
        ((df["chrom1"] == df["chrom2"]) &
         ((df["pos2"] - df["pos1"]).abs() >= cis_distance)).sum() * scale
    ))

    raw_fl: Counter[int] = Counter()
    chroms: set[str] = set()
    raw_chrm = _polars_fraglen_into(raw_fl, chroms, df, chrom_pattern)
    fragment_lengths = (
        {k: int(round(v * scale)) for k, v in raw_fl.items()} if scale > 1.0 else dict(raw_fl)
    )

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
        "fragment_length_chromosomes": sorted(chroms),
        "fragment_count": 2 * total_rows,
        "chrM_fragment_count": int(round(raw_chrm * scale)),
    }


def _polars_fraglen_into(counter: Counter, chroms: set, df, chrom_pattern: str) -> int:
    """Add per-read-end fragment lengths on matching chromosomes to `counter`.

    Each read end is kept or dropped by its own chromosome (chrom1 for end 1, chrom2 for
    end 2). Matching chromosome names are added to `chroms`; returns the number of chrM ends.
    """
    c1 = df["chrom1"].cast(pl.Utf8)
    c2 = df["chrom2"].cast(pl.Utf8)
    k1 = c1.str.contains(chrom_pattern)
    k2 = c2.str.contains(chrom_pattern)
    frag = pl.concat([
        ((df["pos31"] - df["pos51"]).abs() + 1).filter(k1).rename("v"),
        ((df["pos32"] - df["pos52"]).abs() + 1).filter(k2).rename("v"),
    ])
    _value_counts_into(counter, frag)
    chroms.update(c1.filter(k1).unique().to_list())
    chroms.update(c2.filter(k2).unique().to_list())
    return int(c1.is_in(_CHRM_NAMES).sum() + c2.is_in(_CHRM_NAMES).sum())


def _parse_pairs_polars(
    path: Path, columns: list[str], cis_distance: int, progress: bool,
    sample_size: int | None = None,
    chrom_pattern: str = _DEFAULT_FRAGLEN_CHROM_PATTERN,
) -> dict:
    """Fastest path: parse pairs file using polars (SIMD-accelerated CSV reader).

    When sample_size is set, uses a two-pass approach (row count + stride sample).
    Without sampling (the default), reads explicit batches and accumulates counters,
    so peak memory is O(batch size) regardless of file size.
    """
    import io as _io

    needed = ["chrom1", "pos1", "chrom2", "pos2", "read_len1", "read_len2",
              "pos51", "pos52", "pos31", "pos32"]
    int32_cols = [c for c in needed if c not in ("chrom1", "chrom2")]
    schema_ov = {c: pl.Int32 for c in int32_cols}

    def _read_csv_kwargs(source):
        return dict(
            source=source,
            separator="\t",
            comment_prefix="#",
            has_header=False,
            new_columns=columns,
            schema_overrides=schema_ov,
        )

    def _scan():
        return pl.scan_csv(**_read_csv_kwargs(str(path))).select(needed)

    if sample_size is not None:
        # Two lazy scan_csv passes — polars streams in batches so peak memory is
        # O(batch_size + sample_size), never proportional to the full file size.
        # Pass 1: exact row count (streaming, ~1s regardless of file size).
        if progress:
            _progress("  pairs: counting rows...")
        total_rows = _scan().select(pl.len()).collect(streaming=True).item()
        stride = max(1, total_rows // sample_size)

        # Pass 2: stream the file, keep every stride-th row.
        if progress:
            _progress(f"  pairs: {total_rows:,} reads — reading 1/{stride} sample...")
        df = _scan().gather_every(stride).collect()
        scale = total_rows / max(len(df), 1)

        return _polars_metrics_from_df(df, total_rows, scale, cis_distance, progress,
                                       chrom_pattern)

    return _parse_pairs_polars_batched(path, columns, cis_distance, progress,
                                       chrom_pattern=chrom_pattern)


def _iter_polars_batches(path: Path, columns: list[str], needed: list[str], chunk_bytes: int):
    """Yield polars DataFrames of `needed` columns (ints as Int32) with bounded memory.

    The file is read in fixed-size byte chunks cut at line boundaries and each chunk is parsed
    by polars, so peak memory is O(chunk_bytes) whatever the file size and polars version.
    (Lazy scan_csv + group_by/unpivot is not reliably streamed: polars 1.8 used ~50 GB on a
    105 GB pairs file; read_csv_batched memory-maps/reads ahead: ~22 GB peak on a 30 GB file.)
    Columns are read as strings (no per-chunk dtype inference, so numeric chromosome names
    cannot break later chunks) and the integer columns are cast.
    """
    idx = sorted(columns.index(c) for c in needed)
    names = [columns[i] for i in idx]
    int_cols = [c for c in needed if c not in ("chrom1", "chrom2")]
    casts = [pl.col(c).cast(pl.Int32, strict=True) for c in int_cols]

    def _parse(data: bytes):
        df = pl.read_csv(data, separator="\t", has_header=False, columns=idx,
                         infer_schema_length=0, comment_prefix="#")
        df.columns = names
        return df.with_columns(casts)

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:
        rest = b""
        in_header = True
        while True:
            block = fh.read(chunk_bytes)
            data = rest + block
            if not block:
                rest = b""
            else:
                cut = data.rfind(b"\n")
                if cut < 0:
                    rest = data
                    continue
                data, rest = data[:cut + 1], data[cut + 1:]
            if in_header:
                # drop leading header/comment lines (pairs headers precede all data)
                while data.startswith(b"#"):
                    nl = data.find(b"\n")
                    data = b"" if nl < 0 else data[nl + 1:]
                if data:
                    in_header = False
            if data.strip():
                yield _parse(data)
            if not block:
                break


def _value_counts_into(counter: Counter, series) -> None:
    vc = series.value_counts()
    for value, count in zip(vc.to_series(0).to_list(), vc.to_series(1).to_list()):
        counter[int(value)] += int(count)


def _parse_pairs_polars_batched(
    path: Path, columns: list[str], cis_distance: int, progress: bool,
    chunk_bytes: int = 128 << 20,
    chrom_pattern: str = _DEFAULT_FRAGLEN_CHROM_PATTERN,
) -> dict:
    """Exact (no sampling) polars path with memory bounded by chunk_bytes, not file size."""
    needed = ["chrom1", "pos1", "chrom2", "pos2", "read_len1", "read_len2",
              "pos51", "pos52", "pos31", "pos32"]
    total_rows = 0
    cis_lr = 0
    fragment_lengths: Counter[int] = Counter()
    chroms: set[str] = set()
    chrm = 0
    read_lengths: Counter[int] = Counter()
    prev_milestone = 0
    for df in _iter_polars_batches(path, columns, needed, chunk_bytes):
        total_rows += df.height
        cis_lr += int(df.select(
            ((pl.col("chrom1") == pl.col("chrom2")) &
             ((pl.col("pos2") - pl.col("pos1")).abs() >= cis_distance)).sum()
        ).item())
        chrm += _polars_fraglen_into(fragment_lengths, chroms, df, chrom_pattern)
        rl = pl.concat([
            df["read_len1"].filter(df["read_len1"] > 0).rename("v"),
            df["read_len2"].filter(df["read_len2"] > 0).rename("v"),
        ])
        _value_counts_into(read_lengths, rl)
        if progress:
            milestone = (total_rows // (100 * _PROGRESS_INTERVAL)) * (100 * _PROGRESS_INTERVAL)
            if milestone > prev_milestone:
                _progress(f"  pairs: {total_rows:,} reads")
                prev_milestone = milestone

    if progress:
        _progress(f"  pairs: {total_rows:,} reads (100%)", end="\n")

    return {
        "non_dup_reads": total_rows,
        "cis_long_range_pairs": cis_lr,
        "read_length": infer_consensus_read_length(read_lengths),
        "read_length_distribution": dict(sorted(read_lengths.items())),
        "fragment_length_distribution": dict(sorted(fragment_lengths.items())),
        "fragment_length_chromosomes": sorted(chroms),
        "fragment_count": 2 * total_rows,
        "chrM_fragment_count": chrm,
    }


def _chrom_matcher(chrom_pattern: str):
    """Return a cached chromosome-name -> bool function for `chrom_pattern` (regex search)."""
    regex = re.compile(chrom_pattern)
    cache: dict[str, bool] = {}

    def keep(chrom: str) -> bool:
        hit = cache.get(chrom)
        if hit is None:
            hit = cache[chrom] = regex.search(chrom) is not None
        return hit

    keep.matched = lambda: sorted(c for c, hit in cache.items() if hit)
    return keep


def _parse_pairs_pandas(
    path: Path, columns: list[str], cis_distance: int, progress: bool,
    sample_size: int | None = None,
    chrom_pattern: str = _DEFAULT_FRAGLEN_CHROM_PATTERN,
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
    keep = _chrom_matcher(chrom_pattern)
    chrm = 0

    chunks = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=columns,
        usecols=needed,
        chunksize=1_000_000,
        dtype={
            "chrom1": "str", "chrom2": "str",
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

        k1 = chunk.chrom1.map(keep).astype(bool)
        k2 = chunk.chrom2.map(keep).astype(bool)
        frag1 = ((chunk.pos31 - chunk.pos51).abs() + 1)[k1]
        frag2 = ((chunk.pos32 - chunk.pos52).abs() + 1)[k2]
        for length, count in pd.concat([frag1, frag2]).value_counts().items():
            fragment_lengths[int(length)] += int(count)
        chrm += int(chunk.chrom1.isin(_CHRM_NAMES).sum() + chunk.chrom2.isin(_CHRM_NAMES).sum())

    non_dup_reads = total_rows if total_rows is not None else base_row
    scale = non_dup_reads / sampled_rows if sampled_rows < non_dup_reads else 1.0

    if progress:
        _progress(f"  pairs: {non_dup_reads:,} reads (100%)", end="\n")

    if scale > 1.0:
        cis_long_range_pairs = int(round(cis_long_range_pairs * scale))
        fragment_lengths = Counter({k: int(round(v * scale)) for k, v in fragment_lengths.items()})
        chrm = int(round(chrm * scale))

    return {
        "non_dup_reads": non_dup_reads,
        "cis_long_range_pairs": cis_long_range_pairs,
        "read_length": infer_consensus_read_length(read_lengths),
        "read_length_distribution": dict(sorted(read_lengths.items())),
        "fragment_length_distribution": dict(sorted(fragment_lengths.items())),
        "fragment_length_chromosomes": keep.matched(),
        "fragment_count": 2 * non_dup_reads,
        "chrM_fragment_count": chrm,
    }


def _parse_pairs_python(
    path: Path, columns: list[str], cis_distance: int, progress: bool,
    sample_size: int | None = None,
    chrom_pattern: str = _DEFAULT_FRAGLEN_CHROM_PATTERN,
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
    keep = _chrom_matcher(chrom_pattern)
    chrm = 0
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
            if keep(chrom1):
                fragment_lengths[abs(pos31 - pos51) + 1] += 1
            elif chrom1 in _CHRM_NAMES:
                chrm += 1
            if keep(chrom2):
                fragment_lengths[abs(pos32 - pos52) + 1] += 1
            elif chrom2 in _CHRM_NAMES:
                chrm += 1

    non_dup_reads = total_rows if total_rows is not None else data_line
    sampled = data_line // stride  # number of rows actually processed
    scale = non_dup_reads / sampled if sampled < non_dup_reads else 1.0

    if progress:
        _progress(f"  pairs: {non_dup_reads:,} reads (100%)", end="\n")

    if scale > 1.0:
        cis_long_range_pairs = int(round(cis_long_range_pairs * scale))
        fragment_lengths = Counter({k: int(round(v * scale)) for k, v in fragment_lengths.items()})
        chrm = int(round(chrm * scale))

    return {
        "non_dup_reads": non_dup_reads,
        "cis_long_range_pairs": cis_long_range_pairs,
        "read_length": infer_consensus_read_length(read_lengths),
        "read_length_distribution": dict(sorted(read_lengths.items())),
        "fragment_length_distribution": dict(sorted(fragment_lengths.items())),
        "fragment_length_chromosomes": keep.matched(),
        "fragment_count": 2 * non_dup_reads,
        "chrM_fragment_count": chrm,
    }


BACKENDS = ("auto", "polars", "pandas", "python")


def resolve_backend(backend: str = "auto") -> str:
    """Map a requested backend to an available one ('auto': polars > pandas > python)."""
    if backend not in BACKENDS:
        raise ValueError(f"unknown backend {backend!r}; choose from {', '.join(BACKENDS)}")
    if backend == "auto":
        return "polars" if _POLARS_AVAILABLE else ("pandas" if _PANDAS_AVAILABLE else "python")
    if backend == "polars" and not _POLARS_AVAILABLE:
        raise ValueError("backend 'polars' requested but polars is not installed")
    if backend == "pandas" and not _PANDAS_AVAILABLE:
        raise ValueError("backend 'pandas' requested but pandas is not installed")
    return backend


def parse_pairs_file(
    path: Path,
    cis_distance: int = 10_000,
    progress: bool = True,
    sample_size: int | None = None,
    backend: str = "auto",
    chrom_pattern: str = _DEFAULT_FRAGLEN_CHROM_PATTERN,
) -> dict:
    """Parse a pairs/pairsam file and compute QC metrics.

    Uses the fastest available library:
      polars  (install: pip install polars)  — ~30× vs pure Python
      pandas  (install: pip install pandas)  — ~4× vs pure Python
      pure Python                            — always available, no extra deps

    When sample_size is set, non_dup_reads is always exact (full file count) while
    rate/distribution metrics are estimated from a stride sample of sample_size reads.
    Set sample_size=0 to disable sampling and read every row.

    fragment_length_distribution counts only read ends whose own chromosome matches
    `chrom_pattern` (regex search; default: autosomes + X). chrM read ends are counted
    separately as chrM_fragment_count.

    All three backends give identical results and, in exact mode, use memory bounded by their
    batch size rather than the file size. `backend` selects one explicitly ('auto' = fastest
    available).
    """
    if sample_size == 0:
        sample_size = None
    columns, _ = _parse_pairs_header(path)
    backend = resolve_backend(backend)
    print(f"  pairs: backend={backend}", file=sys.stderr)
    if backend == "polars":
        return _parse_pairs_polars(path, columns, cis_distance, progress, sample_size,
                                   chrom_pattern)
    if backend == "pandas":
        return _parse_pairs_pandas(path, columns, cis_distance, progress, sample_size,
                                   chrom_pattern)
    return _parse_pairs_python(path, columns, cis_distance, progress, sample_size,
                               chrom_pattern)


def is_unique_alignment(record, unique_mapq_min: int) -> bool:
    """Determine unique mapping directly from BAM fields."""
    if record.is_unmapped:
        return False
    if record.has_tag("NH"):
        return record.get_tag("NH") == 1
    return record.mapping_quality >= unique_mapq_min


def parse_bam_file(path: Path, unique_mapq_min: int = 20, progress: bool = True) -> dict:
    """Compute read-based metrics from BAM using pysam.

    Reads only the first ``_BAM_SAMPLE_READS`` primary reads (exits early).
    This is fast even for very large BAM files.

    ``total_reads`` and ``unique_reads`` are not computed here because exact
    counts require scanning the entire BAM — use a pairtools stats file via
    ``--stats`` for exact read totals.  ``unique_fraction`` is the fraction of
    sampled primary alignments that are unique (NH==1 when present, else
    MAPQ >= ``unique_mapq_min``).  For a coordinate-sorted BAM the sample comes
    from the start of the first reference sequence(s), so treat it as an estimate.
    """
    try:
        import pysam
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "pysam is required to read BAM files. "
            "Install it with conda or pip, or omit --bam."
        ) from exc

    sampled_total = 0
    sampled_unique = 0
    read_lengths: Counter[int] = Counter()

    if progress:
        _progress(f"  bam:   sampling first {_BAM_SAMPLE_READS:,} reads...")
    with pysam.AlignmentFile(str(path), "rb") as bam_file:
        for record in bam_file.fetch(until_eof=True):
            if record.is_secondary or record.is_supplementary:
                continue
            sampled_total += 1
            if record.query_length:
                read_lengths[record.query_length] += 1
            if is_unique_alignment(record, unique_mapq_min):
                sampled_unique += 1
            if sampled_total >= _BAM_SAMPLE_READS:
                break

    unique_fraction = sampled_unique / sampled_total if sampled_total else 0.0

    if progress:
        _progress(
            f"  bam:   unique fraction {unique_fraction:.3f} "
            f"(from {sampled_total:,}-read sample)",
            end="\n",
        )

    return {
        "total_reads": None,          # requires full scan or stats file; not computed here
        "unique_reads": None,         # populated by build_summary when total_reads is known
        "unique_fraction": unique_fraction,
        "read_length": infer_consensus_read_length(read_lengths),
        "read_length_distribution": dict(sorted(read_lengths.items())),
        "unique_mapq_min": unique_mapq_min,
        "unique_fraction_sample_size": sampled_total,
    }


def parse_stats_file(path: Path) -> dict:
    """Parse a pairtools stats file and return key QC metrics."""
    data: dict[str, int] = {}
    with path.open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                try:
                    data[parts[0]] = int(parts[1])
                except ValueError:
                    pass  # skip non-integer rows (chrom_freq, dist_freq, etc.)

    total = data.get("total")
    total_mapped = data.get("total_mapped")
    total_dups = data.get("total_dups")
    total_nodups = data.get("total_nodups")

    return {
        "total_reads": total,
        "total_mapped": total_mapped,
        "total_dups": total_dups,
        "total_nodups": total_nodups,
        "mapping_rate": total_mapped / total if total else None,
        "duplication_rate": (
            total_dups / total_mapped if (total_dups is not None and total_mapped) else None
        ),
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


def discover_stats(directory: Path, sample_ids: list[str]) -> list[Path | None]:
    """Find ``<sample_id>.stats.txt`` pairtools stats files under ``directory``.

    Returns one entry per sample id (``None`` when no stats file is found).
    Raises if a sample has more than one candidate stats file.
    """
    found: dict[str, list[Path]] = {}
    for path in sorted(directory.rglob("*.stats.txt")):
        if path.is_file():
            found.setdefault(path.name[: -len(".stats.txt")], []).append(path)
    stats_paths: list[Path | None] = []
    for sample_id in sample_ids:
        candidates = found.get(sample_id, [])
        if len(candidates) > 1:
            raise ValueError(
                f"Multiple stats files found for sample {sample_id}: "
                + ", ".join(str(c) for c in candidates)
            )
        stats_paths.append(candidates[0] if candidates else None)
    return stats_paths


def build_summary(
    sample_id: str,
    pairs_metrics: dict,
    bam_metrics: dict | None,
    stats_metrics: dict | None,
    cis_distance: int,
    pairs_path: Path,
    bam_path: Path | None,
    stats_path: Path | None,
    sample_size: int | None = None,
    chrom_pattern: str = _DEFAULT_FRAGLEN_CHROM_PATTERN,
) -> dict:
    """Assemble the final QC summary document."""
    read_length = pairs_metrics["read_length"]
    read_length_source = "pairs"
    if bam_metrics and bam_metrics.get("read_length") is not None:
        read_length = bam_metrics["read_length"]
        read_length_source = "bam"

    # Prefer stats file for all mapping/duplication metrics; fall back to BAM sample.
    if stats_metrics:
        total_reads = stats_metrics["total_reads"]
        total_mapped = stats_metrics["total_mapped"]
        total_dups = stats_metrics["total_dups"]
        mapping_rate = stats_metrics["mapping_rate"]
        duplication_rate = stats_metrics["duplication_rate"]
        total_reads_source = f"stats({stats_path.name})"
    else:
        total_reads = None
        total_mapped = None
        total_dups = None
        mapping_rate = None
        duplication_rate = None
        total_reads_source = "unavailable — provide --stats for exact counts"

    sources: dict = {
        "read_length": read_length_source,
        "total_reads": total_reads_source,
        "total_mapped": total_reads_source,
        "total_dups": total_reads_source,
        "mapping_rate": total_reads_source,
        "duplication_rate": total_reads_source,
        "non_dup_reads": "pairs(exact line count)",
        "cis_long_range_pairs": f"pairs(same chromosome and distance >= {cis_distance})",
        "fragment_length_distribution": (
            "pairs(abs(pos31-pos51)+1, abs(pos32-pos52)+1; each read end kept if its own "
            f"chromosome matches {chrom_pattern!r})"
        ),
        "chrM_fraction": f"pairs(read ends on {'/'.join(_CHRM_NAMES)} / all read ends)",
        "fraction_fragments_le80bp": "fragment_length_distribution(length <= 80 / total)",
        "fraction_fragments_ge_read_length": (
            "fragment_length_distribution(length >= read_length / total); aligned spans cannot "
            "exceed read length except via gapped alignments, so this is the read-length pile-up"
        ),
    }
    unique_fraction = bam_metrics.get("unique_fraction") if bam_metrics else None
    unique_reads = bam_metrics.get("unique_reads") if bam_metrics else None
    if bam_metrics and (unique_fraction is not None or unique_reads is not None):
        rule = f"NH==1 when present, else MAPQ>={bam_metrics.get('unique_mapq_min')}"
        n_sampled = bam_metrics.get("unique_fraction_sample_size")
        scope = f"first {n_sampled:,} primary reads" if n_sampled else "primary reads"
        if unique_fraction is not None:
            sources["unique_fraction"] = f"bam({scope}; {rule})"
        if unique_reads is not None:
            sources["unique_reads"] = f"bam({scope}; {rule})"
    fraglen = pairs_metrics["fragment_length_distribution"]
    nuclear_fragments = sum(fraglen.values())

    def _fraglen_fraction(pred) -> float | None:
        if not nuclear_fragments:
            return None
        return sum(v for k, v in fraglen.items() if pred(k)) / nuclear_fragments

    fragment_count = pairs_metrics.get("fragment_count")
    chrm_count = pairs_metrics.get("chrM_fragment_count")
    chrm_fraction = chrm_count / fragment_count if fragment_count and chrm_count is not None else None

    if sample_size:
        sources["sampling"] = (
            f"stride sample of {sample_size:,} from {pairs_metrics['non_dup_reads']:,} reads"
        )

    precision: dict = {
        "read_length": "exact",
        "total_reads": "exact" if stats_metrics else None,
        "total_mapped": "exact" if stats_metrics else None,
        "total_dups": "exact" if stats_metrics else None,
        "mapping_rate": "exact" if stats_metrics else None,
        "duplication_rate": "exact" if stats_metrics else None,
        "non_dup_reads": "exact",
        "cis_long_range_pairs": "exact",
    }
    if unique_fraction is not None:
        precision["unique_fraction"] = (
            "estimate" if bam_metrics.get("unique_fraction_sample_size") else "exact"
        )

    return {
        "sample_id": sample_id,
        "microc_qc_version": __version__,
        "inputs": {
            "pairs": str(pairs_path),
            "bam": str(bam_path) if bam_path else None,
            "stats": str(stats_path) if stats_path else None,
        },
        "metrics": {
            "read_length": read_length,
            "total_reads": total_reads,
            "total_mapped": total_mapped,
            "total_dups": total_dups,
            "mapping_rate": mapping_rate,
            "duplication_rate": duplication_rate,
            "unique_fraction": unique_fraction,
            "unique_fraction_sample_size": (
                bam_metrics.get("unique_fraction_sample_size") if bam_metrics else None
            ),
            "unique_reads": unique_reads,
            "non_dup_reads": pairs_metrics["non_dup_reads"],
            "cis_long_range_pairs": pairs_metrics["cis_long_range_pairs"],
            "nuclear_fragment_count": nuclear_fragments,
            "chrM_fragment_count": chrm_count,
            "chrM_fraction": chrm_fraction,
            "fraction_fragments_le80bp": _fraglen_fraction(lambda k: k <= 80),
            "fraction_fragments_ge_read_length": (
                _fraglen_fraction(lambda k: k >= read_length) if read_length else None
            ),
        },
        "precision": precision,
        "fragment_length_distribution": fraglen,
        "fragment_length_chromosomes": pairs_metrics.get("fragment_length_chromosomes"),
        "fragment_length_chrom_pattern": chrom_pattern,
        "sources": sources,
    }


def build_batch_summary(
    sample_ids: list[str],
    pairs_paths: list[Path],
    bam_paths: list[Path | None],
    stats_paths: list[Path | None],
    cis_distance: int,
    unique_mapq_min: int,
    sample_size: int | None = None,
    backend: str = "auto",
    chrom_pattern: str = _DEFAULT_FRAGLEN_CHROM_PATTERN,
) -> dict:
    """Assemble per-sample summaries for a batch run."""
    n = len(sample_ids)
    print(f"Processing {n} sample{'s' if n != 1 else ''}...", file=sys.stderr)
    samples = []
    for i, (sample_id, pairs_path, bam_path, stats_path) in enumerate(
        zip(sample_ids, pairs_paths, bam_paths, stats_paths), start=1
    ):
        print(f"[{i}/{n}] {sample_id}", file=sys.stderr)

        # Parse stats file immediately (tiny file, no I/O cost).
        stats_metrics = parse_stats_file(stats_path) if stats_path else None
        if stats_metrics:
            print(
                f"  stats: total={stats_metrics['total_reads']:,}  "
                f"mapped={stats_metrics['total_mapped']:,}  "
                f"dups={stats_metrics['total_dups']:,}  "
                f"mapping_rate={stats_metrics['mapping_rate']:.3f}  "
                f"dup_rate={stats_metrics['duplication_rate']:.3f}",
                file=sys.stderr,
            )

        # Run pairs (and optionally BAM) processing.
        if bam_path:
            with ThreadPoolExecutor(max_workers=2) as pool:
                pairs_future = pool.submit(
                    parse_pairs_file, pairs_path, cis_distance=cis_distance, sample_size=sample_size,
                    backend=backend, chrom_pattern=chrom_pattern,
                )
                bam_future = pool.submit(
                    parse_bam_file, bam_path, unique_mapq_min=unique_mapq_min, progress=False
                )
                pairs_metrics = pairs_future.result()
                bam_metrics = bam_future.result()
        else:
            pairs_metrics = parse_pairs_file(
                pairs_path, cis_distance=cis_distance, sample_size=sample_size, backend=backend,
                chrom_pattern=chrom_pattern,
            )
            bam_metrics = None

        samples.append(
            build_summary(
                sample_id=sample_id,
                pairs_metrics=pairs_metrics,
                bam_metrics=bam_metrics,
                stats_metrics=stats_metrics,
                cis_distance=cis_distance,
                pairs_path=pairs_path,
                bam_path=bam_path,
                stats_path=stats_path,
                sample_size=sample_size,
                chrom_pattern=chrom_pattern,
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
    parser.add_argument("--version", action="version", version=f"microc-qc {__version__}")
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
        help="BAM files matching --pairs (optional; used for read-length cross-check).",
    )
    parser.add_argument(
        "--stats",
        type=Path,
        nargs="+",
        help=(
            "Pairtools stats files matching --pairs order. "
            "Provides exact total_reads, mapping_rate, and duplication_rate."
        ),
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
            "Stride-sample N read pairs from each pairs file for approximate rate and "
            "distribution metrics. Default: 0 (exact — streaming aggregations make "
            "exact as fast and memory-efficient as sampling for any file size)."
        ),
    )
    parser.add_argument(
        "--backend",
        choices=BACKENDS,
        default="auto",
        help=("Pairs parsing backend (default: auto = polars if installed, else pandas, else pure "
              "Python). All give identical output with memory bounded by batch size."),
    )
    parser.add_argument(
        "--fraglen-chrom-pattern",
        default=_DEFAULT_FRAGLEN_CHROM_PATTERN,
        metavar="REGEX",
        help=("Chromosomes whose read ends enter fragment_length_distribution (regex search). "
              "Default: autosomes + X, optional 'chr' prefix and ':start-end' subregion suffix; "
              "excludes chrM, chrY and scaffolds."),
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
        stats_paths: list[Path | None] = discover_stats(args.dir, sample_ids)
        if args.sample_ids is not None:
            raise ValueError("--sample-ids cannot be used with --dir")
    else:
        pairs_paths = args.pairs
        bam_paths = args.bams or [None] * len(pairs_paths)
        stats_paths = args.stats or [None] * len(pairs_paths)

        if len(bam_paths) != len(pairs_paths):
            raise ValueError("--bams must have the same number of entries as --pairs")
        if len(stats_paths) != len(pairs_paths):
            raise ValueError("--stats must have the same number of entries as --pairs")
        if args.sample_ids is not None and len(args.sample_ids) != len(pairs_paths):
            raise ValueError("--sample-ids must match the number of pairs files")

        sample_ids = args.sample_ids
        if sample_ids is None:
            sample_ids = [
                infer_sample_id(pairs_path=p, bam_path=b)
                for p, b in zip(pairs_paths, bam_paths)
            ]

    summary = build_batch_summary(
        sample_ids=sample_ids,
        pairs_paths=pairs_paths,
        bam_paths=bam_paths,
        stats_paths=stats_paths,
        cis_distance=args.cis_distance,
        unique_mapq_min=args.unique_mapq_min,
        sample_size=args.sample_size if args.sample_size != 0 else None,
        backend=args.backend,
        chrom_pattern=args.fraglen_chrom_pattern,
    )

    write_sample_summaries(summary["samples"], args.out)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
