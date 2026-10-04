import importlib.util
import json
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO_ROOT / "microc-qc.py"


spec = importlib.util.spec_from_file_location("microc_qc", MODULE_PATH)
microc_qc = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(microc_qc)

FIXTURE_DIR = REPO_ROOT / "test-output"
needs_fixtures = pytest.mark.skipif(
    not (FIXTURE_DIR / "small-rcmc.mapped.pairs").exists(),
    reason="test-output/small-rcmc fixtures not present",
)


@needs_fixtures
def test_parse_pairs_file_small_rcmc():
    pairs_path = REPO_ROOT / "test-output" / "small-rcmc.mapped.pairs"
    metrics = microc_qc.parse_pairs_file(pairs_path)

    assert metrics["read_length"] == 101
    assert metrics["non_dup_reads"] == 22110
    assert metrics["cis_long_range_pairs"] == 8379
    assert metrics["fragment_count"] == 44220
    assert metrics["fragment_length_distribution"][101] == 27300


@needs_fixtures
def test_build_summary_single_sample_shape():
    pairs_path = REPO_ROOT / "test-output" / "small-rcmc.mapped.pairs"
    bam_path = REPO_ROOT / "test-output" / "small-rcmc.bam"

    summary = microc_qc.build_summary(
        sample_id="small-rcmc",
        pairs_metrics=microc_qc.parse_pairs_file(pairs_path),
        bam_metrics={
            "unique_reads": 43776,
            "read_length": 101,
            "total_reads": 44220,
            "unique_mapq_min": 20,
        },
        stats_metrics=None,
        cis_distance=10_000,
        pairs_path=pairs_path,
        bam_path=bam_path,
        stats_path=None,
    )

    assert summary["sample_id"] == "small-rcmc"
    assert summary["metrics"]["read_length"] == 101
    assert summary["metrics"]["total_reads"] is None  # exact totals now come from --stats
    assert summary["metrics"]["unique_reads"] == 43776
    assert summary["metrics"]["non_dup_reads"] == 22110
    assert summary["metrics"]["cis_long_range_pairs"] == 8379
    assert "pairtools_stats" not in summary
    assert "MAPQ>=20" in summary["sources"]["unique_reads"]


def test_build_batch_summary_wraps_samples():
    pairs_path = REPO_ROOT / "test-output" / "small-rcmc.mapped.pairs"
    bam_path = REPO_ROOT / "test-output" / "small-rcmc.bam"

    original_parse_pairs = microc_qc.parse_pairs_file
    original_parse_bam = microc_qc.parse_bam_file

    try:
        microc_qc.parse_pairs_file = lambda path, cis_distance=10_000, sample_size=None, backend="auto": {
            "non_dup_reads": 10,
            "cis_long_range_pairs": 4,
            "read_length": 101,
            "fragment_length_distribution": {101: 20},
            "fragment_count": 20,
        }
        microc_qc.parse_bam_file = lambda path, unique_mapq_min=20, progress=True: {
            "total_reads": 20,
            "unique_reads": 18,
            "read_length": 101,
            "read_length_distribution": {101: 20},
            "unique_mapq_min": unique_mapq_min,
        }

        batch = microc_qc.build_batch_summary(
            sample_ids=["sample-a", "sample-b"],
            pairs_paths=[pairs_path, pairs_path],
            bam_paths=[bam_path, bam_path],
            stats_paths=[None, None],
            cis_distance=10_000,
            unique_mapq_min=20,
        )
    finally:
        microc_qc.parse_pairs_file = original_parse_pairs
        microc_qc.parse_bam_file = original_parse_bam

    assert list(batch) == ["samples"]
    assert len(batch["samples"]) == 2
    assert batch["samples"][0]["sample_id"] == "sample-a"
    assert batch["samples"][1]["sample_id"] == "sample-b"


def test_write_sample_summaries_single_file(tmp_path):
    sample = {
        "sample_id": "example",
        "metrics": {"total_reads": 10},
        "fragment_length_distribution": {},
        "inputs": {},
        "sources": {},
    }
    output_path = tmp_path / "example.json"

    microc_qc.write_sample_summaries([sample], output_path)

    written = json.loads(output_path.read_text())
    assert written["sample_id"] == "example"


def test_write_sample_summaries_defaults_to_cwd(tmp_path, monkeypatch):
    sample = {
        "sample_id": "example",
        "metrics": {"total_reads": 10},
        "fragment_length_distribution": {},
        "inputs": {},
        "sources": {},
    }
    monkeypatch.chdir(tmp_path)

    microc_qc.write_sample_summaries([sample], None)

    written = json.loads((tmp_path / "example.qc.json").read_text())
    assert written["sample_id"] == "example"


def test_write_sample_summaries_directory(tmp_path):
    samples = [
        {
            "sample_id": "sample-a",
            "metrics": {"total_reads": 10},
            "fragment_length_distribution": {},
            "inputs": {},
            "sources": {},
        },
        {
            "sample_id": "sample-b",
            "metrics": {"total_reads": 20},
            "fragment_length_distribution": {},
            "inputs": {},
            "sources": {},
        },
    ]

    microc_qc.write_sample_summaries(samples, tmp_path / "summaries")

    assert (tmp_path / "summaries" / "sample-a.qc.json").exists()
    assert (tmp_path / "summaries" / "sample-b.qc.json").exists()


@needs_fixtures
def test_discover_samples_finds_matching_pairs_and_bams():
    sample_ids, pairs_paths, bam_paths = microc_qc.discover_samples(REPO_ROOT / "test-output")

    assert sample_ids == ["small-rcmc"]
    assert pairs_paths[0].name == "small-rcmc.mapped.pairs"
    assert bam_paths[0].name == "small-rcmc.bam"


def test_discover_samples_searches_subdirectories(tmp_path):
    nested = tmp_path / "nested" / "deeper"
    nested.mkdir(parents=True)
    bam_path = nested / "example.bam"
    pairs_path = nested / "example.mapped.pairs"
    bam_path.write_text("")
    pairs_path.write_text("")

    sample_ids, pairs_paths, bam_paths = microc_qc.discover_samples(tmp_path)

    assert sample_ids == ["example"]
    assert pairs_paths == [pairs_path]
    assert bam_paths == [bam_path]


@needs_fixtures
def test_parse_bam_file_counts_reads_not_pairs():
    bam_path = REPO_ROOT / "test-output" / "small-rcmc.bam"

    try:
        metrics = microc_qc.parse_bam_file(bam_path)
    except ModuleNotFoundError:
        return

    assert metrics["read_length"] == 101
    assert metrics["total_reads"] is None
    assert metrics["unique_fraction_sample_size"] == 44220  # fixture is smaller than the sample cap
    assert metrics["unique_fraction"] == pytest.approx(43776 / 44220)


@needs_fixtures
def test_parse_pairs_file_python_fallback_matches():
    """Pure-Python path must produce identical output to the fast paths."""
    pairs_path = REPO_ROOT / "test-output" / "small-rcmc.mapped.pairs"
    columns, _ = microc_qc._parse_pairs_header(pairs_path)
    python_metrics = microc_qc._parse_pairs_python(pairs_path, columns, 10_000, False)

    assert python_metrics["read_length"] == 101
    assert python_metrics["non_dup_reads"] == 22110
    assert python_metrics["cis_long_range_pairs"] == 8379
    assert python_metrics["fragment_count"] == 44220
    assert python_metrics["fragment_length_distribution"][101] == 27300


@needs_fixtures
def test_parse_pairs_file_sampling_exact_count_and_close_rate():
    """Stride sampling must return exact non_dup_reads and a close cis-LR rate."""
    pairs_path = REPO_ROOT / "test-output" / "small-rcmc.mapped.pairs"

    # small-rcmc has 22110 rows; sample_size=1000 → stride≈22, ~1005 rows sampled
    metrics = microc_qc.parse_pairs_file(pairs_path, sample_size=1_000, progress=False)

    # Exact count must be preserved regardless of sampling
    assert metrics["non_dup_reads"] == 22110
    assert metrics["fragment_count"] == 2 * 22110

    # Rate metrics may differ from the full-file value, but should be within 20%
    full_cis_lr = 8379
    assert abs(metrics["cis_long_range_pairs"] - full_cis_lr) / full_cis_lr < 0.20


@needs_fixtures
def test_count_data_lines():
    """_count_data_lines must match the actual row count of the pairs file."""
    pairs_path = REPO_ROOT / "test-output" / "small-rcmc.mapped.pairs"
    assert microc_qc._count_data_lines(pairs_path) == 22110


def test_parse_stats_file(tmp_path):
    stats = tmp_path / "s.stats.txt"
    stats.write_text("total\t1000\ntotal_mapped\t800\ntotal_dups\t80\ntotal_nodups\t720\n"
                     "chrom_freq/chr1/chr1\t5\nsummary/frac_cis\t0.5\n")
    m = microc_qc.parse_stats_file(stats)
    assert m["total_reads"] == 1000 and m["total_mapped"] == 800
    assert m["mapping_rate"] == pytest.approx(0.8)
    assert m["duplication_rate"] == pytest.approx(0.1)


def test_discover_stats_matches_sample_ids(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "a" / "s1.stats.txt").write_text("total\t1\n")
    assert microc_qc.discover_stats(tmp_path, ["s1", "s2"]) == [tmp_path / "a" / "s1.stats.txt", None]


def test_build_summary_with_stats_and_bam_sample():
    summary = microc_qc.build_summary(
        sample_id="x",
        pairs_metrics={"non_dup_reads": 10, "cis_long_range_pairs": 4, "read_length": 150,
                       "fragment_length_distribution": {150: 20}, "fragment_count": 20},
        bam_metrics={"total_reads": None, "unique_reads": None, "unique_fraction": 0.9,
                     "read_length": 150, "unique_mapq_min": 20, "unique_fraction_sample_size": 500000},
        stats_metrics={"total_reads": 100, "total_mapped": 80, "total_dups": 8, "total_nodups": 72,
                       "mapping_rate": 0.8, "duplication_rate": 0.1},
        cis_distance=10_000, pairs_path=Path("x.pairs"), bam_path=Path("x.bam"),
        stats_path=Path("x.stats.txt"),
    )
    assert summary["metrics"]["total_reads"] == 100
    assert summary["metrics"]["unique_fraction"] == 0.9
    assert summary["metrics"]["unique_reads"] is None
    assert summary["precision"]["unique_fraction"] == "estimate"
    assert "first 500,000 primary reads" in summary["sources"]["unique_fraction"]


def _write_synthetic_pairs(path, n=25_000, seed=7):
    """Small .pairs file with the columns microc-qc needs (no external fixtures required)."""
    import random
    rng = random.Random(seed)
    chroms = ["chr1", "chr2", "1", "X"]  # include a numeric-looking chromosome name
    cols = ["readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2", "pair_type",
            "read_len1", "read_len2", "pos51", "pos52", "pos31", "pos32"]
    with open(path, "w") as fh:
        fh.write("## pairs format v1.0\n#columns: " + " ".join(cols) + "\n")
        for i in range(n):
            c1 = rng.choice(chroms)
            c2 = c1 if rng.random() < 0.8 else rng.choice(chroms)
            p1 = rng.randint(1, 5_000_000)
            p2 = p1 + rng.randint(0, 200_000) if c1 == c2 else rng.randint(1, 5_000_000)
            rl1 = rng.choice([0, 50, 150, 150, 150])
            rl2 = rng.choice([150, 150, 66])
            f1 = rng.randint(20, 400)
            f2 = rng.randint(20, 400)
            fh.write("\t".join(map(str, [
                f"r{i}", c1, p1, c2, p2, "+", "-", "UU", rl1, rl2,
                p1, p2, p1 + f1 - 1, p2 - f2 + 1])) + "\n")


@pytest.mark.parametrize("backend", ["polars", "pandas"])
def test_backends_match_python_reference(tmp_path, backend):
    """Exact mode: polars (batched) and pandas must equal the pure-Python reference."""
    if backend == "polars" and not microc_qc._POLARS_AVAILABLE:
        pytest.skip("polars not installed")
    if backend == "pandas" and not microc_qc._PANDAS_AVAILABLE:
        pytest.skip("pandas not installed")
    pairs = tmp_path / "synthetic.pairs"
    _write_synthetic_pairs(pairs)
    columns, _ = microc_qc._parse_pairs_header(pairs)
    reference = microc_qc._parse_pairs_python(pairs, columns, 10_000, False)
    got = microc_qc.parse_pairs_file(pairs, progress=False, backend=backend)
    assert got == reference
    assert got["non_dup_reads"] == 25_000
    assert sum(got["fragment_length_distribution"].values()) == 2 * 25_000


def test_polars_batched_is_independent_of_chunk_size(tmp_path):
    if not microc_qc._POLARS_AVAILABLE:
        pytest.skip("polars not installed")
    pairs = tmp_path / "synthetic.pairs"
    _write_synthetic_pairs(pairs, n=10_000)
    columns, _ = microc_qc._parse_pairs_header(pairs)
    a = microc_qc._parse_pairs_polars_batched(pairs, columns, 10_000, False, chunk_bytes=4093)
    b = microc_qc._parse_pairs_polars_batched(pairs, columns, 10_000, False, chunk_bytes=1 << 30)
    assert a == b


def test_resolve_backend_validation():
    assert microc_qc.resolve_backend("auto") in ("polars", "pandas", "python")
    assert microc_qc.resolve_backend("python") == "python"
    with pytest.raises(ValueError):
        microc_qc.resolve_backend("duckdb")


def test_cli_backend_flag(tmp_path):
    pairs = tmp_path / "synthetic.pairs"
    _write_synthetic_pairs(pairs, n=2_000)
    out = tmp_path / "s.qc.json"
    rc = microc_qc.main(["--pairs", str(pairs), "--sample-ids", "s", "--backend", "python",
                         "--out", str(out)])
    assert rc == 0
    data = json.loads(out.read_text())
    blob = json.dumps(data)
    assert "fragment_length_distribution" in blob


def test_polars_batched_reads_gzip(tmp_path):
    if not microc_qc._POLARS_AVAILABLE:
        pytest.skip("polars not installed")
    import gzip, shutil
    pairs = tmp_path / "synthetic.pairs"
    _write_synthetic_pairs(pairs, n=3_000)
    gz = tmp_path / "synthetic.pairs.gz"
    with open(pairs, "rb") as src, gzip.open(gz, "wb") as dst:
        shutil.copyfileobj(src, dst)
    columns, _ = microc_qc._parse_pairs_header(pairs)
    assert (microc_qc._parse_pairs_polars_batched(gz, columns, 10_000, False, chunk_bytes=2048)
            == microc_qc._parse_pairs_python(pairs, columns, 10_000, False))
