import importlib.util
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO_ROOT / "microc-qc.py"


spec = importlib.util.spec_from_file_location("microc_qc", MODULE_PATH)
microc_qc = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(microc_qc)


def test_parse_pairs_file_small_rcmc():
    pairs_path = REPO_ROOT / "test-output" / "small-rcmc.mapped.pairs"
    metrics = microc_qc.parse_pairs_file(pairs_path)

    assert metrics["read_length"] == 101
    assert metrics["non_dup_reads"] == 22110
    assert metrics["cis_long_range_pairs"] == 8379
    assert metrics["fragment_count"] == 44220
    assert metrics["fragment_length_distribution"][101] == 27300


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
        cis_distance=10_000,
        pairs_path=pairs_path,
        bam_path=bam_path,
    )

    assert summary["sample_id"] == "small-rcmc"
    assert summary["metrics"]["read_length"] == 101
    assert summary["metrics"]["total_reads"] == 44220
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
        microc_qc.parse_pairs_file = lambda path, cis_distance=10_000: {
            "non_dup_reads": 10,
            "cis_long_range_pairs": 4,
            "read_length": 101,
            "fragment_length_distribution": {101: 20},
            "fragment_count": 20,
        }
        microc_qc.parse_bam_file = lambda path, unique_mapq_min=20: {
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


def test_parse_bam_file_counts_reads_not_pairs():
    bam_path = REPO_ROOT / "test-output" / "small-rcmc.bam"

    try:
        metrics = microc_qc.parse_bam_file(bam_path)
    except ModuleNotFoundError:
        return

    assert metrics["read_length"] == 101
    assert metrics["total_reads"] == 44220
    assert metrics["unique_reads"] == 43776
