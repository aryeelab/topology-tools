"""Tests for preprocess-hic (nf/preprocess-hic).

Unit tests of the helper scripts run anywhere with Python 3. The end-to-end test runs the Nextflow pipeline
(local profile) on synthetic data from make_test_data.py and needs the preprocess-hic env's bin on PATH
(nextflow, bwa-mem2, bwa, pairtools, samtools, seqkit, cooler, pairix); it is skipped otherwise.

  PATH=<env>/bin:$PATH JAVA_HOME=<env> python -m pytest -v tests/preprocess-hic
"""
import gzip
import json
import os
import shutil
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
NF = os.path.join(REPO, "nf", "preprocess-hic")
BIN = os.path.join(NF, "bin")
sys.path.insert(0, BIN)
sys.path.insert(0, HERE)

import digest_genome  # noqa: E402
import hic_frag_stats  # noqa: E402
import make_test_data  # noqa: E402


# ---------------------------------------------------------------- digest_genome.py
def test_digest_dpnii_hinfi(tmp_path):
    #          0         1         2         3
    #          0123456789012345678901234567890123
    seq = "AAAGATCAAAAGAATCAAAAgatcAAGNNTCAAA"
    (tmp_path / "g.fa").write_text(">c1 x\n" + seq[:20] + "\n" + seq[20:] + "\n>other\nACGT\n")
    (tmp_path / "s").write_text(f"c1\t{len(seq)}\n")
    out = tmp_path / "f.bed"
    digest_genome.main(["--fasta", str(tmp_path / "g.fa"), "--chrom-sizes", str(tmp_path / "s"),
                        "--enzymes", "DpnII,HinfI", "-o", str(out)])
    # DpnII ^GATC at 3 and 20 (lower case matches); HinfI G^ANTC at 11+1; GANNTC with N in the genome: no site
    assert out.read_text().split("\n")[:-1] == ["c1\t0\t3", "c1\t3\t12", "c1\t12\t20", f"c1\t20\t{len(seq)}"]


def test_digest_matches_independent_cutter(tmp_path):
    import random
    rng = random.Random(1)
    seq = "".join(rng.choice("ACGT") for _ in range(20000))
    (tmp_path / "g.fa").write_text(">c\n" + seq + "\n")
    (tmp_path / "s").write_text("c\t20000\n")
    out = tmp_path / "f.bed"
    digest_genome.main(["--fasta", str(tmp_path / "g.fa"), "--chrom-sizes", str(tmp_path / "s"),
                        "--enzymes", "HinfI,DpnII", "-o", str(out)])
    starts = [int(x.split("\t")[1]) for x in out.read_text().splitlines()]
    assert starts == [0] + make_test_data.cuts(seq)


def test_digest_length_mismatch(tmp_path):
    (tmp_path / "g.fa").write_text(">c\nACGTGATC\n")
    (tmp_path / "s").write_text("c\t9\n")
    with pytest.raises(SystemExit):
        digest_genome.main(["--fasta", str(tmp_path / "g.fa"), "--chrom-sizes", str(tmp_path / "s"), "--enzymes", "DpnII"])


def test_digest_unknown_enzyme(tmp_path):
    with pytest.raises(SystemExit):
        digest_genome.main(["--fasta", "x", "--chrom-sizes", "y", "--enzymes", "Foo"])


# ---------------------------------------------------------------- hic_frag_stats.py
HDR = "## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type rfrag1 rfrag_start1 rfrag_end1 rfrag2 rfrag_start2 rfrag_end2\n"


def row(c1, p1, c2, p2, s1, s2, f1, f2):
    return f"r\t{c1}\t{p1}\t{c2}\t{p2}\t{s1}\t{s2}\tUU\t{f1}\t0\t1\t{f2}\t0\t1\n"


def test_frag_stats_categories(tmp_path):
    lines = [row("c", 10, "c", 20, "+", "-", 5, 5),   # dangling end
             row("c", 10, "c", 20, "-", "+", 5, 5),   # self-circle
             row("c", 10, "c", 20, "+", "+", 5, 5),   # same fragment, same strand
             row("c", 10, "c", 90, "+", "-", 5, 6),   # religation
             row("c", 10, "c", 90, "-", "+", 5, 6),   # adjacent, other
             row("c", 10, "c", 900, "+", "-", 5, 9),  # valid
             row("c", 10, "d", 9, "+", "-", 5, 1),    # trans
             row("!", 0, "c", 9, "-", "-", -1, 1)]    # one side unmapped
    p = tmp_path / "x.pairs"
    p.write_text(HDR + "".join(lines))
    n = hic_frag_stats.count_pairs(open(p))
    assert n == {"total": 8, "mapped": 7, "trans": 1, "cis": 6, "dangling_end": 1, "self_circle": 1,
                 "same_fragment_same_strand": 1, "religation": 1, "adjacent_other": 1, "cis_valid": 1}
    a, b = tmp_path / "a", tmp_path / "b"
    a.write_text("".join(f"{k}\t{v}\n" for k, v in n.items()))
    b.write_text("".join(f"{k}\t{v}\n" for k, v in n.items()))
    assert hic_frag_stats.merge([a, b]) == {k: 2 * v for k, v in n.items()}


def test_frag_stats_needs_restrict_columns(tmp_path):
    p = tmp_path / "x.pairs"
    p.write_text("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n")
    with pytest.raises(SystemExit):
        hic_frag_stats.count_pairs(open(p))


# ---------------------------------------------------------------- end to end
TOOLS = ["nextflow", "bwa-mem2", "bwa", "pairtools", "samtools", "seqkit", "cooler", "pairix", "bgzip"]
need_tools = pytest.mark.skipif(any(shutil.which(t) is None for t in TOOLS),
                                reason="needs the preprocess-hic env on PATH: " + ",".join(TOOLS))


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    d = tmp_path_factory.mktemp("e2e")
    data = d / "data"
    make_test_data.main(str(data))
    fa = str(data / "genome.fa")
    for cmd in (["bwa", "index", fa], ["bwa-mem2", "index", fa]):
        subprocess.run(cmd, check=True, capture_output=True)
    out = d / "out"
    out.mkdir()
    cmd = ["nextflow", "run", os.path.join(NF, "main.nf"), "-profile", "local",
           "--sample_id", "test", "--bwa_index", fa, "--genome_fasta", fa,
           "--chrom_sizes", str(data / "chrom.sizes"),
           "--fq_input", str(data / "fastq" / "*_{1,2}.fastq.gz"),
           "--outdir", str(out), "--tmpdir", str(d / "tmp"), "--chunk_reads", "200",
           "--mcool_resolutions", "1000,5000,20000", "--assembly", "synthetic", "--git_commit", "test"]
    r = subprocess.run(cmd, cwd=str(out), capture_output=True, text=True)
    print(r.stdout[-5000:], r.stderr[-5000:])
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return {"out": out, "data": data, "expected": json.loads((data / "expected.json").read_text())}


def read_pairs(path):
    header, rows = [], []
    with gzip.open(path, "rt") as f:
        for line in f:
            (header if line.startswith("#") else rows).append(line.rstrip("\n"))
    cols = [h for h in header if h.startswith("#columns:")][0].split()[1:]
    return header, [dict(zip(cols, r.split("\t"))) for r in rows]


def classify(p):
    if p["chrom1"] != p["chrom2"]:
        return "trans"
    d = int(p["rfrag2"]) - int(p["rfrag1"])
    o = p["strand1"] + p["strand2"]
    if d == 0:
        return {"+-": "dangling_end", "-+": "self_circle"}.get(o, "same_fragment_same_strand")
    if d == 1:
        return "religation" if o == "+-" else "adjacent_other"
    return "cis_valid"


@need_tools
def test_e2e_pairs(run):
    header, pairs = read_pairs(run["out"] / "test.pairs.gz")
    cats = {}
    for p in pairs:
        name_cat = p["readID"].split("_", 1)[1]
        cats[name_cat] = cats.get(name_cat, 0) + 1
        assert classify(p) == ("cis_valid" if name_cat == "dup" else name_cat), p
        assert int(p["mapq1"]) >= 1 and int(p["mapq2"]) >= 1
    n = run["expected"]["categories"]
    assert cats.get("cis_valid", 0) + cats.get("dup", 0) == n["cis_valid"]   # one of each duplicate kept
    for c in ["trans", "dangling_end", "self_circle", "religation"]:
        assert cats.get(c, 0) == n[c], c
    for c in ["multi", "unmapped", "half_mapped"]:
        assert c not in cats
    parse = [h for h in header if "ID:pairtools_parse" in h][0]
    assert "--min-mapq 1" in parse and "--walks-policy 5unique" in parse and "mapq,pos5,pos3" in parse
    assert any("ID:bwa-mem2" in h and "-SP5M" in h for h in header)
    assert any("ID:pairtools_restrict" in h for h in header)
    cols = [h for h in header if h.startswith("#columns:")][0].split()[1:]
    for c in ["mapq1", "mapq2", "pos51", "pos31", "rfrag1", "rfrag_start1", "rfrag_end2"]:
        assert c in cols
    # pairix index: query the trans block
    r = subprocess.run(["pairix", str(run["out"] / "test.pairs.gz"), "chrA|chrB"], capture_output=True, text=True, check=True)
    assert len(r.stdout.splitlines()) == n["trans"]


@need_tools
def test_e2e_stats(run):
    n = run["expected"]["categories"]
    s = dict(l.split("\t")[:2] for l in (run["out"] / "test.dedup.stats.txt").read_text().splitlines() if "\t" in l)
    assert int(s["total"]) == run["expected"]["read_pairs"]
    assert int(s["total_dups"]) == n["dup"]
    assert int(s["total_unmapped"]) == n["unmapped"]
    nodup = n["cis_valid"] + n["trans"] + n["dangling_end"] + n["self_circle"] + n["religation"]
    assert int(s["total_nodups"]) == nodup
    r = dict(l.split("\t") for l in (run["out"] / "test.restrict_stats.tsv").read_text().splitlines())
    assert {k: int(v) for k, v in r.items()} == {
        "total": nodup, "mapped": nodup, "trans": n["trans"], "cis": nodup - n["trans"],
        "dangling_end": n["dangling_end"], "self_circle": n["self_circle"], "same_fragment_same_strand": 0,
        "religation": n["religation"], "adjacent_other": 0, "cis_valid": n["cis_valid"]}
    qc = json.loads((run["out"] / "test.hic_qc.json").read_text())
    assert qc["bam"]["lossless_check"] == "pass"
    assert qc["read_pairs"] == run["expected"]["read_pairs"]
    assert qc["restriction_fragments"]["frac_dangling_end"] == round(n["dangling_end"] / nodup, 6)


@need_tools
def test_e2e_bam_lossless(run):
    bam = str(run["out"] / "test.bam")
    assert os.path.exists(bam + ".bai")
    names = subprocess.run(f"samtools view -F 0x900 {bam} | cut -f1 | sort | uniq -c", shell=True,
                           capture_output=True, text=True, check=True).stdout.split()
    counts = dict(zip(names[1::2], map(int, names[0::2])))
    assert len(counts) == run["expected"]["read_pairs"]       # every read pair, every pair type
    assert set(counts.values()) == {2}                         # both reads of each
    hdr = subprocess.run(["samtools", "view", "-H", bam], capture_output=True, text=True, check=True).stdout
    assert "SO:coordinate" in hdr
    yt = subprocess.run(f"samtools view {bam} | grep -o 'Yt:Z:[A-Z]*' | sort -u", shell=True,
                        capture_output=True, text=True).stdout.split()
    assert {"Yt:Z:UU", "Yt:Z:NN", "Yt:Z:MU"} <= set(yt) or {"Yt:Z:UU", "Yt:Z:NN", "Yt:Z:UM"} <= set(yt)


@need_tools
def test_e2e_mcool(run):
    import cooler
    mc = str(run["out"] / "test.mapq_30.mcool")
    assert cooler.fileops.list_coolers(mc) == ["/resolutions/1000", "/resolutions/5000", "/resolutions/20000"]
    n = run["expected"]["categories"]
    nodup = n["cis_valid"] + n["trans"] + n["dangling_end"] + n["self_circle"] + n["religation"]
    for res in (1000, 5000, 20000):
        c = cooler.Cooler(f"{mc}::/resolutions/{res}")
        assert c.info["sum"] == nodup                           # all unique reads have MAPQ 60
        assert "weight" in c.bins().columns
    meta = cooler.Cooler(f"{mc}::/resolutions/1000").info["metadata"]
    assert meta["min_mapq"] == 30 and meta["pipeline"] == "preprocess-hic"


@need_tools
def test_e2e_provenance(run):
    prov = json.loads((run["out"] / "test.preprocess-hic.provenance.json").read_text())
    assert prov["success"] is True and prov["pipeline"] == "preprocess-hic"
    assert prov["params"]["bwa_mem_options"] == "-SP5M" and prov["params"]["min_mapq"] == 1
    assert len(prov["inputs"]) == 4
    assert any(v.startswith("pairtools") and v.endswith("1.1.3") for v in prov["tool_versions"])
    assert prov["qc"]["bam"]["lossless_check"] == "pass"
