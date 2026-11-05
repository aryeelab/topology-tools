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
def data(tmp_path_factory):
    data = tmp_path_factory.mktemp("e2e_data")
    make_test_data.main(str(data))
    fa = str(data / "genome.fa")
    for cmd in (["bwa", "index", fa], ["bwa-mem2", "index", fa]):
        subprocess.run(cmd, check=True, capture_output=True)
    return data


# pairtools parse --min-mapq: 1 (default; MAPQ 0 sides are M, their pairs leave the pairs file) and 0 (kept as U).
@pytest.fixture(scope="module", params=[1, 0], ids=["minmapq1", "minmapq0"])
def run(request, data, tmp_path_factory):
    mq = request.param
    d = tmp_path_factory.mktemp(f"e2e_mapq{mq}")
    fa = str(data / "genome.fa")
    out = d / "out"
    out.mkdir()
    cmd = ["nextflow", "run", os.path.join(NF, "main.nf"), "-profile", "local",
           "--sample_id", "test", "--bwa_index", fa, "--genome_fasta", fa,
           "--chrom_sizes", str(data / "chrom.sizes"),
           "--fq_input", str(data / "fastq" / "*_{1,2}.fastq.gz"),
           "--outdir", str(out), "--tmpdir", str(d / "tmp"), "--chunk_reads", "200",
           "--mcool_resolutions", "1000,5000,20000", "--assembly", "synthetic", "--git_commit", "test",
           "--min_mapq", str(mq)]
    r = subprocess.run(cmd, cwd=str(out), capture_output=True, text=True)
    print(r.stdout[-5000:], r.stderr[-5000:])
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    exp = json.loads((data / "expected.json").read_text())
    n = exp["categories"]
    nodup = n["cis_valid"] + n["trans"] + n["dangling_end"] + n["self_circle"] + n["religation"]
    # with --min-mapq 0 the pairs with one side in the duplicated segment are kept (as UU, that side mapq 0)
    return {"out": out, "data": data, "expected": exp, "mq": mq, "nodup_unique": nodup,
            "nodup": nodup + (n["multi"] if mq == 0 else 0)}


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
        if name_cat == "multi":
            assert min(int(p["mapq1"]), int(p["mapq2"])) == 0 and p["pair_type"] == "UU", p
            continue
        assert classify(p) == ("cis_valid" if name_cat == "dup" else name_cat), p
        assert int(p["mapq1"]) >= 1 and int(p["mapq2"]) >= 1
    n = run["expected"]["categories"]
    assert cats.get("cis_valid", 0) + cats.get("dup", 0) == n["cis_valid"]   # one of each duplicate kept
    for c in ["trans", "dangling_end", "self_circle", "religation"]:
        assert cats.get(c, 0) == n[c], c
    assert cats.get("multi", 0) == (n["multi"] if run["mq"] == 0 else 0)
    for c in ["unmapped", "half_mapped"]:
        assert c not in cats
    parse = [h for h in header if "ID:pairtools_parse" in h][0]
    assert f"--min-mapq {run['mq']}" in parse and "--walks-policy 5unique" in parse and "mapq,pos5,pos3" in parse
    assert any("ID:bwa-mem2" in h and "-SP5M" in h for h in header)
    assert any("ID:pairtools_restrict" in h for h in header)
    cols = [h for h in header if h.startswith("#columns:")][0].split()[1:]
    for c in ["mapq1", "mapq2", "pos51", "pos31", "rfrag1", "rfrag_start1", "rfrag_end2"]:
        assert c in cols
    # pairix index: query the trans block
    r = subprocess.run(["pairix", str(run["out"] / "test.pairs.gz"), "chrA|chrB"], capture_output=True, text=True, check=True)
    trans_multi = sum(1 for p in pairs if p["readID"].endswith("_multi") and p["chrom1"] != p["chrom2"])
    assert len(r.stdout.splitlines()) == n["trans"] + trans_multi


@need_tools
def test_e2e_stats(run):
    n = run["expected"]["categories"]
    s = dict(l.split("\t")[:2] for l in (run["out"] / "test.dedup.stats.txt").read_text().splitlines() if "\t" in l)
    assert int(s["total"]) == run["expected"]["read_pairs"]
    assert int(s["total_dups"]) == n["dup"]
    assert int(s["total_unmapped"]) == n["unmapped"]
    nodup = run["nodup"]
    assert int(s["total_nodups"]) == nodup
    types = {k.split("/", 1)[1] for k in s if k.startswith("pair_types/")}
    assert (not any("M" in t for t in types)) if run["mq"] == 0 else ({"MU", "UM"} & types)
    r = {k: int(v) for k, v in (l.split("\t") for l in (run["out"] / "test.restrict_stats.tsv").read_text().splitlines())}
    extra = n["multi"] if run["mq"] == 0 else 0   # multi pairs: valid cis or trans, by bwa's random primary
    assert r["total"] == r["mapped"] == nodup and r["trans"] + r["cis"] == nodup
    assert (r["dangling_end"], r["self_circle"], r["religation"]) == (n["dangling_end"], n["self_circle"], n["religation"])
    assert r["same_fragment_same_strand"] == 0 and r["adjacent_other"] == 0
    assert r["cis_valid"] + r["trans"] == n["cis_valid"] + n["trans"] + extra
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
    if run["mq"] == 0:
        assert {"Yt:Z:UU", "Yt:Z:NN"} <= set(yt) and not any("M" in t[5:] for t in yt)
    else:
        assert {"Yt:Z:UU", "Yt:Z:NN", "Yt:Z:MU"} <= set(yt) or {"Yt:Z:UU", "Yt:Z:NN", "Yt:Z:UM"} <= set(yt)


@need_tools
def test_e2e_mcool(run):
    import cooler
    mc = str(run["out"] / "test.mapq_30.mcool")
    assert cooler.fileops.list_coolers(mc) == ["/resolutions/1000", "/resolutions/5000", "/resolutions/20000"]
    for res in (1000, 5000, 20000):
        c = cooler.Cooler(f"{mc}::/resolutions/{res}")
        assert c.info["sum"] == run["nodup_unique"]   # unique reads have MAPQ 60; MAPQ 0 pairs never reach the mcool
        assert "weight" in c.bins().columns
    meta = cooler.Cooler(f"{mc}::/resolutions/1000").info["metadata"]
    assert meta["min_mapq"] == 30 and meta["pairs_min_mapq"] == run["mq"] and meta["pipeline"] == "preprocess-hic"


@need_tools
def test_e2e_provenance(run):
    prov = json.loads((run["out"] / "test.preprocess-hic.provenance.json").read_text())
    assert prov["success"] is True and prov["pipeline"] == "preprocess-hic"
    assert prov["params"]["bwa_mem_options"] == "-SP5M" and prov["params"]["min_mapq"] == run["mq"]
    assert len(prov["inputs"]) == 4
    assert any(v.startswith("pairtools") and v.endswith("1.1.3") for v in prov["tool_versions"])
    assert prov["qc"]["bam"]["lossless_check"] == "pass"


@need_tools
def test_e2e_pool(run, tmp_path):
    import cooler
    mc = str(run["out"] / "test.mapq_30.mcool")
    pool = tmp_path / "pool"
    subprocess.run([sys.executable, os.path.join(BIN, "pool_mcools.py"), "--pool-id", "P", "--outdir", str(pool),
                    "--nproc", "2", "--version", "test", mc, mc], check=True)
    out = pool / "P.mapq_30.mcool"
    one = cooler.Cooler(f"{mc}::/resolutions/1000").info["sum"]
    assert cooler.fileops.list_coolers(str(out)) == ["/resolutions/1000", "/resolutions/5000", "/resolutions/20000"]
    for res in (1000, 5000, 20000):
        c = cooler.Cooler(f"{out}::/resolutions/{res}")
        assert c.info["sum"] == 2 * one and "weight" in c.bins().columns
    prov = json.loads((pool / "P.mapq_30.pool.provenance.json").read_text())
    assert prov["sum_check"] == "pass" and prov["min_mapq"] == 30 and len(prov["inputs"]) == 2
    # never overwrites
    r = subprocess.run([sys.executable, os.path.join(BIN, "pool_mcools.py"), "--pool-id", "P", "--outdir", str(pool), mc, mc])
    assert r.returncode != 0
