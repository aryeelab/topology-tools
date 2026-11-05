# preprocess-hic: Hi-C Nextflow pipeline

Processes paired-end restriction-enzyme Hi-C FASTQs (e.g. Arima Hi-C+, DpnII + HinfI) on a Slurm cluster
into a lossless BAM, MAPQ ≥ 1 contact pairs, deduplication and restriction-fragment QC, and a balanced
`.mcool`. It follows the 4DN / Open2C recipe (`bwa mem -SP5M` + pairtools, as in 4DN's docker-4dn-hic and
Open2C's distiller) and reuses the job layout of `microc-preprocess` (`MICROC-NF.md`).

## Files

| File | Purpose |
|------|---------|
| `preprocess-hic-submit.sh` | Submit one sample (all its runs) |
| `preprocess-hic-pool.sh` | Pool replicate `.mcool` files into one balanced `.mcool` (Slurm job) |
| `nf/preprocess-hic/main.nf` | Pipeline (DSL2) |
| `nf/preprocess-hic/nextflow.config` | Version, parameters, resources, `local` test profile |
| `nf/preprocess-hic/bin/` | Helpers on the tasks' `PATH`: `digest_genome.py`, `hic_frag_stats.py`, `hic_qc_summary.py`, `pool_mcools.py`, `split_pairs_by_chrom1.py` (shared with microc-preprocess) |
| `envs/preprocess-hic.yml` | Build spec of the env |
| `envs/preprocess-hic-<version>.linux-64.explicit.txt` | Exact env of each release (install from this) |
| `tests/preprocess-hic/` | Unit tests and an end-to-end test on synthetic data with known answers |
| `skills/hic-preprocessing/SKILL.md` | Agent skill: the recipe, outputs, checks and comparisons, independent of Nextflow |

## Setup (once per user)

Install a released version (tag `preprocess-hic/v<version>`) and its env; see `RELEASING.md`. The env goes
to `~/miniforge3/envs/preprocess-hic-<version>`, where the scripts look for it (override with
`PREPROCESS_HIC_ENV`). It holds Nextflow 24.10.6, OpenJDK 17, bwa-mem2 2.2.1, bwa 0.7.18, pairtools 1.1.3,
samtools 1.21, seqkit 2.14.0, cooler 0.12.1 and pairix 0.3.9 (the first six as in microc-preprocess 1.1.2).

The default aligner, bwa-mem2, needs its index files next to the FASTA (`bwa-mem2 index genome.fa` adds
`.0123` and `.bwt.2bit.64`); `ALIGNER=bwa` uses a bwa index. Build envs as batch jobs.

## Running

```bash
SAMPLEID=GSM8517989_KT_rep1 \
BWA_INDEX=/path/to/genome.fna CHROM_SIZES=/path/to/genome.chrom.sizes ASSEMBLY=mm10 \
FASTQ_GLOB="/path/fastq/{SRR30681381,SRR30681382}_{1,2}.fastq.gz" \
OUTDIR=/scratch/hic/GSM8517989_KT_rep1 TMPDIR_NF=/scratch/hic/tmp \
bash preprocess-hic-submit.sh
```

| Setting | Default | |
|---|---|---|
| `SAMPLEID`, `BWA_INDEX`, `CHROM_SIZES`, `FASTQ_GLOB`, `OUTDIR`, `TMPDIR_NF` | required | `FASTQ_GLOB` must match every run of the sample (`fromFilePairs`); runs are merged before dedup |
| `DATASET_DIR` | empty | copy outputs to `<DS>/bam`, `<DS>/pairs` (pairs and mcool) and `<DS>/qc` instead of `OUTDIR` (1.0.0 put the mcool in `<DS>/cool`) |
| `FRAG_BED` | empty | precomputed restriction fragments (BED); otherwise digest `GENOME_FASTA` (default: the `BWA_INDEX` FASTA) with `ENZYMES` |
| `ENZYMES` | `DpnII,HinfI` | Arima Hi-C+; see `digest_genome.py` for the others |
| `MCOOL_MIN_MAPQ` | `30` | MAPQ (both sides) of the pairs binned into the `.mcool` |
| `MCOOL_RESOLUTIONS` | 1 kb … 10 Mb (4DN's list) | first = base resolution |
| `ASSEMBLY` | empty | `cooler --assembly` label |
| `ALIGNER` | `bwa-mem2` | or `bwa` |
| `CHUNK_READS` | 30,000,000 | reads per FASTQ chunk |

## What it does

1. **Split** each FASTQ pair into chunks (`seqkit split2`).
2. **Align and parse** each chunk (one Slurm job each):
   `bwa-mem2 mem -SP5M | pairtools parse --min-mapq 1 --walks-policy 5unique --add-columns mapq,pos5,pos3,dist_to_5,dist_to_3,read_len | pairtools split`.
   - `-SP5M`: no mate rescue or pairing, 5'-most alignment primary, split hits flagged secondary (as 4DN);
     bwa's default score threshold (`-T 30`). pairtools 1.1.3 uses every SAM record of a read (secondary and
     supplementary alike), so `-M` changes only the flags in the BAM, not the pairs.
   - `--min-mapq 1` (pairtools' default, as 4DN): a side with MAPQ 0 is `M` (multi-mapped). The `mapq1/mapq2`
     columns let you filter further (e.g. `pairtools select '(mapq1>=30) and (mapq2>=30)'`).
   - `--walks-policy 5unique`: for a ligation walk, the 5'-most unique pair (as microc-preprocess).
     `--max-inter-align-gap` is pairtools' default (20).
   - `pairtools split` writes every SAM record of every read pair, of every pair type, to a coordinate-sorted chunk
     BAM (with pairtools' `Yt` pair-type tag), and light pairs (no SAM columns), sorted and split per chrom1.
3. **Merge** (one job per sample):
   - the chunk BAMs of all runs are merged into the **lossless BAM** (`samtools merge`), indexed; `samtools flagstat`
     must show 2 × read pairs primary records, or the job fails;
   - per chrom1, 8 at a time: `pairtools merge | pairtools dedup | pairtools restrict -f <fragments>`; dedup keeps
     pairs with both sides mapped (UU, UR, RU) and drops unmapped/multi-mapped pairs from the pairs (they stay in
     the BAM and the stats); the blocks are concatenated in sort order into a bgzipped `.pairs.gz` and indexed with
     `pairix`; stats merged with `pairtools stats --merge`.
   - restriction-fragment QC (`hic_frag_stats.py`, QC only, nothing removed): each deduplicated pair is classed with
     HiC-Pro's categories from `rfrag1/rfrag2` and the strands: dangling end (same fragment, inward `+-`),
     self-circle (same fragment, outward `-+`), same fragment same strand, religation (adjacent fragments, inward),
     adjacent other, valid cis, trans.
4. **mcool**: `pairtools select` (MAPQ ≥ `MCOOL_MIN_MAPQ` on both sides) `| cooler cload pairs` at the base resolution,
   `cooler zoomify --balance` (ICE, cooler defaults) at every resolution. The cooler metadata records the pipeline,
   version, sample and MAPQ filter.

The digest (`digest_genome.py`): fragments from cut to cut (DpnII `^GATC`, HinfI `G^ANTC`, union of sites), in
chrom-sizes order. Mouse GRCm38 with DpnII+HinfI: 14.6M fragments (median 122 bp), 1.6 min, 1.3 GB. Digest once and pass
`FRAG_BED`, or let each run digest.

## Outputs

```
DATASET_DIR/                            (or OUTDIR)
  bam/{sid}.bam{,.bai}                  # lossless: every alignment of every read pair, coordinate-sorted, Yt tags
  pairs/{sid}.pairs.gz{,.px2}           # deduplicated UU/UR/RU pairs, MAPQ >= 1, with mapq, pos5/pos3, dist_to_5/3,
                                        #   read_len and rfrag columns; pairix index
  qc/{sid}.dedup.stats.txt              # pairtools dedup stats (pre-dedup totals)
  qc/{sid}.restrict_stats.tsv           # restriction-fragment categories of the deduplicated pairs
  qc/{sid}.flagstat.tsv                 # samtools flagstat of the BAM
  qc/{sid}.hic_qc.json                  # summary: mapping, duplication, cis/trans, cis >1 kb / >20 kb, fragment
                                        #   artefact fractions, BAM lossless check
  pairs/{sid}.mapq_{Q}.mcool            # balanced multi-resolution cooler
OUTDIR/
  {sid}.preprocess-hic.provenance.json  # version + commit, params, FASTQs (path/size/mtime), tool versions,
                                        #   Slurm master job/node, times, QC summary
  {sid}.tool_versions.txt
  {sid}.trace.<timestamp>.txt           # Nextflow trace per launch
```

Unlike microc-preprocess, the BAM is not reduced to the deduplicated pairs: it keeps unmapped, multi-mapped and
duplicate reads, and duplicates are not flagged in it (dedup works on the pairs; the duplicate read IDs are the pairs
in the BAM that are not in `.pairs.gz` among UU/UR/RU).

## Pooling replicates

```bash
POOL_ID=KT OUTDIR=/scratch/hic/pooled bash preprocess-hic-pool.sh \
    rep1/pairs/rep1.mapq_30.mcool rep2/pairs/rep2.mapq_30.mcool
```
Checks the inputs come from preprocess-hic with the same MAPQ filter, chromosomes and resolutions, sums them at
the base resolution (`cooler merge`), coarsens and balances (`cooler zoomify --balance`), checks the pooled total
equals the sum of the inputs, and writes `OUTDIR/<POOL_ID>.mapq_<Q>.mcool` and `.pool.provenance.json`. Never
overwrites.

## Resources

| Process | CPUs | Memory | Notes |
|---|---|---|---|
| `split_fastq` | 4 | 16 GB | one per FASTQ pair |
| `digest_genome` | 1 | 8 GB | once per run unless `FRAG_BED` |
| `hic_align` | 6 | 28 GB | bwa-mem2 index ~19 GB for human/mouse |
| `hic_merge` | 8 | 96 GB | 8 chromosome pipelines at once; each `pairtools restrict` loads the whole fragment BED (mouse DpnII+HinfI, 14.6M fragments: ~6.5 GB peak, ~1:45 per load) |
| `hic_mcool` | 8 | 32 GB | cload + zoomify + balance at every resolution |

## Tests

```bash
PATH=<env>/bin:$PATH JAVA_HOME=<env> python -m pytest -v tests/preprocess-hic
```
The end-to-end test builds a 700 kb synthetic genome (with a duplicated segment) and 720 read pairs in two runs whose
names carry their true class (valid cis, trans, dangling end, self-circle, religation, duplicate, multi-mapped,
unmapped, half-mapped), runs the pipeline with `-profile local` and small chunks, and checks every output against
the known answers.

## Re-running, bad nodes

As microc-preprocess: re-run the same command to resume (`-resume`); failing nodes are appended to
`OUTDIR/failed_nodes_hic_<SAMPLEID>.txt` and excluded from later jobs.
