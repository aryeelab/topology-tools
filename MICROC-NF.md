# Micro-C Nextflow Pipeline

Processes paired-end Micro-C FASTQ files through alignment, deduplication, and output of `.mapped.pairs` + BAM on a Slurm cluster.

## Files

| File | Purpose |
|------|---------|
| `microc-nf-submit.sh` | Submit script — the only file a user edits |
| `nf/microc_preprocess.nf` | Pipeline definition (DSL1) |
| `nf/nextflow.config` | Resource settings and Slurm configuration |
| `nf/bin/split_pairs_by_chrom1.py` | Writes a chunk's sorted pairs as one file per chrom1 (on the tasks' `PATH`) |
| `envs/microc-preprocess.yml` | Build spec for the pipeline's conda env |
| `envs/microc-preprocess-<version>.linux-64.explicit.txt` | Exact env each release was tested in (install from this) |
| `microc-env.yml` | Legacy env (pairtools 0.3.0); only for reproducing pairs made before 1.0.0 |

## Setup (once per user)

Install a released version (tag `microc-preprocess/v<version>`) and its env; see `RELEASING.md`,
"Installing a released version". The env goes to `~/miniforge3/envs/microc-preprocess-<version>`,
where the submit script looks for it (override with `MICROC_ENV=/path/to/env`). It holds Nextflow
24.10.6, OpenJDK 17, bwa-mem2 2.2.1, BWA 0.7.18, pairtools 1.1.3, samtools and seqkit. No other software is needed.

From 1.1.0 the default aligner is bwa-mem2, which needs its own index files next to the FASTA:
`bwa-mem2 index genome.fa` adds `genome.fa.0123` and `genome.fa.bwt.2bit.64` (human: ~1 h, ~72 GB RAM).
Its `.amb/.ann/.pac` are the same as bwa's, so one prefix can serve both aligners. A `.tar.gz` index must
contain them too. `ALIGNER=bwa` runs bwa 0.7.18 with a bwa index, as up to 1.0.3.

Build envs as a batch job: on network storage one can take tens of minutes.

## Running

1. Edit the **User configuration** block at the top of `microc-nf-submit.sh`:

   ```bash
   SAMPLEID=221014_LPS141_MicroC_RepB          # sample identifier                       (required)
   BWA_INDEX=/path/to/hg38.tgz                 # BWA index tarball or indexed FASTA       (required)
   CHROM_SIZES=/path/to/hg38.chrom.sizes       # chromosome sizes file                    (required)
   FASTQ_GLOB="/path/to/fastq/${SAMPLEID}/*R{1,2}*.fastq.gz"                             # (required)
   OUTDIR=/path/to/scratch/${SAMPLEID}         # run dir: provenance, manifest (shared FS) (required)
   TMPDIR_NF=/path/to/scratch/tmp              # temp dir (shared filesystem)             (required)
   DATASET_DIR=/path/to/<dataset>              # dataset dir; empty = keep outputs in OUTDIR
   ASSAY=microc                                # microc | rcmc
   MAPQ=30                                     # pairtools parse --min-mapq
   CHUNK_READS=30000000                        # reads per FASTQ chunk
   ALIGNER=bwa-mem2                            # bwa-mem2 | bwa
   DEDUP_BACKEND=scipy                         # pairtools dedup --backend: scipy | cython
   ```

2. Submit:

   ```bash
   bash microc-nf-submit.sh
   ```

Each of these can also be set in the environment instead (`SAMPLEID=… bash microc-nf-submit.sh`).
There are no site defaults: the script stops if a required setting is empty. It checks the env has pairtools 1.1.3, manages the Slurm master job, and passes all paths to Nextflow. No other configuration is needed.

## What it does

1. Splits each FASTQ pair into 30M-read chunks in its own Slurm job (`split_fastq`: `seqkit split2`, chunks
   gzipped at level 1). Up to 1.0.2 Nextflow's `splitFastq` did this inside the 4 GB master job, which was
   killed on full-size libraries: the page cache of the uncompressed chunk writes counts against the job's
   memory limit. The master now only orchestrates.
2. Aligns each chunk in parallel: `bwa-mem2 mem -5SP -T0 | pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 --add-columns pos5,pos3,dist_to_5,dist_to_3,read_len | pairtools split`.
   The pairs (without SAM columns) are sorted and written as one file per chrom1 (`nf/bin/split_pairs_by_chrom1.py`);
   the alignments of UU/UR/RU pairs go to a coordinate-sorted chunk BAM.
3. Merge job: each chrom1 is merged across chunks and deduplicated on its own (`pairtools merge | pairtools dedup
   --mark-dups`, 8 at a time); the results are concatenated in sort order and the stats combined
   (`pairtools stats --merge`). Meanwhile the chunk BAMs are merged; the BAM is then filtered to the read IDs of the
   deduplicated pairs (`samtools view -N`: keep the kept read IDs or drop the duplicates', whichever list is
   shorter; from 1.1.1) and indexed.
4. Copies pairs, BAM and dedup stats to `DATASET_DIR/{pairs,bam,qc}` (or `OUTDIR`), never overwriting
   existing files, and writes provenance to `OUTDIR`

Intermediate files (FASTQ chunks, per-chunk pairs, merged pairs) go to the Nextflow work dir,
`OUTDIR/work` (or `$NXF_WORK`). For a deep library this is hundreds of GB: put `OUTDIR` on scratch and delete
`work/` once the outputs are checked.

The parse options follow the FOCI-3D looping-factors paper (aryeelab/foci-3d-looping-factors, `Protocol.Rmd`),
except MAPQ (default 30; the paper used 40; this pipeline used 20 with pairtools 0.3.0 before 1.0.0)
and the extra `dist_to_5,dist_to_3,read_len` columns, which microc-qc uses. The MAPQ a pairs file was made with
is in its header (`#samheader: @PG ID:pairtools_parse … --min-mapq N`).

From 1.0.0 the pipeline uses pairtools 1.1.3 instead of 0.3.0. On 4M read pairs of Caov-3 at MAPQ 30
(one shared bwa alignment), every pair kept by both versions was identical, 0.06% were kept by only one
(a different duplicate chosen), and microc-qc's ≤80 bp fraction was 13.10% vs 13.11%. pairtools 1.x
labels walk-rescued pairs `UU` where 0.3.0 used `uu`/`Uu`/`uU`; filter on `pair_type` with that in mind.

### Changes in 1.1.0 (speed)
On a 182M read-pair RCMC library a sample took 2 h 36 min instead of 6 h 49 min with 1.0.3; the merge job
took 20 min instead of 2 h 14 min.
- **bwa-mem2** replaces bwa: identical pairs and stats on the full library, align jobs ~2x faster.
- **Dedup on pairs without SAM columns**: in 1.0.3 dedup parsed every pair's SAM records (~4x slower); the
  BAM is now built from per-chunk sorted alignments, filtered to the kept read IDs: the same records as 1.0.3.
- **Dedup per chrom1 in parallel** (duplicates always share chrom1). scipy's dedup clusters pairs in chunks of
  10,000, and splitting by chromosome moves those chunk boundaries: on the test library 1,078 of 57.6M kept pairs
  (0.002%) differ from 1.0.3 (another duplicate kept, or a pair counted as duplicate), 1,064 of them in the deeply
  covered capture regions; 34 fewer pairs in total. Everything else is identical.
- The default dedup backend is scipy again (cython is not transitive and slower on these pairs).
- Inputs given to pairtools as `.gz` paths are decompressed with pbgzip when it is installed; its threads spin on
  locks (system time > CPU time) without being faster, so the pipeline pipes data to pairtools or uses lz4.

Align jobs (6 CPUs each) run in parallel up to your cluster's limits; a 400-CPU QOS allows ~66.

## Outputs

```
DATASET_DIR/                       (or OUTDIR when DATASET_DIR is empty)
  pairs/{sample_id}.mapped.pairs   # deduplicated contact pairs
  bam/{sample_id}.bam{,.bai}       # coordinate-sorted BAM
  qc/{sample_id}.stats.txt         # pairtools dedup statistics (pre-dedup counts; pass to microc-qc --stats)
OUTDIR/
  {sample_id}.preprocess.provenance.json  # pipeline version + commit, params, input FASTQs, tool versions, run times
  {sample_id}.tool_versions.txt
  {sample_id}.foci_manifest.tsv           # one FOCI-3D batch manifest row (only with DATASET_DIR)
  {sample_id}.trace.<YYYYMMDD-HHMMSS>.txt # Nextflow trace, one per launch (from 1.1.2)
```

The trace has one row per task attempt: Slurm job id, submit/start/complete times, queue and run time,
CPUs, `%cpu`, peak RSS and virtual memory, and I/O (`rchar`/`wchar`, `read_bytes`/`write_bytes`),
with durations in ms and sizes in bytes (`raw`). Keep it when you delete `work/`: it is the record of
how long each step took and what it used, for sizing later runs.

Outputs are copied, not symlinked: earlier versions used Nextflow's default `symlink` publish mode, so files in
`OUTDIR` were links into the Nextflow work directory on scratch. Check old runs with `ls -l` before deleting
`work/`.

The manifest row is `dataset_dir  pairs  sample_id  assay  stats` (tab-separated), the input format of a
FOCI-3D counts + QC batch run; concatenate rows for several samples. Pass `qc/{sample_id}.stats.txt` to
`microc-qc.py --stats` (it holds pre-dedup totals; stats regenerated from `.mapped.pairs` do not).

## Re-running / resuming

Each sample runs Nextflow in its own `OUTDIR` (cache, `.nextflow.log`, master log `<SAMPLEID>.master.<jobid>.log`, trace), so several samples can run at once from one install; before 1.0.2 a second concurrent sample failed with `Unable to acquire lock on session`.

Re-run the same command (same `OUTDIR`, so the same work dir). Nextflow's `-resume` flag is set by default — completed tasks are skipped. From 1.0.3 the chunks are task outputs, so a resumed run reuses finished `split_fastq` and `microc_align` tasks.

## Bad node exclusion

On job failure the ERR trap appends the node hostname to:
```
OUTDIR/failed_nodes_microc_{SAMPLEID}.txt
```
All subsequent jobs (master and children) exclude those nodes. To manually exclude a node:
```bash
echo node04 >> /path/to/OUTDIR/failed_nodes_microc_{SAMPLEID}.txt
```
To clear the blacklist, delete or empty the file.

## Resource settings

| Process | CPUs | Memory | Notes |
|---------|------|--------|-------|
| `split_fastq` | 4 | 16 GB | one per FASTQ pair; ~100 s per 30M read pairs |
| `microc_align` | 6 | 28 GB | bwa-mem2 dominates (~19 GB for a human index; bwa ~7 GB); ~1.1–1.8 h per 30M read pairs, by node |
| `mergepairs` | 8 | 64 GB | 8 chromosome dedups at once; `samtools view -N` holds the shorter of the kept / duplicate read-ID lists (~50–75 B per read name; at most half the mapped pairs). 182M read-pair library: 14 GB peak, 20 min |

Settings are in `nf/nextflow.config` under `withName:split_fastq`, `withName:microc_align` and `withName:mergepairs`. The master job (Nextflow itself) needs 4 GB.

## Troubleshooting

**"conda env … not found"** — install the release env (Setup above) or set `MICROC_ENV`.

**`…0123 not found`** — build the bwa-mem2 index (`bwa-mem2 index <fasta>`, Setup above) or set `ALIGNER=bwa`.

**All jobs land on a bad node** — pre-seed the exclusion file before submitting (see Bad node exclusion above).

**`pairtools dedup` crashes with `AttributeError: np.int`** — this is pairtools 0.3.0 with numpy ≥1.24 (legacy env only). `microc-env.yml` pins numpy to 1.23.5.

**Jobs stay PENDING** — check `squeue -u $USER`. If too many jobs are queued, the 400-CPU QOS limit may be reached; they will start as running jobs finish.
