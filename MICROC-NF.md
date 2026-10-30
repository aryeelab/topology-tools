# Micro-C Nextflow Pipeline

Processes paired-end Micro-C FASTQ files through alignment, deduplication, and output of `.mapped.pairs` + BAM on a Slurm cluster.

## Files

| File | Purpose |
|------|---------|
| `microc-nf-submit.sh` | Submit script — the only file a user edits |
| `nf/microc_preprocess.nf` | Pipeline definition (DSL1) |
| `nf/nextflow.config` | Resource settings and Slurm configuration |
| `envs/microc-preprocess.yml` | Build spec for the pipeline's conda env |
| `envs/microc-preprocess-<version>.linux-64.explicit.txt` | Exact env each release was tested in (install from this) |
| `microc-env.yml` | Legacy env (pairtools 0.3.0); only for reproducing pairs made before 1.0.0 |

## Setup (once per user)

Install a released version (tag `microc-preprocess/v<version>`) and its env; see `RELEASING.md`,
"Installing a released version". The env goes to `~/miniforge3/envs/microc-preprocess-<version>`,
where the submit script looks for it (override with `MICROC_ENV=/path/to/env`). It holds Nextflow
24.10.6, OpenJDK 17, BWA 0.7.18, pairtools 1.1.3 and samtools. No other software is needed.

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
   ```

2. Submit:

   ```bash
   bash microc-nf-submit.sh
   ```

Each of these can also be set in the environment instead (`SAMPLEID=… bash microc-nf-submit.sh`).
There are no site defaults: the script stops if a required setting is empty. It checks the env has pairtools 1.1.3, manages the Slurm master job, and passes all paths to Nextflow. No other configuration is needed.

## What it does

1. Splits input FASTQs into 30M-read chunks (`splitFastq`)
2. Aligns each chunk in parallel: `bwa mem -5SP -T0 | pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 --add-columns pos5,pos3,dist_to_5,dist_to_3,read_len | pairtools sort`
3. Merges all chunks: `pairtools merge | dedup --mark-dups | split`
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
```

Outputs are copied, not symlinked: earlier versions used Nextflow's default `symlink` publish mode, so files in
`OUTDIR` were links into the Nextflow work directory on scratch. Check old runs with `ls -l` before deleting
`work/`.

The manifest row is `dataset_dir  pairs  sample_id  assay  stats` (tab-separated), the input format of a
FOCI-3D counts + QC batch run; concatenate rows for several samples. Pass `qc/{sample_id}.stats.txt` to
`microc-qc.py --stats` (it holds pre-dedup totals; stats regenerated from `.mapped.pairs` do not).

## Re-running / resuming

Re-run the same command (same `OUTDIR`, so the same work dir). Nextflow's `-resume` flag is set by default — completed tasks are skipped. Note: `splitFastq` creates new temp files each run, so align jobs are not cached across submissions; only `mergepairs` benefits from resume.

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
| `microc_align` | 6 | 10 GB | BWA dominates; parse nproc scaled to cores/2 |
| `mergepairs` | 8 | 100 GB | pairtools dedup is I/O-bound; memory for k-way merge |

Settings are in `nf/nextflow.config` under `withName:microc_align` and `withName:mergepairs`.

## Troubleshooting

**"conda env … not found"** — install the release env (Setup above) or set `MICROC_ENV`.

**All jobs land on a bad node** — pre-seed the exclusion file before submitting (see Bad node exclusion above).

**`pairtools dedup` crashes with `AttributeError: np.int`** — this is pairtools 0.3.0 with numpy ≥1.24 (legacy env only). `microc-env.yml` pins numpy to 1.23.5.

**Jobs stay PENDING** — check `squeue -u $USER`. If too many jobs are queued, the 400-CPU QOS limit may be reached; they will start as running jobs finish.
