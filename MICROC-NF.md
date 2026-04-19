# Micro-C Nextflow Pipeline

Processes paired-end Micro-C FASTQ files through alignment, deduplication, and output of `.mapped.pairs` + BAM on a Slurm cluster.

## Files

| File | Purpose |
|------|---------|
| `microc-nf-submit.sh` | Submit script — the only file a user edits |
| `nf/microc_preprocess.nf` | Pipeline definition (DSL1) |
| `nf/nextflow.config` | Resource settings and Slurm configuration |
| `microc-env.yml` | Pinned conda environment (`name: microc`) |

## Setup (once per user)

```bash
conda env create -f microc-env.yml
```

Requires conda (miniforge/mamba). Installs Nextflow, BWA, pairtools 0.3.0, samtools, and numpy 1.23.5 into an env named `microc`. No other software is needed.

## Running

1. Edit the **User configuration** block at the top of `microc-nf-submit.sh`:

   ```bash
   SAMPLEID=221014_LPS141_MicroC_RepB          # sample identifier
   BWA_INDEX=/path/to/hg38.tgz                 # BWA index tarball
   CHROM_SIZES=/path/to/hg38.chrom.sizes       # chromosome sizes file
   FASTQ_GLOB="/path/to/fastq/${SAMPLEID}/*R{1,2}*.fastq.gz"
   OUTDIR=/cluster/aryeelab/${USER}/${SAMPLEID} # output (shared filesystem)
   TMPDIR_NF=/cluster/aryeelab/${USER}/tmp      # temp dir (shared filesystem)
   ```

2. Submit:

   ```bash
   bash microc-nf-submit.sh
   ```

The script auto-detects the `microc` conda env, manages the Slurm master job, and passes all paths to Nextflow. No other configuration is needed.

## What it does

1. Splits input FASTQs into 30M-read chunks (`splitFastq`)
2. Aligns each chunk in parallel: `bwa mem | pairtools parse | pairtools sort`
3. Merges all chunks: `pairtools merge | dedup | split`
4. Outputs final files to `OUTDIR`

Up to ~66 align jobs run in parallel (400 CPU QOS limit ÷ 6 CPUs/job).

## Outputs

```
OUTDIR/
  {sample_id}.mapped.pairs    # deduplicated contact pairs
  {sample_id}.bam             # coordinate-sorted BAM
  {sample_id}.bam.bai         # BAM index
  {sample_id}.stats.txt       # pairtools dedup statistics
```

## Re-running / resuming

Re-run the same command. Nextflow's `-resume` flag is set by default — completed tasks are skipped. Note: `splitFastq` creates new temp files each run, so align jobs are not cached across submissions; only `mergepairs` benefits from resume.

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

**"microc conda environment not found"** — run `conda env create -f microc-env.yml`.

**All jobs land on a bad node** — pre-seed the exclusion file before submitting (see Bad node exclusion above).

**`pairtools dedup` crashes with `AttributeError: np.int`** — numpy version is wrong. The `microc` env pins numpy to 1.23.5; recreate the env from `microc-env.yml`.

**Jobs stay PENDING** — check `squeue -u $USER`. If too many jobs are queued, the 400-CPU QOS limit may be reached; they will start as running jobs finish.
