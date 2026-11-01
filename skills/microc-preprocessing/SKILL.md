---
name: microc-preprocessing
description: Preprocess paired-end Micro-C or Region Capture Micro-C (RCMC) FASTQs into deduplicated .mapped.pairs, a sorted BAM and pairtools dedup stats with the topology-tools microc-preprocess recipe (bwa mem -5SP, pairtools parse MAPQ 30 / walks-policy 5unique with pos5, pos3, dist and read_len columns, merge, dedup, split). Covers the software versions, every step and option, parallelising by chunk, the outputs, checks, and comparing pairs made with older recipes (MAPQ 20/40, pairtools 0.3.0 vs 1.x). Independent of workflow manager and scheduler.
---

# Micro-C preprocessing (FASTQ → pairs)

The recipe of the `microc-preprocess` component of https://github.com/aryeelab/topology-tools
(releases are tags `microc-preprocess/vX.Y.Z`; this skill describes 1.0.3). It is written so you can
run it with whatever workflow manager and scheduler you use.

## Software
| Tool | Version | Used for |
|---|---|---|
| bwa | 0.7.18 | alignment |
| pairtools | 1.1.3 | parse, sort, merge, dedup, split |
| samtools | 1.21 | BAM sort and index |
| seqkit | 2.14.0 | splitting FASTQs into chunks (`seqkit split2`) |

Exact environment of a release: `envs/microc-preprocess-<version>.linux-64.explicit.txt` in the repo
(`conda create -p <env> --file <that file>`); build spec `envs/microc-preprocess.yml`. Install into a
directory you control, one env per version. Record `pairtools --version`, bwa's `Version:` line and
`samtools --version | head -1` with every output.

Inputs: paired FASTQs (all lanes of a sample), a BWA index of the genome, and a chromosome sizes file
for the same genome (its order sets the pairs header and sort order).

## Steps
**1. Align and parse** (per FASTQ pair, or per chunk of it; see Parallelising):
```bash
bwa mem -5SP -T0 -t $T $INDEX R1.fastq.gz R2.fastq.gz \
 | pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 \
     --add-columns pos5,pos3,dist_to_5,dist_to_3,read_len \
     --nproc-in $((T/2)) --nproc-out $((T/2)) --chroms-path $CHROM_SIZES \
 | pairtools sort --nproc $T --tmpdir $TMP -o chunk.pairsam.gz
```
- `-5SP`: 5'-most alignment primary, no mate rescue or pairing (Hi-C style). `-T0`: report all alignments
  and let pairtools filter.
- `--min-mapq 30`: the default of this recipe (see Comparing below). `--walks-policy 5unique` keeps the
  5'-most unique pair of a ligation walk.
- `pos5,pos3` are required by `foci-3d count`; `dist_to_5,dist_to_3,read_len` by microc-qc.

**2. Merge, deduplicate, split** (once per sample, over all chunks and lanes):
```bash
pairtools merge --tmpdir $TMP --nproc 4 -o merge.pairs.gz chunk*.pairsam.gz
pairtools dedup merge.pairs.gz --nproc-in 4 --nproc-out 4 --mark-dups --output-stats $SID.stats.txt \
 | pairtools split --nproc-in 4 --nproc-out 4 --output-pairs $SID.mapped.pairs --output-sam - \
 | samtools view -bS -@4 | samtools sort -@4 -o $SID.bam
samtools index $SID.bam
```
Deduplicate only after merging every lane and chunk of a sample: duplicates span lanes.

## Parallelising
Split each FASTQ pair into chunks of equal read count (the release default is 30,000,000 reads;
`seqkit split2 -1 R1 -2 R2 -s 30000000 -j 4 --compress-level 1 -O chunks`), run
step 1 on each chunk in parallel, then step 2 once. bwa estimates the insert-size distribution per batch,
so different chunk sizes change a few alignments per 10,000 pairs (~99.99% identical); keep the chunk size
fixed within a project. Write chunks compressed: uncompressed chunks of a deep library are hundreds of GB
of writes. Resources used by the reference implementation: splitting 4 CPUs (~100 s per 30M read pairs);
step 1 6 CPUs / 10 GB per chunk;
step 2 8 CPUs / 100 GB (k-way merge and dedup of a deep library). `$TMP` must have room for the merged
pairs.

## Outputs
| File | Content |
|---|---|
| `<SID>.mapped.pairs` | deduplicated pairs (pairtools 1.x format, header with `@PG` lines) |
| `<SID>.bam` + `.bai` | coordinate-sorted alignments of the kept pairs |
| `<SID>.stats.txt` | dedup stats with **pre-dedup** totals: pass to `microc-qc.py --stats` |

Stats regenerated later from `.mapped.pairs` (`pairtools stats`) describe deduplicated pairs only
(no total reads, mapping or duplication rates). Keep the original.

## Checks
1. Every FASTQ of the sample was used: no empty glob, both mates of every lane.
2. `grep -m1 'ID:pairtools_parse' <SID>.mapped.pairs` shows `--min-mapq 30`, `--walks-policy 5unique`,
   the add-columns, and pairtools 1.1.3.
3. The pairs file ends in a newline and its last line has all columns (truncated copies happen).
4. `total_nodups` in the stats is close to the number of data lines in the pairs file.

## Comparing with older pairs
- The MAPQ a pairs file was made with is in its header (`--min-mapq N`). Check before pooling or
  comparing samples. MAPQ 30 drops ~2% of pairs relative to 20.
- The FOCI-3D looping-factors protocol (aryeelab/foci-3d-looping-factors, `Protocol.Rmd`) uses the same
  steps at MAPQ 40 and without `dist_to_5,dist_to_3,read_len`. microc-preprocess before 1.0.0 used
  pairtools 0.3.0 at MAPQ 20.
- Same alignment, MAPQ 30: pairtools 1.1.3 and 0.3.0 keep identical pairs except ~0.06% (a different
  duplicate chosen); the ≤80 bp fragment fraction is unchanged (13.10% vs 13.11% on a Caov-3 test).
- pairtools 1.x labels walk-rescued pairs `UU`; 0.3.0 used `uu`/`Uu`/`uU`. foci-3d and microc-qc do not
  filter on `pair_type`; your own filters may.
- pairtools 0.3.0 needs numpy < 1.24 (`np.int` error in dedup otherwise).
