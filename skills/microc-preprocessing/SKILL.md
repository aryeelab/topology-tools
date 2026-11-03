---
name: microc-preprocessing
description: Preprocess paired-end Micro-C or Region Capture Micro-C (RCMC) FASTQs into deduplicated .mapped.pairs, a sorted BAM and pairtools dedup stats with the topology-tools microc-preprocess recipe (bwa-mem2 mem -5SP, pairtools parse MAPQ 30 / walks-policy 5unique with pos5, pos3, dist and read_len columns, per-chromosome merge and dedup, BAM from per-chunk sorted alignments). Covers the software versions, every step and option, parallelising by chunk, the outputs, checks, and comparing pairs made with older recipes (MAPQ 20/40, pairtools 0.3.0 vs 1.x). Independent of workflow manager and scheduler.
---

# Micro-C preprocessing (FASTQ → pairs)

The recipe of the `microc-preprocess` component of https://github.com/aryeelab/topology-tools
(releases are tags `microc-preprocess/vX.Y.Z`; this skill describes 1.1.0). It is written so you can
run it with whatever workflow manager and scheduler you use.

## Software
| Tool | Version | Used for |
|---|---|---|
| bwa-mem2 | 2.2.1 | alignment (default from 1.1.0) |
| bwa | 0.7.18 | alignment up to 1.0.3 (same alignments, ~2x slower) |
| pairtools | 1.1.3 | parse, sort, merge, dedup, split |
| samtools | 1.21 | BAM sort, merge, filter, index |
| seqkit | 2.14.0 | splitting FASTQs into chunks (`seqkit split2`) |

Exact environment of a release: `envs/microc-preprocess-<version>.linux-64.explicit.txt` in the repo
(`conda create -p <env> --file <that file>`); build spec `envs/microc-preprocess.yml`. Install into a
directory you control, one env per version. Record `pairtools --version`, `bwa-mem2 version` (or bwa's `Version:` line) and
`samtools --version | head -1` with every output.

Inputs: paired FASTQs (all lanes of a sample), a bwa-mem2 index of the genome (`bwa-mem2 index genome.fa`
adds `.0123` and `.bwt.2bit.64`; human: ~1 h, ~72 GB RAM; `.amb/.ann/.pac` are shared with a bwa index), and a chromosome sizes file
for the same genome (its order sets the pairs header and sort order).

## Steps
**1. Align, parse, split pairs from alignments** (per FASTQ pair, or per chunk of it; see Parallelising):
```bash
bwa-mem2 mem -5SP -T0 -t $T $INDEX R1.fastq.gz R2.fastq.gz \
 | pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 \
     --add-columns pos5,pos3,dist_to_5,dist_to_3,read_len --chroms-path $CHROM_SIZES \
 | pairtools split --output-pairs chunk.unsorted.pairs.lz4 --output-sam - \
 | samtools view -u -e '[Yt]=="UU" || [Yt]=="UR" || [Yt]=="RU"' - \
 | samtools sort -@2 -m 1G -T $TMP/chunk -o chunk.bam
pairtools sort --nproc $T --memory 2G --tmpdir $TMP chunk.unsorted.pairs.lz4 \
 | split_pairs_by_chrom1.py chunk 2    # repo nf/bin/: chunk.blk_<chrom1>.pairs.gz, chunk.blocks.tsv
```
- `-5SP`: 5'-most alignment primary, no mate rescue or pairing (Hi-C style). `-T0`: report all alignments
  and let pairtools filter. bwa-mem2 gives the same alignments as `bwa mem` 0.7.18 (identical pairs on a
  full library) about twice as fast, but holds ~19 GB for a human index.
- `--min-mapq 30`: the default of this recipe (see Comparing below). `--walks-policy 5unique` keeps the
  5'-most unique pair of a ligation walk.
- `pos5,pos3` are required by `foci-3d count`; `dist_to_5,dist_to_3,read_len` by microc-qc.
- Pairs and alignments go separate ways: dedup on pairs without SAM columns is ~4x faster. Only the
  alignments of `UU/UR/RU` pairs (the types dedup can keep; pairtools tags each record `Yt`) go to the BAM.
- Don't give pairtools `.gz` paths when pbgzip is installed: pairtools then decompresses with pbgzip, whose
  threads spin on locks (system time above CPU time, no speed gain). Pipe the data in, or use `.lz4`.

**2. Merge and deduplicate, one chrom1 at a time** (once per sample, over all chunks and lanes):
```bash
# for each chrom1 block b, several at once:
pairtools merge --tmpdir $TMP --nproc 1 chunk*.blk_$b.pairs.gz \
 | pairtools dedup --backend scipy --mark-dups --output-stats blk/$b.stats -o blk/$b.pairs
# $SID.mapped.pairs: header of the first block, then the data lines of every block in chrom1 order
# (C collation, as pairtools sort): LC_ALL=C sort -t$'\t' -k2,2 -u chunk*.blocks.tsv gives the order
pairtools stats --merge blk/*.stats -o $SID.stats.txt
```
- Duplicates always share chrom1, so the blocks are independent. scipy dedup clusters pairs in chunks of
  10,000, though, and splitting moves the chunk boundaries: ~2 per 100,000 kept pairs differ from one
  whole-sample dedup, mostly where coverage is very deep (capture regions). For exact agreement with
  microc-preprocess <= 1.0.3, merge all blocks into one `pairtools dedup` (~6x slower).
- Deduplicate only after merging every lane and chunk of a sample: duplicates span lanes.

**3. BAM** (the merge can run during step 2):
```bash
samtools merge -@4 -c -p -l 1 -o all.bam chunk*.bam
grep -v '^#' $SID.mapped.pairs | cut -f1 > nodup_ids.txt
samtools view -@8 -N nodup_ids.txt -b -o $SID.bam all.bam
samtools index $SID.bam
```
The result holds the alignments of the kept pairs, as in earlier releases (same records; reads at the same
position may be in a different order). The read-ID list needs ~3 GB of memory per 60M kept pairs.

## Parallelising
Split each FASTQ pair into chunks of equal read count (the release default is 30,000,000 reads;
`seqkit split2 -1 R1 -2 R2 -s 30000000 -j 4 --compress-level 1 -O chunks`), run
step 1 on each chunk in parallel, then step 2 once. bwa estimates the insert-size distribution per batch,
so different chunk sizes change a few alignments per 10,000 pairs (~99.99% identical); keep the chunk size
fixed within a project. Write chunks compressed: uncompressed chunks of a deep library are hundreds of GB
of writes. Resources used by the reference implementation: splitting 4 CPUs (~100 s per 30M read pairs);
step 1 6 CPUs / 28 GB per chunk (bwa-mem2 ~19 GB; ~1.1–1.8 h per 30M read pairs); steps 2–3 8 CPUs / 64 GB
(eight blocks at once; 182M read pairs: 20 min, 14 GB peak). The work dir needs room for the chunk BAMs
and the merged BAM (several times the final BAM).

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
- microc-preprocess 1.1.0 vs 1.0.3 (182M read-pair RCMC library): alignments identical (bwa-mem2 vs bwa);
  0.002% of kept pairs differ through the per-chromosome dedup (34 fewer in total, nearly all in capture
  regions); duplication rate 61.04298% vs 61.04301%. Treat the two as the same recipe.
- pairtools 0.3.0 needs numpy < 1.24 (`np.int` error in dedup otherwise).
