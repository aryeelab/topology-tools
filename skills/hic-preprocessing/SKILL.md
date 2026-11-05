---
name: hic-preprocessing
description: Preprocess paired-end restriction-enzyme Hi-C FASTQs (e.g. Arima Hi-C+, DpnII + HinfI) into a lossless BAM, contact pairs (MAPQ >= 1, or >= 0 to keep multi-mapping reads), dedup and restriction-fragment QC (dangling ends, self-circles, religation) and a balanced .mcool with the topology-tools preprocess-hic recipe (4DN/Open2C: bwa-mem2 mem -SP5M, pairtools parse --min-mapq 1 --walks-policy 5unique with mapq column, per-chromosome merge/dedup, pairtools restrict, pairix, cooler). Covers versions, every step and option, outputs, checks, and pooling replicates into one .mcool. Independent of workflow manager and scheduler.
---

# Hi-C preprocessing (FASTQ → BAM, pairs, mcool)

The recipe of the `preprocess-hic` component of https://github.com/aryeelab/topology-tools (tags
`preprocess-hic/vX.Y.Z`; this skill describes 1.0.0). The repo's Nextflow + SLURM implementation and its
resources are documented in `PREPROCESS-HIC.md`; the steps below work with any workflow manager.
For Micro-C use `microc-preprocessing` instead (no restriction fragments, MAPQ 30, deduplicated BAM).

## Software
| Tool | Version | Used for |
|---|---|---|
| bwa-mem2 | 2.2.1 | alignment (bwa 0.7.18 gives the same alignments) |
| pairtools | 1.1.3 | parse, sort, merge, dedup, restrict, select, split |
| samtools | 1.21 | BAM sort, merge, index, flagstat |
| seqkit | 2.14.0 | FASTQ chunks |
| cooler | 0.12.1 | cload, zoomify, balance, merge |
| pairix | 0.3.9 | `.pairs.gz` index (the binary prints 0.3.8) |

Exact env of a release: `envs/preprocess-hic-<version>.linux-64.explicit.txt` (`conda create -p <env> --file`;
behind a conda mirror, first rewrite its `https://conda.anaconda.org/` URLs to the mirror's). One env per version.

Inputs: paired FASTQs of every run/lane of a sample; a bwa-mem2 (or bwa) index; chromosome sizes of the same
genome (their order sets the pairs and cooler order); restriction fragments as BED, or the genome FASTA to digest.

## Steps
**0. Digest** (once per genome and enzyme set): `nf/preprocess-hic/bin/digest_genome.py --fasta G.fa
--chrom-sizes G.sizes --enzymes DpnII,HinfI -o frags.bed`: fragments from cut to cut, union of the enzymes' sites
(DpnII `^GATC`, HinfI `G^ANTC`; others listed in the script). Mouse GRCm38: 14.6M fragments, median 122 bp.

**1. Per chunk** (split FASTQs into chunks of equal read count; 30M reads by default):
```bash
bwa-mem2 mem -SP5M -t $T $INDEX R1.fq.gz R2.fq.gz \
 | pairtools parse --min-mapq 1 --walks-policy 5unique --chroms-path $SIZES \
     --add-columns mapq,pos5,pos3,dist_to_5,dist_to_3,read_len \
 | pairtools split --output-pairs chunk.pairs.lz4 --output-sam - \
 | samtools sort -o chunk.bam
pairtools sort chunk.pairs.lz4 > chunk.sorted.pairs   # the pipeline writes one file per chrom1
```
- `-SP5M`: no mate rescue or pairing, 5'-most alignment primary, split hits flagged secondary; bwa's default
  `-T 30`. pairtools 1.1.3 uses every SAM record of a read (secondary and supplementary), so `-M` changes only BAM flags.
- `--min-mapq 1` (default, as 4DN): a side with MAPQ 0 becomes `M` and its pair leaves the pairs. `--min-mapq 0`
  (`MIN_MAPQ=0`, from 1.1.0) keeps it as `U` at bwa's randomly chosen primary, with `mapq` 0. Filter on `mapq1/mapq2` later.
- `--walks-policy 5unique`: the 5'-most unique pair of a ligation walk. `--max-inter-align-gap` default (20).
- The chunk BAM gets every record of every read pair, of every pair type, with pairtools' `Yt` tag.

**2. Per sample** (after all runs and chunks; duplicates span lanes and runs):
```bash
samtools merge -c -p -o $SID.bam chunk*.bam && samtools index $SID.bam
# per chrom1 in parallel (duplicates always share chrom1), then concatenate in sort order:
pairtools merge chunk*.chr1.pairs | pairtools dedup --backend scipy --output-stats chr1.stats \
 | pairtools restrict -f frags.bed > chr1.pairs
cat header chr*.pairs | bgzip > $SID.pairs.gz && pairix -p pairs $SID.pairs.gz
pairtools stats --merge chr*.stats -o $SID.dedup.stats.txt
```
dedup keeps pairs with both sides mapped (UU, UR, RU) and drops unmapped/multi pairs from the pairs only.
`pairtools restrict` loads the whole BED with `np.genfromtxt` (14.6M fragments: ~6.5 GB, ~1:45 per process):
budget for it when running chromosomes in parallel.

**3. Restriction-fragment QC** (`hic_frag_stats.py`, nothing removed): each deduplicated pair, from `rfrag1/rfrag2`
and strands, in HiC-Pro's categories: dangling end (same fragment, `+-`), self-circle (same fragment, `-+`),
same fragment same strand, religation (adjacent fragments, `+-`), adjacent other, valid cis, trans.

**4. mcool**:
```bash
bgzip -dc $SID.pairs.gz | pairtools select '(mapq1>=30) and (mapq2>=30)' \
 | cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 $SIZES:1000 - base.cool
cooler zoomify --balance -r 1000,2000,5000,...,10000000 -o $SID.mapq_30.mcool base.cool
```
The pipeline records `min_mapq`, the pairs' `pairs_min_mapq`, sample and version in the cooler metadata.

**Pooling replicates** (`pool_mcools.py` / `preprocess-hic-pool.sh`): `cooler merge` the base resolutions
(same MAPQ filter, chromosomes and resolutions), `cooler zoomify --balance`; the pooled sum must equal the inputs'.
Pool after per-replicate dedup, never by deduplicating replicates together.

## Outputs
| File | Content |
|---|---|
| `<SID>.bam` + `.bai` | lossless: every alignment of every read pair, coordinate-sorted, `Yt` tags; duplicates not flagged |
| `<SID>.pairs.gz` + `.px2` | deduplicated UU/UR/RU pairs, MAPQ ≥ `--min-mapq`, mapq, pos5/pos3, dist_to_5/3, read_len, rfrag columns |
| `<SID>.dedup.stats.txt` | pairtools dedup stats with pre-dedup totals (keep the original) |
| `<SID>.restrict_stats.tsv` | restriction-fragment categories of the deduplicated pairs |
| `<SID>.hic_qc.json` | mapping, duplication, cis/trans, cis ≥1 kb / ≥20 kb, fragment artefact fractions, BAM check |
| `<SID>.mapq_<Q>.mcool` | balanced multi-resolution cooler (default Q = 30, 1 kb–10 Mb) |

## Checks
1. `hic_qc.json` `bam.lossless_check` is `pass`: BAM primary records = 2 × read pairs (the pipeline fails otherwise).
2. The pairs header has `-SP5M` on the bwa `@PG`, `--min-mapq` (1 or 0) `--walks-policy 5unique` and the add-columns on parse,
   pairtools 1.1.3, and a `pairtools_restrict` `@PG`.
3. `total_nodups` in the dedup stats = data lines in `.pairs.gz` = `restrict_stats` `total`.
4. The mcool metadata `min_mapq` is the threshold you meant; pool only mcools with the same one.
