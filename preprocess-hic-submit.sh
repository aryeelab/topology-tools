#!/bin/bash
# preprocess-hic-submit.sh — submit the Hi-C preprocessing Nextflow pipeline (nf/preprocess-hic) on Slurm
#
# One sample per submission: all runs/lanes of the sample are merged before deduplication.
# Requirements: the preprocess-hic conda env for this release, built from
#   envs/preprocess-hic-<version>.linux-64.explicit.txt (see RELEASING.md).
#
# Usage: SAMPLEID=... BWA_INDEX=... CHROM_SIZES=... FASTQ_GLOB=... OUTDIR=... TMPDIR_NF=... bash preprocess-hic-submit.sh

set -euo pipefail

# ============================================================
# User configuration — edit these, or set them in the environment. The first six are required.
# ============================================================
SAMPLEID="${SAMPLEID:-}"                 # sample identifier
BWA_INDEX="${BWA_INDEX:-}"               # indexed FASTA (bwa and/or bwa-mem2 index files next to it), or a .tar.gz
CHROM_SIZES="${CHROM_SIZES:-}"           # chromosome sizes of the same genome; its order sets pairs and cooler order
# glob over every FASTQ of the sample, e.g. "/path/fastq/{SRR1,SRR2}_{1,2}.fastq.gz" or "/path/fastq/S1_*_R{1,2}_001.fastq.gz";
# files are paired by run name (<run>_1/_2 or _R1/_R2, optional _001), and each run needs exactly one R1 and one R2
FASTQ_GLOB="${FASTQ_GLOB:-}"
OUTDIR="${OUTDIR:-}"                     # run dir (shared FS): provenance, versions, trace; outputs unless DATASET_DIR
TMPDIR_NF="${TMPDIR_NF:-}"               # temp dir for pairtools sort/merge (shared FS)

# Optional: dataset directory (<DS>): BAM to <DS>/bam, pairs and mcool to <DS>/pairs, QC to <DS>/qc.
DATASET_DIR="${DATASET_DIR:-}"

# Restriction fragments for pairtools restrict (QC only): a precomputed BED (FRAG_BED), or an in-silico digest
# of GENOME_FASTA (default: the BWA_INDEX FASTA) with ENZYMES (default Arima Hi-C+: DpnII,HinfI).
FRAG_BED="${FRAG_BED:-}"
GENOME_FASTA="${GENOME_FASTA:-}"
ENZYMES="${ENZYMES:-DpnII,HinfI}"

MIN_MAPQ="${MIN_MAPQ:-1}"               # pairtools parse --min-mapq: 1 (sides with MAPQ 0 become M and their pairs
                                         # leave the pairs file) or 0 (keep them: bwa's random primary of a multi-mapping
                                         # read becomes a U side with mapq 0; filter on mapq1/mapq2 downstream)
MCOOL_MIN_MAPQ="${MCOOL_MIN_MAPQ:-30}"   # MAPQ (both sides) of the pairs binned into the .mcool
MCOOL_RESOLUTIONS="${MCOOL_RESOLUTIONS:-1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000}"
ASSEMBLY="${ASSEMBLY:-}"                 # cooler --assembly label, e.g. mm10
ALIGNER="${ALIGNER:-bwa-mem2}"           # bwa-mem2 (needs <BWA_INDEX>.0123 and .bwt.2bit.64) | bwa
CHUNK_READS="${CHUNK_READS:-30000000}"   # reads per FASTQ chunk for parallel alignment

PIPELINE_VERSION=1.1.1
PREPROCESS_HIC_ENV=${PREPROCESS_HIC_ENV:-${HOME}/miniforge3/envs/preprocess-hic-${PIPELINE_VERSION}}
# ============================================================

missing=""
for v in SAMPLEID BWA_INDEX CHROM_SIZES FASTQ_GLOB OUTDIR TMPDIR_NF; do
    [ -n "${!v}" ] || missing="${missing} ${v}"
done
if [ -n "${missing}" ]; then
    echo "Error: set${missing} (edit the User configuration block or export them)."
    exit 1
fi

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export PREPROCESS_HIC_ENV
E="${PREPROCESS_HIC_ENV}/bin"
if [ ! -x "${E}/pairtools" ]; then
    echo "Error: conda env ${PREPROCESS_HIC_ENV} not found."
    echo "Create it with: conda create -p ${PREPROCESS_HIC_ENV} --file ${REPO}/envs/preprocess-hic-${PIPELINE_VERSION}.linux-64.explicit.txt"
    exit 1
fi
if ! "${E}/pairtools" --version | grep -q ' 1\.1\.3$'; then
    echo "Error: ${PREPROCESS_HIC_ENV} does not have pairtools 1.1.3 (preprocess-hic ${PIPELINE_VERSION} needs it)."
    exit 1
fi
for t in seqkit cooler pairix bgzip samtools nextflow; do
    [ -x "${E}/${t}" ] || { echo "Error: ${PREPROCESS_HIC_ENV} has no ${t}."; exit 1; }
done

case "${ALIGNER}" in
    bwa-mem2)
        [ -x "${E}/bwa-mem2" ] || { echo "Error: ${PREPROCESS_HIC_ENV} has no bwa-mem2."; exit 1; }
        case "${BWA_INDEX}" in *.tar.gz|*.tgz) ;; *)
            for e in 0123 bwt.2bit.64; do
                [ -s "${BWA_INDEX}.${e}" ] || { echo "Error: ${BWA_INDEX}.${e} not found. Build it (bwa-mem2 index ${BWA_INDEX}) or set ALIGNER=bwa."; exit 1; }
            done ;;
        esac ;;
    bwa) ;;
    *) echo "Error: ALIGNER must be bwa-mem2 or bwa (got '${ALIGNER}')."; exit 1;;
esac
if [[ ! "${MCOOL_MIN_MAPQ}" =~ ^[0-9]+$ ]]; then echo "Error: MCOOL_MIN_MAPQ must be an integer."; exit 1; fi
if [[ ! "${MIN_MAPQ}" =~ ^[0-9]+$ ]]; then echo "Error: MIN_MAPQ must be an integer."; exit 1; fi

abspath() { case "$1" in /*) echo "$1" ;; *) echo "${PWD}/$1" ;; esac; }
WORKDIR="${NXF_WORK:-${OUTDIR}/work}"
mkdir -p "${OUTDIR}" "${TMPDIR_NF}" "${WORKDIR}"
OUTDIR=$(abspath "${OUTDIR}"); TMPDIR_NF=$(abspath "${TMPDIR_NF}"); WORKDIR=$(abspath "${WORKDIR}")
BWA_INDEX=$(abspath "${BWA_INDEX}"); CHROM_SIZES=$(abspath "${CHROM_SIZES}"); FASTQ_GLOB=$(abspath "${FASTQ_GLOB}")

FRAG_OPTS=""
if [ -n "${FRAG_BED}" ]; then
    FRAG_OPTS="--frag_bed '$(abspath "${FRAG_BED}")'"
else
    if [ -z "${GENOME_FASTA}" ]; then
        case "${BWA_INDEX}" in *.tar.gz|*.tgz) echo "Error: BWA_INDEX is a tarball: set GENOME_FASTA or FRAG_BED."; exit 1;; esac
        GENOME_FASTA="${BWA_INDEX}"
    fi
    [ -s "${GENOME_FASTA}" ] || { echo "Error: GENOME_FASTA ${GENOME_FASTA} not found."; exit 1; }
    FRAG_OPTS="--genome_fasta '$(abspath "${GENOME_FASTA}")' --enzymes '${ENZYMES}'"
fi

DS_OPTS=""
if [ -n "${DATASET_DIR}" ]; then
    DATASET_DIR=$(abspath "${DATASET_DIR}")
    mkdir -p "${DATASET_DIR}/pairs" "${DATASET_DIR}/bam" "${DATASET_DIR}/qc"
    DS_OPTS="--pairsdir '${DATASET_DIR}/pairs' --bamdir '${DATASET_DIR}/bam' --statsdir '${DATASET_DIR}/qc' --coolsdir '${DATASET_DIR}/pairs'"
fi
ASM_OPT=""
[ -n "${ASSEMBLY}" ] && ASM_OPT="--assembly '${ASSEMBLY}'"

GIT_COMMIT=$(cd "${REPO}" && git rev-parse HEAD 2>/dev/null || echo unknown)
if [ -n "$(cd "${REPO}" && git status --porcelain -- nf/preprocess-hic preprocess-hic-submit.sh 2>/dev/null)" ]; then
    GIT_COMMIT="${GIT_COMMIT}-dirty"
fi

FAILED_NODES_FILE="${OUTDIR}/failed_nodes_hic_${SAMPLEID}.txt"
touch "${FAILED_NODES_FILE}"
EXCLUDE_NODES=$(sort -u "${FAILED_NODES_FILE}" | tr '\n' ',' | sed 's/,$//')
EXCLUDE_OPT=""
[ -n "${EXCLUDE_NODES}" ] && EXCLUDE_OPT="--exclude=${EXCLUDE_NODES}"

# Nextflow runs in OUTDIR: per-sample cache, .nextflow.log, master log and trace there; the checkout is never written to.
sbatch --mem=4G ${EXCLUDE_OPT} -J "hic-${SAMPLEID}" -o "${OUTDIR}/${SAMPLEID}.master.%j.log" --wrap="
  export PREPROCESS_HIC_ENV=${PREPROCESS_HIC_ENV}
  export JAVA_HOME=${PREPROCESS_HIC_ENV}
  export JAVA_CMD=${PREPROCESS_HIC_ENV}/bin/java
  cd '${OUTDIR}'
  ${E}/nextflow run -resume \
    -w '${WORKDIR}' \
    -c '${REPO}/nf/preprocess-hic/nextflow.config' \
    '${REPO}/nf/preprocess-hic/main.nf' \
    --sample_id '${SAMPLEID}' \
    --bwa_index '${BWA_INDEX}' \
    --chrom_sizes '${CHROM_SIZES}' \
    --fq_input '${FASTQ_GLOB}' \
    --outdir '${OUTDIR}' \
    --tmpdir '${TMPDIR_NF}' \
    --failed_nodes_file '${FAILED_NODES_FILE}' \
    --aligner ${ALIGNER} \
    --chunk_reads ${CHUNK_READS} \
    --min_mapq ${MIN_MAPQ} \
    --mcool_min_mapq ${MCOOL_MIN_MAPQ} \
    --mcool_resolutions '${MCOOL_RESOLUTIONS}' \
    --git_commit '${GIT_COMMIT}' \
    ${FRAG_OPTS} ${ASM_OPT} ${DS_OPTS}
"
