#!/bin/bash
# microc-nf-submit.sh — submit the Micro-C Nextflow pipeline on a Slurm cluster
#
# Requirements: the 'microc' conda environment must be installed.
#   conda env create -f microc-env.yml
#
# Usage: bash microc-nf-submit.sh

set -euo pipefail

# ============================================================
# User configuration — edit these for your site/sample
# ============================================================
SAMPLEID=221014_LPS141_MicroC_RepB

# Reference genome files
BWA_INDEX=/aryeelab/users/mark/nextflowmicroc/hg38.tgz
CHROM_SIZES=/aryeelab/users/mark/nextflowmicroc/hg38.chrom.sizes

# Input FASTQs — glob pattern passed to Nextflow's fromFilePairs
FASTQ_GLOB="/aryeelab/data/johnstonelab/2022_LPS_MicroC/fastq/${SAMPLEID}/*R{1,2}*.fastq.gz"

# Output directory and temp directory (both must be on a shared filesystem)
OUTDIR=/cluster/aryeelab/${USER}/${SAMPLEID}
TMPDIR_NF=/cluster/aryeelab/${USER}/tmp
# ============================================================

# Resolve the repo root from the script's location
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Auto-detect the microc conda environment
export MICROC_ENV
MICROC_ENV=$(conda env list 2>/dev/null | awk '$1=="microc"{print $NF}')
if [ -z "${MICROC_ENV}" ]; then
    echo "Error: 'microc' conda environment not found."
    echo "Create it with: conda env create -f ${REPO}/microc-env.yml"
    exit 1
fi

export JAVA_HOME="${MICROC_ENV}"
export JAVA_CMD="${MICROC_ENV}/bin/java"
NEXTFLOW="${MICROC_ENV}/bin/nextflow"

mkdir -p "${OUTDIR}" "${TMPDIR_NF}"

# Per-sample node blacklist — written only on job failure (ERR trap in pipeline).
# Persists across re-submissions so bad nodes stay excluded until manually cleared.
FAILED_NODES_FILE="${OUTDIR}/failed_nodes_microc_${SAMPLEID}.txt"
touch "${FAILED_NODES_FILE}"

# Build --exclude flag for the master job from the same blacklist
EXCLUDE_NODES=$(sort -u "${FAILED_NODES_FILE}" | tr '\n' ',' | sed 's/,$//')
EXCLUDE_OPT=""
[ -n "${EXCLUDE_NODES}" ] && EXCLUDE_OPT="--exclude=${EXCLUDE_NODES}"

sbatch --mem=4G ${EXCLUDE_OPT} --wrap="
  export MICROC_ENV=${MICROC_ENV}
  export JAVA_HOME=${MICROC_ENV}
  export JAVA_CMD=${MICROC_ENV}/bin/java
  cd ${REPO}
  ${NEXTFLOW} run -resume \
    -c nf/nextflow.config \
    nf/microc_preprocess.nf \
    --sample_id '${SAMPLEID}' \
    --bwa_index '${BWA_INDEX}' \
    --chrom_sizes '${CHROM_SIZES}' \
    --fq_input '${FASTQ_GLOB}' \
    --outdir '${OUTDIR}' \
    --tmpdir '${TMPDIR_NF}' \
    --failed_nodes_file '${FAILED_NODES_FILE}'
"
