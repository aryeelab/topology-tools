#!/bin/bash
# microc-nf-submit.sh — submit the Micro-C Nextflow pipeline on a Slurm cluster
#
# Requirements: the microc-preprocess conda env for this release, built from
#   envs/microc-preprocess-<version>.linux-64.explicit.txt (see RELEASING.md).
#
# Usage: bash microc-nf-submit.sh

set -euo pipefail

# ============================================================
# User configuration — edit these for your site/sample, or set them in the environment.
# The first six are required; there are no site defaults.
# ============================================================
SAMPLEID="${SAMPLEID:-}"                 # sample identifier, e.g. 221014_LPS141_MicroC_RepB

# Reference genome files
BWA_INDEX="${BWA_INDEX:-}"               # BWA index: .tgz tarball, or the indexed FASTA (prefix of .bwt etc.)
CHROM_SIZES="${CHROM_SIZES:-}"           # chromosome sizes file for the same genome

# Input FASTQs — glob pattern passed to Nextflow's fromFilePairs, e.g. "/path/to/fastq/${SAMPLEID}/*R{1,2}*.fastq.gz"
# (keep the {1,2} default out of ${FASTQ_GLOB:-...}: its "}" would end the expansion)
FASTQ_GLOB="${FASTQ_GLOB:-}"

# Nextflow run directory and temp directory (both must be on a shared filesystem)
OUTDIR="${OUTDIR:-}"                     # provenance, manifest, failed-node list; also outputs if DATASET_DIR is empty
                                         # Nextflow work dir: ${NXF_WORK:-OUTDIR/work} (large; delete after checking outputs)
TMPDIR_NF="${TMPDIR_NF:-}"               # temp dir for pairtools sort/merge

# Optional: dataset directory (<DS>). Pairs, BAM and dedup stats are copied to <DS>/{pairs,bam,qc}
# and a one-row manifest for the FOCI-3D batch pipeline is written; empty = keep everything in OUTDIR.
DATASET_DIR="${DATASET_DIR:-}"
ASSAY="${ASSAY:-microc}"                 # microc | rcmc (manifest column)

# pairtools parse --min-mapq
MAPQ="${MAPQ:-30}"

# Reads per FASTQ chunk for parallel alignment
CHUNK_READS="${CHUNK_READS:-30000000}"

# Conda env for this release (override with MICROC_ENV=/path/to/env)
PIPELINE_VERSION=1.0.1
MICROC_ENV=${MICROC_ENV:-${HOME}/miniforge3/envs/microc-preprocess-${PIPELINE_VERSION}}
# ============================================================

missing=""
for v in SAMPLEID BWA_INDEX CHROM_SIZES FASTQ_GLOB OUTDIR TMPDIR_NF; do
    [ -n "${!v}" ] || missing="${missing} ${v}"
done
if [ -n "${missing}" ]; then
    echo "Error: set${missing} (edit the User configuration block or export them)."
    exit 1
fi
case "${ASSAY}" in microc|rcmc) ;; *) echo "Error: ASSAY must be microc or rcmc (got '${ASSAY}')."; exit 1;; esac

# Resolve the repo root from the script's location
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

export MICROC_ENV
if [ ! -x "${MICROC_ENV}/bin/pairtools" ]; then
    echo "Error: conda env ${MICROC_ENV} not found."
    echo "Create it with: conda create -p ${MICROC_ENV} --file ${REPO}/envs/microc-preprocess-${PIPELINE_VERSION}.linux-64.explicit.txt"
    exit 1
fi
if ! "${MICROC_ENV}/bin/pairtools" --version | grep -q ' 1\.1\.3$'; then
    echo "Error: ${MICROC_ENV} does not have pairtools 1.1.3 (microc-preprocess ${PIPELINE_VERSION} needs it)."
    exit 1
fi

export JAVA_HOME="${MICROC_ENV}"
export JAVA_CMD="${MICROC_ENV}/bin/java"
NEXTFLOW="${MICROC_ENV}/bin/nextflow"

WORKDIR="${NXF_WORK:-${OUTDIR}/work}"
mkdir -p "${OUTDIR}" "${TMPDIR_NF}" "${WORKDIR}"

GIT_COMMIT=$(cd "${REPO}" && git rev-parse HEAD 2>/dev/null || echo unknown)
DS_OPTS=""
if [ -n "${DATASET_DIR}" ]; then
    mkdir -p "${DATASET_DIR}/pairs" "${DATASET_DIR}/bam" "${DATASET_DIR}/qc"
    DS_OPTS="--dataset_dir '${DATASET_DIR}' --assay '${ASSAY}' --pairsdir '${DATASET_DIR}/pairs' --bamdir '${DATASET_DIR}/bam' --statsdir '${DATASET_DIR}/qc'"
fi

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
    -w '${WORKDIR}' \
    -c nf/nextflow.config \
    nf/microc_preprocess.nf \
    --sample_id '${SAMPLEID}' \
    --bwa_index '${BWA_INDEX}' \
    --chrom_sizes '${CHROM_SIZES}' \
    --fq_input '${FASTQ_GLOB}' \
    --outdir '${OUTDIR}' \
    --tmpdir '${TMPDIR_NF}' \
    --failed_nodes_file '${FAILED_NODES_FILE}' \
    --mapq ${MAPQ} \
    --chunk_reads ${CHUNK_READS} \
    --git_commit '${GIT_COMMIT}' \
    ${DS_OPTS}
"
