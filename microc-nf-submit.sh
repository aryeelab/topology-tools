#!/bin/bash
# microc-nf-submit.sh — submit the Micro-C Nextflow pipeline on the Aryee lab Slurm cluster
#
# Edit SAMPLEID before running. All other paths should work as-is for the
# Johnstone lab 2022 LPS Micro-C dataset.
#
# Usage: bash microc-nf-submit.sh

REPO=/aryeelab/users/martin/projects/topology-tools

# --- Edit this ---
export SAMPLEID=221014_LPS141_MicroC_RepB
# -----------------

export SINGULARITY_CACHE_DIR=/cluster/aryeelab/mark/singularityimages/
export JAVA_HOME=/homes9/martin/miniforge3/envs/microc
export JAVA_CMD=/homes9/martin/miniforge3/envs/microc/bin/java

# Per-sample node blacklist — only written on job failure (ERR trap in pipeline)
FAILED_NODES_FILE=/cluster/aryeelab/martin/failed_nodes_microc_${SAMPLEID}.txt
touch ${FAILED_NODES_FILE}

# Build --exclude flag for the master job from the same exclusion list
EXCLUDE_NODES=$(sort -u ${FAILED_NODES_FILE} | tr '\n' ',' | sed 's/,$//')
EXCLUDE_OPT=""
if [ -n "${EXCLUDE_NODES}" ]; then
    EXCLUDE_OPT="--exclude=${EXCLUDE_NODES}"
fi

sbatch --mem=4G ${EXCLUDE_OPT} --wrap="
  export SAMPLEID=${SAMPLEID}
  export JAVA_HOME=/homes9/martin/miniforge3/envs/microc
  export JAVA_CMD=/homes9/martin/miniforge3/envs/microc/bin/java
  cd ${REPO}
  /homes9/martin/miniforge3/envs/microc/bin/nextflow run -resume \
    -c nf/nextflow.config \
    nf/microc_preprocess.nf
"
