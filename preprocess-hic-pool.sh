#!/bin/bash
# preprocess-hic-pool.sh — pool replicate .mcool files of preprocess-hic into one balanced .mcool (Slurm job)
#
# Usage: POOL_ID=KT OUTDIR=/path/out bash preprocess-hic-pool.sh rep1.mapq_30.mcool rep2.mapq_30.mcool [...]
# Writes OUTDIR/<POOL_ID>.mapq_<Q>.mcool and its .pool.provenance.json (see nf/preprocess-hic/bin/pool_mcools.py).
set -euo pipefail

POOL_ID="${POOL_ID:-}"
OUTDIR="${OUTDIR:-}"
NPROC="${NPROC:-8}"
MEM="${MEM:-32G}"
PIPELINE_VERSION=1.1.1
PREPROCESS_HIC_ENV=${PREPROCESS_HIC_ENV:-${HOME}/miniforge3/envs/preprocess-hic-${PIPELINE_VERSION}}

[ -n "${POOL_ID}" ] && [ -n "${OUTDIR}" ] || { echo "Error: set POOL_ID and OUTDIR."; exit 1; }
[ $# -ge 2 ] || { echo "Error: give at least two .mcool files."; exit 1; }
[ -x "${PREPROCESS_HIC_ENV}/bin/cooler" ] || { echo "Error: conda env ${PREPROCESS_HIC_ENV} not found."; exit 1; }
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
abspath() { case "$1" in /*) echo "$1" ;; *) echo "${PWD}/$1" ;; esac; }
OUTDIR=$(abspath "${OUTDIR}")
IN=""
for m in "$@"; do [ -s "$m" ] || { echo "Error: $m not found."; exit 1; }; IN="${IN} '$(abspath "$m")'"; done
GIT_COMMIT=$(cd "${REPO}" && git rev-parse HEAD 2>/dev/null || echo unknown)
[ -n "$(cd "${REPO}" && git status --porcelain -- nf/preprocess-hic 2>/dev/null)" ] && GIT_COMMIT="${GIT_COMMIT}-dirty"
mkdir -p "${OUTDIR}"
# The helper is copied into OUTDIR, so the job never reads a script that a later checkout could change.
cp "${REPO}/nf/preprocess-hic/bin/pool_mcools.py" "${OUTDIR}/.pool_mcools.${POOL_ID}.py"

sbatch -c "${NPROC}" --mem="${MEM}" -J "hicpool-${POOL_ID}" -o "${OUTDIR}/${POOL_ID}.pool.%j.log" --wrap="
  export PATH=${PREPROCESS_HIC_ENV}/bin:\$PATH
  python '${OUTDIR}/.pool_mcools.${POOL_ID}.py' --pool-id '${POOL_ID}' --outdir '${OUTDIR}' --nproc ${NPROC} \
    --git-commit '${GIT_COMMIT}' --version '${PIPELINE_VERSION}' ${IN}
"
