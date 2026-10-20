#!/usr/bin/env bash
# Portable SLURM launcher: no Nextflow or MPIfR partition assumptions.
# Example: sbatch --array=1-9 scripts/slurm_inception_lite.sh /path/config.json host
# Set INCEPTION_SIF=/path/image.sif and second arg sif for SingularityCE.
#SBATCH --job-name=inception-lite
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=inception_%A_%a.log

set -euo pipefail
if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch --array=1-N $0 CONFIG_JSON [host|sif] [yes|no]" >&2
  exit 2
fi
CONFIG="$(readlink -f "$1")"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"
ENV_MODE="$( [ $# -ge 2 ] && echo "$2" || echo host )"
FOLD_FLAGS=()
if [[ $# -ge 3 ]]; then
  case "$3" in
    yes|no) FOLD_FLAGS=(--perform-folds "$3") ;;
    *) echo "Third argument must be yes or no" >&2; exit 2 ;;
  esac
fi
ARRAY_INDEX="$(printenv SLURM_ARRAY_TASK_ID || echo 1)"
printf -v CASE "case_%04d" "$ARRAY_INDEX"
echo "Running $CASE with $CONFIG ($ENV_MODE)"

if [[ "$ENV_MODE" == "host" ]]; then
  python3 -m inception_lite --config "$CONFIG" --case "$CASE" "${FOLD_FLAGS[@]}"
elif [[ "$ENV_MODE" == "sif" ]]; then
  if [[ -z "$(printenv INCEPTION_SIF 2>/dev/null || true)" ]]; then
    echo "Set INCEPTION_SIF to a real .sif path" >&2
    exit 2
  fi
  BIND_ARGS=()
  if [[ -n "$(printenv INCEPTION_BIND 2>/dev/null || true)" ]]; then
    BIND_ARGS=(--bind "$INCEPTION_BIND")
  fi
  singularity exec --cleanenv "${BIND_ARGS[@]}" "$INCEPTION_SIF" \
    python -m inception_lite --config "$CONFIG" --case "$CASE" "${FOLD_FLAGS[@]}"
else
  echo "Second argument must be host or sif" >&2
  exit 2
fi
