#!/usr/bin/env bash
# Submit one window per Slurm job using:
#   sbatch --cpus-per-task=16 --output=/.../W1_slurm_%j.out \
#          scripts/slurm_m62_weekend_single_window.sh 1
# Choose partition, node and time limit on the sbatch command after
# inspecting the cluster scheduler. This script never deletes old output.
set -Eeuo pipefail

WINDOW="${1:?Usage: slurm_m62_weekend_single_window.sh WINDOW_ID(1..4)}"
case "$WINDOW" in 1|2|3|4) ;; *) echo "Invalid window: $WINDOW" >&2; exit 2;; esac

WORK="$HOME/PULSEJET_validation/stage6_m62"
SIF="$HOME/PULSEJET_validation/images/inception_lite.sif"
CFG="$WORK/weekend_configs/m62_weekend_20_w${WINDOW}.json"
INPUT="$WORK/m62_15m_windows/W${WINDOW}_M62_background.fil"
RESULTS="$WORK/m62_15m_injections/W${WINDOW}"
START=$(date +%s)

on_exit() {
  code=$?
  echo "=== WINDOW W${WINDOW} FINAL exit=${code} at $(date -Is); elapsed=$(( $(date +%s)-START )) seconds ==="
}
trap on_exit EXIT

echo "=== Starting W${WINDOW} on host $(hostname) ==="
echo "Slurm job ID: ${SLURM_JOB_ID:-none}; cpus requested: ${SLURM_CPUS_PER_TASK:-unknown}"
date -Is

for path in "$SIF" "$CFG" "$INPUT"; do
  test -s "$path" || { echo "ERROR: input missing or empty: $path" >&2; exit 3; }
done

if [ -e "$RESULTS" ]; then
  echo "ERROR: Refusing to overwrite existing output directory: $RESULTS" >&2
  exit 4
fi

# Requirement: match INCEPTION-Lite config ncpu=16 with Slurm allocation.
if [ "${SLURM_CPUS_PER_TASK:-0}" -lt 16 ]; then
  echo "ERROR: Request >=16 Slurm CPUs per task to match injection config." >&2
  exit 5
fi

# The baseline input is read-only; injector writes a new output file and
# auto-folds five truth models via the real PulsarX executable.
# PYTHONNOUSERSITE prevents conflicting ~/.local pandas/NumPy leaking in.
singularity exec "$SIF" env PYTHONNOUSERSITE=1 \
  python3 -m inception_lite \
    --config "$CFG" \
    --perform-folds yes

echo "=== W${WINDOW} output products ==="
find "$RESULTS" -maxdepth 4 -type f -printf '%p (%s bytes)\n'
