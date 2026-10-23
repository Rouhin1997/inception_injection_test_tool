#!/usr/bin/env bash
# Submit one window per Slurm job using:
#   sbatch --cpus-per-task=30 --output=/.../W1_slurm_%j.out \
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

# Use ALL the CPUs actually assigned to this job. Keep source config unchanged
# and save a unique resolved copy for reproducible Monday comparisons.
CPUS="${SLURM_CPUS_PER_TASK:-0}"
if [ "$CPUS" -lt 16 ] || [ "$CPUS" -gt 64 ]; then
  echo "ERROR: Request 16-64 Slurm CPUs per task (received $CPUS)." >&2
  exit 5
fi
JOB_CFG="$WORK/weekend_configs/m62_w${WINDOW}_job${SLURM_JOB_ID:?}_${CPUS}cpu.json"
if [ -e "$JOB_CFG" ]; then
  echo "ERROR: Resolved job config exists: $JOB_CFG" >&2
  exit 6
fi
python3 - "$CFG" "$JOB_CFG" "$CPUS" <<'PY'
import json, sys
from pathlib import Path
source, dest, cpus = sys.argv[1], Path(sys.argv[2]), int(sys.argv[3])
cfg = json.loads(Path(source).read_text())
if cfg.get("mode") != "together" or len(cfg.get("pulsars", [])) != 5:
    raise SystemExit("Unexpected injection configuration; abort")
cfg["ncpu"] = cpus
with dest.open("x") as f:
    json.dump(cfg, f, indent=2)
    f.write("\n")
print("Reproducible job configuration:", dest, "with", cpus, "injection workers", flush=True)
PY

# The baseline input is read-only; injector writes a new output file and
# auto-folds five truth models via the real PulsarX executable.
# PYTHONNOUSERSITE prevents conflicting ~/.local pandas/NumPy leaking in.
singularity exec "$SIF" env PYTHONNOUSERSITE=1 \
  python3 -m inception_lite \
    --config "$JOB_CFG" \
    --perform-folds yes

echo "=== W${WINDOW} output products ==="
find "$RESULTS" -maxdepth 4 -type f -printf '%p (%s bytes)\n'
