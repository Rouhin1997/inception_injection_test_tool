#!/usr/bin/env bash
# Stage 10 one five-source binary injection per Slurm array task.
# NO overwrite; original background windows are read-only.
set -Eeuo pipefail
ROOT="${PULSEJET_VALIDATION_ROOT:-/home/rnag/PULSEJET_validation}"
PLAN="$ROOT/stage10_ucb100"
SIF="${INCEPTION_SIF:-$ROOT/images/inception_lite.sif}"
BATCH="${SLURM_ARRAY_TASK_ID:-${1:-}}"
if [[ ! "$BATCH" =~ ^[0-9]+$ ]] || (( 10#$BATCH < 1 || 10#$BATCH > 20 )); then
  echo "Expected Slurm array task 1..20" >&2; exit 2
fi
printf -v NUM "%02d" "$((10#$BATCH))"
CFG="$PLAN/configs/batch_$NUM.json"
OUT="$PLAN/injections/batch_$NUM/case_0001"
META="$PLAN/logs/batch_$NUM.provenance.txt"
START=$(date +%s)
done_report() {
  rc=$?
  echo "=== UCB100 batch $NUM exit=$rc at $(date -Is); elapsed=$(( $(date +%s)-START ))s ==="
}
trap done_report EXIT

[[ -s "$SIF" && -s "$CFG" ]] || { echo "Missing injection SIF/config" >&2; exit 3; }
[[ ! -e "$OUT" && ! -e "$META" ]] || { echo "Output/provenance already exists; refusing overwrite: $OUT $META" >&2; exit 3; }
command -v singularity >/dev/null || { echo "Singularity unavailable" >&2; exit 3; }
CPUS="${SLURM_CPUS_PER_TASK:-0}"
[[ -n "${SLURM_JOB_ID:-}" && "$CPUS" -ge 48 ]] || {
  echo "Require an allocated Slurm job with >=48 CPUs; found $CPUS" >&2; exit 3;
}
# Basic independent validation before starting hours of injection.
python3 - "$CFG" "$PLAN" <<'PY'
import json,sys,os
from pathlib import Path
cfg=json.loads(Path(sys.argv[1]).read_text())
assert cfg["mode"]=="together" and len(cfg["pulsars"])==5
assert cfg["fold"]["enabled"] is False
assert cfg["ncpu"] <= int(os.environ.get("SLURM_CPUS_PER_TASK", "0"))
assert len({p["ID"] for p in cfg["pulsars"]}) == 5
for p in cfg["pulsars"]:
    assert "A0" not in p and "A1" not in p
    if "binary_period" in p:
        assert p["binary_period"]>0 and p["T0"]>50000
        assert all(k in p for k in ("M1","M2","inc","ecc"))
source=Path(cfg["filterbank"])
assert source.is_file()
size=source.stat().st_size
vfs=os.statvfs(sys.argv[2])
free=vfs.f_bavail*vfs.f_frsize
assert free > int(size*1.5), "Not enough free space for new filterbank + temporary files"
assert not (Path(cfg["output_dir"])/"case_0001").exists()
print("Preflight OK; source bytes",size,"free bytes",free,flush=True)
PY
mkdir -p "$PLAN/logs"
{
  printf 'UTC launch: '; date -u '+%Y-%m-%dT%H:%M:%SZ'
  echo "batch=$NUM node=$(hostname) slurm_job=${SLURM_JOB_ID} allocated_cpus=$CPUS"
  echo "injection_workers=48 config=$CFG"
  sha256sum "$SIF" "$CFG"
  echo "Model: full Keplerian binary where supplied; NO polynomial A0/A1 concurrently"
  echo "PulsarX automatic folds OFF: unsupported for full binary model"
} > "$META"

echo "=== Stage10 UCB100 batch $NUM on $(hostname) ==="
echo "Config: $CFG"
echo "Output dir: $PLAN/injections/batch_$NUM"
if [[ -x /usr/bin/time ]]; then
  /usr/bin/time -v singularity exec --cleanenv "$SIF" env PYTHONNOUSERSITE=1 \
    python3 -m inception_lite --config "$CFG" --perform-folds no
else
  singularity exec --cleanenv "$SIF" env PYTHONNOUSERSITE=1 \
    python3 -m inception_lite --config "$CFG" --perform-folds no
fi

python3 - "$CFG" <<'PY'
import json,sys
from pathlib import Path
cfg=json.loads(Path(sys.argv[1]).read_text())
case=Path(cfg["output_dir"])/"case_0001"
manifest=case/"run_manifest.json"
assert manifest.is_file(), "Missing completed run_manifest.json"
run=json.loads(manifest.read_text())
assert run["status"]=="success" and len(run["pulsars"])==5
fil=Path(run["output_fil"])
assert fil.is_file() and fil.stat().st_size>0
pars=list(case.glob("*.par"))
assert len(pars)==5, "Expected 5 orbital/isolated ephemeris .par files, got %d" % len(pars)
assert len(list(case.glob("report_*.json")))==1
assert not run["folding"]
print("PASS: five injections, five .par files, output .fil, manifest and report",fil,flush=True)
PY
echo "completed_utc=$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$META"
echo "SUCCESS UCB100 batch $NUM at $(date -Is)"
