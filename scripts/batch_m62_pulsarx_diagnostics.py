#!/usr/bin/env python3
"""Generate and audit all 20 existing M62 INCEPTION-Lite PulsarX truth folds.

Reads existing .ar files and injection manifests only. Does not fold, inject,
modify source files, or assess PULSEJET detection performance.
Writes corrected three-panel PNGs and a CSV of descriptive profile contrasts.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from inception_lite.fold_plots import read_dedispersed_pulsarx_archive
from scripts.plot_pulsarx_diagnostic import plot


def peak_over_mad(cube, weights):
    """Descriptive profile contrast; NOT calibrated signal-to-noise."""
    total_weight = np.sum(weights)
    if total_weight <= 0:
        raise ValueError("No positive archive weights")
    p = np.einsum("tcb,tc->b", cube, weights) / total_weight
    p = p - np.median(p)
    scatter = 1.4826 * np.median(np.abs(p))
    if not np.isfinite(scatter) or scatter <= 0:
        raise ValueError("Zero/non-finite phase scatter in integrated profile")
    return float(np.max(p) / scatter), int(np.argmax(p))


def original_calibrated_cube(archive):
    """Read the unaligned PulsarX archive for a before/after comparison."""
    with fits.open(archive, memmap=False) as hdul:
        sub = hdul["SUBINT"]
        data = np.asarray(sub.data["DATA"], dtype=np.float64)
        if data.ndim != 4 or data.shape[1] != 1:
            raise ValueError("Expected PulsarX total-intensity SUBINT data")
        nt, _, nc, _ = data.shape
        scl = np.asarray(sub.data["DAT_SCL"], dtype=float).reshape(nt, nc)
        offs = np.asarray(sub.data["DAT_OFFS"], dtype=float).reshape(nt, nc)
        return data[:, 0] * scl[:, :, None] + offs[:, :, None]


FIELDS = [
    "id", "window", "status", "error", "dm_truth", "dm_archive",
    "period_ms", "acc_m_s2", "jerk_m_s3", "snr_target",
    "raw_peak_mad", "aligned_peak_mad", "gain",
    "aligned_peak_phase", "nsub", "nchan", "nbin", "valid_weight_fraction", "png"
]


def process_one(case, outdir, pulse):
    source_id = pulse["ID"]
    fold_dir = case / "folds" / source_id
    files = sorted(fold_dir.glob("*.ar"))
    if len(files) != 1:
        raise ValueError(f"Expected exactly one .ar in {fold_dir}, found {len(files)}")
    archive = files[0]

    folded = read_dedispersed_pulsarx_archive(archive)
    dm_truth = float(pulse.get("DM", 0))
    dm_archive = folded["dm"]
    if not np.isclose(dm_archive, dm_truth, rtol=0, atol=1e-5):
        raise ValueError(
            f"Archive DM {dm_archive} differs from injection manifest DM {dm_truth}"
        )

    weights = folded["weights"]
    raw = original_calibrated_cube(archive)
    aligned = folded["intensity"]
    before, _ = peak_over_mad(raw, weights)
    after, peak_bin = peak_over_mad(aligned, weights)
    nt, nc, nb = aligned.shape

    png = outdir / (source_id + ".png")
    plot(archive, png, source_id, case / "run_manifest.json")

    period = float(pulse.get("P0", 1 / float(pulse["F0"]))) if "P0" not in pulse else float(pulse["P0"])
    return {
        "id": source_id, "window": source_id.split("_")[1],
        "status": "PLOTTED", "error": "",
        "dm_truth": dm_truth, "dm_archive": dm_archive,
        "period_ms": 1e3 * period,
        "acc_m_s2": float(pulse.get("A0", 0)),
        "jerk_m_s3": float(pulse.get("A1", 0)),
        "snr_target": float(pulse["SNR"]),
        "raw_peak_mad": round(before, 3),
        "aligned_peak_mad": round(after, 3),
        "gain": round(after / before, 3) if before > 0 else "",
        "aligned_peak_phase": round(peak_bin / nb, 4),
        "nsub": nt, "nchan": nc, "nbin": nb,
        "valid_weight_fraction": round(np.count_nonzero(weights) / weights.size, 4),
        "png": str(png)
    }


def run(root, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    results = []
    nplanned = 0
    for i in range(1, 5):
        case = root / f"W{i}" / "case_0001"
        manifest_path = case / "run_manifest.json"
        if not manifest_path.is_file():
            raise FileNotFoundError(f"Missing run manifest: {manifest_path}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        pulsars = manifest.get("pulsars", [])
        if len(pulsars) != 5:
            raise ValueError(f"Expected 5 injected pulsars in {manifest_path}; found {len(pulsars)}")
        for pulse in pulsars:
            source_id = pulse["ID"]
            nplanned += 1
            try:
                result = process_one(case, outdir, pulse)
                print(
                    f"{source_id}: raw={result['raw_peak_mad']:.2f} "
                    f"aligned={result['aligned_peak_mad']:.2f} "
                    f"gain={result['gain']:.2f}x"
                )
            except Exception as exc:
                result = {"id": source_id, "window": f"W{i}",
                          "status": "ERROR", "error": str(exc)}
                print(f"{source_id}: ERROR: {exc}", file=sys.stderr)
            results.append(result)

    report = outdir / "m62_20_truth_fold_qc.csv"
    with report.open("w", newline="", encoding="utf-8") as fd:
        writer = csv.DictWriter(fd, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    plotted = sum(r["status"] == "PLOTTED" for r in results)
    print(f"Summary: {plotted}/{nplanned} archives plotted; report: {report}")
    print("Peak/MAD and gain are descriptive profile contrasts, NOT calibrated S/N.")
    print("Every plot still requires visual inspection; weak injections may not give large contrasts.")
    if nplanned != 20 or plotted != 20:
        return 1
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True,
                        help="Directory holding W1/case_0001 through W4/case_0001")
    parser.add_argument("--out", type=Path, default=None,
                        help="Fresh output folder for PNGs/CSV; defaults to ROOT/pulsarx_diagnostics_corrected")
    args = parser.parse_args()
    outdir = args.out or args.root / "pulsarx_diagnostics_corrected"
    raise SystemExit(run(args.root, outdir))
