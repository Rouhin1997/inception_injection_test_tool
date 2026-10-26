#!/usr/bin/env python3
"""Three-panel scientific diagnostic from an actual PulsarX folded PSRFITS .ar.

Uses the full SUBINT DATA cube, DAT_SCL, DAT_OFFS, DAT_WTS and DAT_FREQ.
No injected model/profile is synthesized; the fold itself must already exist.
Each time/channel cell is baseline-subtracted (median over phase) and divided
by its robust phase noise for display, so plotted values are relative S/N-like
intensities, NOT flux or calibrated pulsar detection significance.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def weighted_mean_profiles(cube, weights, axis):
    """Mean across one or both (subint, frequency) dimensions."""
    if axis == "frequency":
        denominator = weights.sum(axis=1)[:, None]
        numerator = np.einsum("tc,tcb->tb", weights, cube)
    elif axis == "time":
        denominator = weights.sum(axis=0)[:, None]
        numerator = np.einsum("tc,tcb->cb", weights, cube)
    elif axis == "all":
        denominator = weights.sum()
        numerator = np.einsum("tc,tcb->b", weights, cube)
    else:
        raise ValueError("Invalid reduction")
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)


def read_injection_truth(archive: Path, source_id: str, manifest_path: Path | None = None):
    """Load *injected* parameters from the saved run manifest, never fit data.

    An INCEPTION-Lite archive is case_0001/folds/PSR_ID/PSR_ID*.ar;
    the original injection run manifest lives three parents above it.
    """
    case_dir = archive.parent.parent.parent
    path = manifest_path if manifest_path is not None else case_dir / "run_manifest.json"
    if not path.is_file():
        if manifest_path is not None:
            raise FileNotFoundError(f"Requested injection manifest not found: {path}")
        return None
    obj = json.loads(path.read_text(encoding="utf-8"))
    matches = [p for p in obj.get("pulsars", []) if p.get("ID") == source_id]
    if len(matches) != 1:
        raise ValueError(
            f"Manifest {path} contains {len(matches)} source(s) named {source_id!r}; "
            "refusing to label plot with unrelated truth parameters"
        )
    pulse = matches[0]
    resolved = case_dir / "resolved_injection.json"
    global_config = {}
    if resolved.is_file():
        global_config = json.loads(resolved.read_text(encoding="utf-8")).get("psr_global", {})
    if "P0" in pulse:
        period = float(pulse["P0"])
        frequency = 1.0 / period
    else:
        frequency = float(pulse["F0"])
        period = 1.0 / frequency
    if not period > 0:
        raise ValueError("Injected period/frequency is not positive")
    get = lambda key, default=None: pulse.get(key, global_config.get(key, default))
    return {
        "id": source_id,
        "period": period,
        "frequency": frequency,
        "dm": float(get("DM", 0)),
        "a": float(get("A0", 0)),
        "jerk": float(get("A1", 0)),
        "snr": float(pulse["SNR"]),
        "duty": float(get("duty_cycle", 0)),
        "frame": str(get("frame", "unspecified")),
        "accep": get("ACCEPOCH", None),
        "manifest": str(path),
    }


def plot(archive: Path, output: Path, label: str, manifest: Path | None = None):
    truth = read_injection_truth(archive, archive.parent.name, manifest)
    with fits.open(archive, memmap=False) as hdul:
        sub = hdul["SUBINT"]
        data = np.asarray(sub.data["DATA"], dtype=np.float64)
        if data.ndim != 4:
            raise ValueError(f"Expected SUBINT DATA shape [time,pol,chan,phase], got {data.shape}")
        nsub, npol, nchan, nbin = data.shape
        if npol != 1:
            raise ValueError(f"Expected total-intensity 1-pol archive, received npol={npol}")
        if nsub < 2 or nchan < 2 or nbin < 8:
            raise ValueError("Not enough samples/channels/bins for a 3-panel diagnostic")

        # In PSRFITS, integer-quantized DATA are converted to real intensity
        # using DATA * DAT_SCL + DAT_OFFS for each subint/polarization/channel.
        scales = np.asarray(sub.data["DAT_SCL"], dtype=np.float64).reshape(nsub, npol, nchan)
        offsets = np.asarray(sub.data["DAT_OFFS"], dtype=np.float64).reshape(nsub, npol, nchan)
        weights = np.asarray(sub.data["DAT_WTS"], dtype=np.float64).reshape(nsub, nchan)
        frequencies = np.asarray(sub.data["DAT_FREQ"], dtype=np.float64).reshape(nsub, nchan)
        tsub = np.asarray(sub.data["TSUBINT"], dtype=np.float64).ravel()
        tmid = np.asarray(sub.data["OFFS_SUB"], dtype=np.float64).ravel()

        intensity = data[:, 0, :, :] * scales[:, 0, :, None] + offsets[:, 0, :, None]

    # Median phase baseline and robust off-pulse-like phase scatter for each
    # time/channel cell; source duty cycle is expected to be small.
    baseline = np.nanmedian(intensity, axis=-1, keepdims=True)
    residual = intensity - baseline
    robust_sigma = 1.4826 * np.nanmedian(np.abs(residual), axis=-1, keepdims=True)
    fallback_sigma = np.nanstd(residual, axis=-1, keepdims=True)
    sigma = np.where(robust_sigma > 1.0e-12, robust_sigma, fallback_sigma)

    valid = np.isfinite(weights) & (weights > 0) & np.isfinite(sigma[:, :, 0]) & (sigma[:, :, 0] > 1e-12)
    weights = np.where(valid, weights, 0.0)
    if weights.sum() <= 0:
        raise ValueError("No unmasked nonconstant folded time/channel cells found")
    normalized = np.divide(residual, sigma, out=np.zeros_like(residual), where=sigma > 1e-12)
    normalized = np.nan_to_num(normalized, nan=0.0, posinf=0.0, neginf=0.0)

    profile = weighted_mean_profiles(normalized, weights, "all")
    time_phase = weighted_mean_profiles(normalized, weights, "frequency")
    freq_phase = weighted_mean_profiles(normalized, weights, "time")

    if not np.all(np.isfinite(profile)) or np.ptp(profile) <= 0:
        raise ValueError("Folded profile is not finite or shows no phase variation")

    # Sort frequency channels in ascending MHz regardless of foff sign.
    freq_mhz = np.nanmedian(frequencies, axis=0)
    order = np.argsort(freq_mhz)
    freq_mhz = freq_mhz[order]
    freq_phase = freq_phase[order, :]
    dfreq = np.median(np.diff(freq_mhz))
    freq_lo = freq_mhz[0] - dfreq / 2.0
    freq_hi = freq_mhz[-1] + dfreq / 2.0

    # SUBINT OFFS_SUB is the subintegration midpoint relative to observation
    # reference time; plot elapsed seconds from the beginning of this fold.
    if np.all(np.isfinite(tmid)) and np.all(np.diff(tmid) > 0):
        first = tmid[0] - tsub[0] / 2.0
        last = tmid[-1] + tsub[-1] / 2.0
        t_lo, t_hi = 0.0, last - first
    else:
        t_lo, t_hi = 0.0, np.sum(tsub)

    two_cycle_phase = np.arange(2 * nbin) / nbin
    doubled_profile = np.tile(profile, 2)
    doubled_time = np.tile(time_phase, (1, 2))
    doubled_frequency = np.tile(freq_phase, (1, 2))

    plt.rcParams.update({"font.size": 10, "figure.dpi": 120, "savefig.dpi": 170})
    fig = plt.figure(figsize=(12, 11.3))
    gs = GridSpec(4, 1, figure=fig, height_ratios=[0.89, 1.05, 1.55, 1.55], hspace=0.33)
    axmeta = fig.add_subplot(gs[0])
    ax0 = fig.add_subplot(gs[1])
    ax1 = fig.add_subplot(gs[2], sharex=ax0)
    ax2 = fig.add_subplot(gs[3], sharex=ax0)

    # Metadata are injected *truth* taken from run_manifest.json, not
    # recovered values or an S/N estimated from the folded intensity.
    axmeta.set_xlim(0, 1)
    axmeta.set_ylim(0, 1)
    axmeta.set_facecolor("#f6f8fb")
    axmeta.set_xticks([])
    axmeta.set_yticks([])
    for spine in axmeta.spines.values():
        spine.set_color("#d3dce7")
    if truth is not None:
        epoch_text = (
            f'{truth["frame"]}; ACCEPOCH={float(truth["accep"]):g}'
            if truth["accep"] is not None
            else truth["frame"]
        )
        entries = [
            ("Spin period P₀", f'{truth["period"] * 1e3:.9f} ms'),
            ("Spin frequency f₀", f'{truth["frequency"]:.6f} Hz'),
            ("Additional DM", f'{truth["dm"]:g} pc cm⁻³'),
            ("Acceleration a₀", f'{truth["a"]:+g} m s⁻²'),
            ("Jerk j₀", f'{truth["jerk"]:+g} m s⁻³'),
            ("Injected S/N target", f'{truth["snr"]:g} (not measured)'),
            ("Duty cycle", f'{100 * truth["duty"]:g}%'),
            ("Reference epoch", epoch_text),
        ]
        axmeta.text(
            0.02, 0.92, "INJECTION TRUTH  ·  Not search-recovered measurements",
            ha="left", va="top", fontsize=9.5, fontweight="bold", color="#26394d",
        )
        for idx, (key, value) in enumerate(entries):
            col, row = idx % 4, idx // 4
            x = 0.02 + 0.245 * col
            y = 0.68 - 0.40 * row
            axmeta.text(x, y, key, ha="left", va="center", fontsize=8.5, color="#53677e")
            axmeta.text(x, y - 0.16, value, ha="left", va="center",
                        fontsize=10.1, color="#15263b", fontweight="semibold")
    else:
        axmeta.text(0.02, 0.72, "INJECTION TRUTH NOT AVAILABLE", fontsize=10,
                    color="#7b3f2c", fontweight="bold")
        axmeta.text(0.02, 0.37, "Provide --manifest or use the original INCEPTION-Lite fold directory.",
                    fontsize=9.5, color="#53677e")

    ax0.plot(two_cycle_phase, doubled_profile, color="#283d6a", lw=1.7)
    ax0.set_ylabel("Relative folded intensity")
    ax0.set_title(f"{label}  |  PulsarX known-parameter truth fold", fontsize=13, pad=9)
    ax0.grid(alpha=0.20)
    ax0.set_xlim(0, 2)
    ax0.axvline(1, color="gray", alpha=0.4, lw=0.8)

    def dynamic_limits(array):
        finite = array[np.isfinite(array)]
        if not len(finite):
            raise ValueError("Invalid folded waterfall")
        lo, hi = np.percentile(finite, [1, 99.5])
        if hi <= lo:
            lo, hi = float(np.min(finite)), float(np.max(finite))
        if hi <= lo:
            hi = lo + 1
        return lo, hi

    p1 = ax1.imshow(doubled_time, origin="lower", aspect="auto",
                    extent=[0, 2, t_lo / 60.0, t_hi / 60.0],
                    cmap="magma", vmin=dynamic_limits(doubled_time)[0],
                    vmax=dynamic_limits(doubled_time)[1], interpolation="nearest")
    ax1.set_ylabel("Time from fold start (min)")
    ax1.axvline(1, color="white", alpha=0.25, lw=0.6)
    cb1 = fig.colorbar(p1, ax=ax1, fraction=0.018, pad=0.012)
    cb1.set_label("Relative intensity")

    p2 = ax2.imshow(doubled_frequency, origin="lower", aspect="auto",
                    extent=[0, 2, freq_lo, freq_hi],
                    cmap="magma", vmin=dynamic_limits(doubled_frequency)[0],
                    vmax=dynamic_limits(doubled_frequency)[1], interpolation="nearest")
    ax2.set_ylabel("Observing frequency (MHz)")
    ax2.set_xlabel("Pulse phase (cycles; repeated twice)")
    ax2.axvline(1, color="white", alpha=0.25, lw=0.6)
    cb2 = fig.colorbar(p2, ax=ax2, fraction=0.018, pad=0.012)
    cb2.set_label("Relative intensity")

    fig.text(0.06, 0.008,
             f"Real PulsarX PSRFITS SUBINT data: {nsub} subintegrations × {nchan} channels × {nbin} bins. "
             "DAT_SCL/OFFS/WTS applied; per-cell median baseline and robust-noise normalization. "
             "Not flux calibrated.",
             ha="left", va="bottom", fontsize=8, color="dimgray")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", dpi=170)
    plt.close(fig)
    print(f"Saved diagnostic: {output}")
    print(f"Archive: {archive}")
    print(f"SUBINT: {nsub} x {nchan} x {nbin}; unmasked fraction: {np.count_nonzero(valid)/valid.size:.3f}")
    print(f"Frequency range: {freq_lo:.2f} to {freq_hi:.2f} MHz; folded duration: {t_hi / 60:.3f} min")
    print(f"Peak relative folded intensity: {float(np.max(profile)):.4f} (NOT calibrated S/N)")
    return output


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, required=True, help="PulsarX .ar archive")
    p.add_argument("--output", type=Path, required=True, help="Output 3-panel PNG path")
    p.add_argument("--label", type=str, default=None, help="Plot label (injection truth is always matched by archive parent directory ID)")
    p.add_argument("--manifest", type=Path, default=None,
                   help="Optional injection run_manifest.json; defaults to case_0001/run_manifest.json")
    args = p.parse_args()
    plot(args.archive, args.output, args.label or args.archive.parent.name, args.manifest)
