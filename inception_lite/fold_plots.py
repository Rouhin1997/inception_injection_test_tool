"""Portable, headless profile diagnostics from a real PulsarX .ar archive.

PulsarX can omit native PNGs in headless containers. Its archive writer also
restores frequency-dependent dispersion phase offsets: channel profiles must
be re-aligned before summing, even though PulsarX folded the data correctly.
Never synthesize the injected pulse model to make these plots.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np


def restore_pulsarx_dispersion(intensity, frequencies, periods, dm):
    """Undo channel phase offsets added by PulsarX when writing its .ar archive.

    Parameters:
      intensity: calibrated cube [subintegration, channel, pulse-phase bin]
      frequencies: SUBINT DAT_FREQ [subintegration, channel], in MHz
      periods: SUBINT PERIOD [subintegration], in seconds
      dm: archive's SUBINT DM in pc/cm^3

    This is channel *phase alignment for PulsarX folded archives*, not
    raw-sample dedispersion. Never apply it twice or to an archive that
    has already been dedispersed by another utility.
    """
    cube = np.asarray(intensity, dtype=np.float64)
    freqs = np.asarray(frequencies, dtype=np.float64)
    periods = np.asarray(periods, dtype=np.float64)
    if cube.ndim != 3:
        raise ValueError("Expected folded cube [subint, channel, phase]")
    nt, nc, nb = cube.shape
    if freqs.shape != (nt, nc) or periods.shape != (nt,):
        raise ValueError("Frequency or period array does not match folded cube")
    if not np.isfinite(dm) or not np.all(np.isfinite(freqs)) or np.any(freqs <= 0):
        raise ValueError("Invalid archive DM or DAT_FREQ values")
    if not np.all(np.isfinite(periods)) or np.any(periods <= 0):
        raise ValueError("Invalid SUBINT folding periods")

    # PulsarX's ArchiveWriter::get_dedispersion() applies a positive
    # DM-dependent channel roll when exporting the PSRFITS .ar. Undo it.
    # Match the writer's per-subintegration frequency convention via PERIOD.
    k_dm = 4148.808  # s MHz^2 per pc cm^-3
    fmax = np.max(freqs, axis=1, keepdims=True)
    delays = k_dm * float(dm) * (freqs**-2 - fmax**-2)
    shifts = np.rint(delays * nb / periods[:, None]).astype(np.int64)
    indices = (np.arange(nb)[None, None, :] + shifts[:, :, None]) % nb
    return np.take_along_axis(cube, indices, axis=2)


def read_dedispersed_pulsarx_archive(archive):
    """Load and align a PulsarX PSRFITS .ar, honoring scales, offsets, weights.

    The returned dictionary holds aligned channel profiles, not a fitted
    detection statistic. Handles the original PulsarX archive writer's
    dispersive phase offsets, not arbitrary PSRFITS archives.
    """
    from astropy.io import fits

    with fits.open(archive, memmap=False) as hdul:
        sub = hdul["SUBINT"]
        data = np.asarray(sub.data["DATA"], dtype=np.float64)
        if data.ndim != 4 or data.shape[1] != 1:
            raise ValueError(
                f"Expected total-intensity archive [time, 1, chan, phase], got {data.shape}"
            )
        nt, _, nc, nb = data.shape
        if nt == 0 or nc == 0 or nb < 2:
            raise ValueError("PulsarX archive has empty or invalid folded data")
        scales = np.asarray(sub.data["DAT_SCL"], dtype=np.float64).reshape(nt, nc)
        offsets = np.asarray(sub.data["DAT_OFFS"], dtype=np.float64).reshape(nt, nc)
        weights = np.asarray(sub.data["DAT_WTS"], dtype=np.float64).reshape(nt, nc)
        frequencies = np.asarray(sub.data["DAT_FREQ"], dtype=np.float64).reshape(nt, nc)
        periods = np.asarray(sub.data["PERIOD"], dtype=np.float64).reshape(nt)
        tsub = np.asarray(sub.data["TSUBINT"], dtype=np.float64).reshape(nt)
        tmid = np.asarray(sub.data["OFFS_SUB"], dtype=np.float64).reshape(nt)
        if "DM" not in sub.header:
            raise ValueError("SUBINT DM is missing; cannot align archive channels safely")
        dm = float(sub.header["DM"])
        calibrated = data[:, 0, :, :] * scales[:, :, None] + offsets[:, :, None]

    aligned = restore_pulsarx_dispersion(calibrated, frequencies, periods, dm)
    weights = np.where(np.isfinite(weights) & (weights > 0), weights, 0.0)
    if not np.all(np.isfinite(aligned)):
        raise ValueError("PulsarX archive contains non-finite calibrated intensities")
    if weights.sum() <= 0:
        raise ValueError("PulsarX archive has no valid channel weights")
    return dict(intensity=aligned, weights=weights, frequencies=frequencies,
                periods=periods, tsub=tsub, tmid=tmid, dm=dm)


def plot_folded_archive(archive, output, *, label="PulsarX folded pulse"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    folded = read_dedispersed_pulsarx_archive(archive)
    cube = folded["intensity"]
    weights = folded["weights"]
    nb = cube.shape[-1]
    baseline = np.median(cube, axis=-1, keepdims=True)
    profile = np.einsum("tcb,tc->b", cube - baseline, weights) / weights.sum()
    if not np.all(np.isfinite(profile)) or np.ptp(profile) <= 0:
        raise ValueError("PulsarX archive has no finite phase variation")

    phase = np.arange(nb) / nb
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.step(phase, profile, where="mid", linewidth=1.8)
    ax.set(xlabel="Pulse phase (cycles)",
           ylabel="Baseline-subtracted archive intensity",
           title=label, xlim=(0, 1))
    ax.grid(alpha=0.2)
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160)
    plt.close(fig)
    return output
