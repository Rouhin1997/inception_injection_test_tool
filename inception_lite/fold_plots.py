"""Portable headless pulse-profile diagnostic from a real PulsarX .ar archive.

PulsarX native PlotX PNG output is not guaranteed in headless containers.
This module produces a minimal phase-profile PNG from the folded PSRFITS
archive, *not* by regenerating the injected model.
"""
from __future__ import annotations

from pathlib import Path


def plot_folded_archive(archive, output, *, label="PulsarX folded pulse"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.io import fits

    with fits.open(archive, memmap=False) as hdul:
        ext = hdul["SUBINT"]
        nbin = int(ext.header["NBIN"])
        if nbin < 2:
            raise ValueError("PulsarX archive has invalid NBIN")
        raw = np.asarray(ext.data["DATA"], dtype=np.float64)
        if raw.size == 0 or raw.size % nbin:
            raise ValueError("PulsarX archive has empty or malformed DATA")
        if raw.shape[-1] != nbin:
            raw = raw.reshape(-1, nbin)
        profile = np.mean(raw.reshape(-1, nbin), axis=0)
        # PulsarX archives may store quantized DATA with DAT_SCL/OFFS.
        # For an eye-check, relative phase/profile shape is sufficient.
        if not np.all(np.isfinite(profile)):
            raise ValueError("PulsarX archive has non-finite folded profile")
        if np.ptp(profile) == 0:
            raise ValueError("PulsarX archive has a flat folded profile")
    phase = np.arange(nbin) / nbin
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.step(phase, profile, where="mid", linewidth=1.8)
    ax.set(xlabel="Pulse phase (cycles)", ylabel="Folded archive intensity",
           title=label, xlim=(0, 1))
    ax.grid(alpha=0.2)
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160)
    plt.close(fig)
    return output
