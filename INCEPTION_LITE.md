# INCEPTION-Lite: Nextflow-free pulsar injection

The INCEPTION-Lite runner uses the **original upstream INCEPTION Python
injector**, unchanged, but does not require Nextflow, MPIfR-specific paths,
MPI or a GPU. Develop only on dev/pulsejet-injections; keep main clean.

## Basic use

Install Python 3.11 and scientific dependencies:

    python -m pip install -r requirements-inception-lite.txt

Make a small noise-only, 8-bit SIGPROC filterbank template:

    python -m inception_lite.make_filterbank --output /scratch/noise.fil \
      --samples 65536 --nchans 8 --tsamp 0.001 --seed 771

Copy configs/inception_lite_single.json and replace the filterbank and
output_dir placeholder paths with your actual paths. Then:

    python -m inception_lite --config my_single.json --dry-run
    python -m inception_lite --config my_single.json

To choose automated **PulsarX truth folding** on each run without editing
the JSON configuration, use:

    python -m inception_lite --config my_single.json --perform-folds yes
    python -m inception_lite --config my_single.json --perform-folds no

The spelling `--perform_folds` is also accepted. The command-line choice
takes precedence over `fold.enabled` in the JSON. Omitting it uses the
existing JSON value (default: no folding). `yes` means inject, then
automatically call PulsarX to fold each known injected source and save
its profile plot, archive and log under `folds/<pulsar-id>/`.
`no` means inject without running PulsarX. `--dry-run` always only
lists planned cases; it never generates data or folds anything.

**Important:** Folding requires a working `psrfold_fil` executable and
a valid `fold.template` in the chosen environment. Our small CPU-only
INCEPTION-Lite Docker image currently tests folding orchestration with
a fake executable, but does *not* include PulsarX itself; the expanded
PULSEJET development container does include PulsarX. We must complete
that image integration and a real PulsarX smoke test before `yes`
can be used inside the injector-only SIF.

The input and injected output remain SIGPROC `.fil` format.
Search-mode PSRFITS `.fits` input is not directly supported; it
needs a validated conversion to SIGPROC or a separately implemented
PSRFITS read/write path.

To create a parameter grid, use configs/inception_lite_grid.json, or
configs/inception_lite_random.json for seeded random parameter draws.
For a single case in a SLURM array:

    python -m inception_lite --config my_grid.json --case case_0004

The generic SLURM launcher also accepts an optional third argument
to override folding for every array task:

    sbatch --array=1-9 scripts/slurm_inception_lite.sh my_grid.json sif yes

Use `no` instead of `yes` to skip folding; omit the third argument
to follow the JSON setting. The SIF still needs PulsarX when folding
is enabled.

Modes: together (many pulsars in ONE filterbank), separate (one per
filterbank), grid (Cartesian product of parameters), random (reproducible
parameter draws). The runner refuses to overwrite preexisting case
directories. You must supply an explicit global seed.

## Output directories

Each case_0001/ folder contains the injected SIGPROC .fil, the resolved
parameter JSON, the original INCEPTION report and folding input files,
inject.log, and run_manifest.json with input/output names and random
seed. With optional PulsarX folding it also has:

    case_0001/folds/PSR_ID/
        PSR_ID_phase_truth.cands
        fold.log
        [PulsarX-generated PNG figures and .ar archives]

The standalone upstream injector itself DOES NOT automatically fold.
Its Nextflow orchestration separately invokes PulsarX, PRESTO and DSPSR
and collects their plotting products for visual inspection.

## Optional automatic PulsarX folding

Set the fold property in your JSON to:

    "fold": {
      "enabled": true,
      "binary": "psrfold_fil",
      "template": "/usr/local/share/pulsarx/meerkat_fold.template",
      "nbin": 64
    }

Requires a real psrfold_fil executable and compatible template. Initial
supported case: simple topocentric acceleration and jerk with zero
relative velocity at the midpoint, no intrinsic frequency derivatives,
PINT ephemeris or full binary orbital model. The runner refuses unsupported
cases rather than silently generating incorrect fold inputs.

The wrapper produces a separate seven-column candidate file consistent
with the actual INCEPTION first-order Doppler phase. It does not change
the original .candfile, whose tiny acceleration-squared F2 contribution
differs from the injected phase model. F1 is explicitly included, so
acceleration is not applied twice in PulsarX. The fold epoch is the
midpoint of the full injected filterbank, not an arbitrary PULSEJET
FFT segment.

## Physical conventions

With topocentric injection, no extra intrinsic spin derivatives,
constant acceleration a and jerk j, define tau from the acceleration
reference epoch. The original injector implements:

    velocity(tau) = a*tau + j*tau*tau/2
    phase(tau) = f0*[tau - a*tau*tau/(2c) - j*tau**3/(6c)]
    F1(reference) = -f0*a/c
    F2(reference) = -f0*j/c

Use frame="topo", ACCEPOCH=0.5 and do NOT explicitly set PEPOCH for
our initial tests. In the upstream source, spin_ref and ACCEPOCH
both shift the time argument; assigning PEPOCH=0.5 simultaneously
could cause a double epoch offset. The provided configs avoid this.

PULSEJET's Phase A time transformation comes from an exact inversion
of a truncated arrival-time polynomial. INCEPTION instead uses the
first-order Doppler model. Their difference starts at order (v/c)^2;
we validate signs and epochs and quantify the remainder in tests.

## Validation performed in GitHub Actions

- Parameter planning, deterministic random seeds, file safety
- Valid SIGPROC filterbank generation and output header checks
- Independent symbolic Doppler phase, first/second frequency derivatives
- Real end-to-end injections for stationary and moving pulsars
- Re-running with equal seed gives byte-identical filterbanks
- Coherent profile recovery from a stationary injected source
- Deliberately strong jerk waveform demodulation compared with wrong sign
- PulsarX fold command handoff tested with a fake executable

A fake folding executable does NOT prove real PulsarX folding: actual
plots and archives must be verified after building/running the image.
Likewise, PULSEJET GPU candidate recovery remains a separate step.

## Deployment

Dockerfile.inception-lite defines a CPU-only Python environment; it
can be converted into a read-only SIF without creating a sandbox.
scripts/slurm_inception_lite.sh launches individual cases on any SLURM
cluster with Python or SingularityCE, without Nextflow or fixed MPIfR
partitions. Pass source and data paths through writable mounts.


## Published combined image: real injection + PulsarX folding (2026-10-20)

An integrated INCEPTION-Lite / PulsarX Docker image has been **built, tested
and published** under the existing Docker Hub repository
`rouhin1997/pulsejet`, with **separate injector-specific tags**:

    rouhin1997/pulsejet:inception-lite-5a9cd1c9b252ffd1797575267c4528fa5b8e81c5
    rouhin1997/pulsejet:inception-lite-latest

Published from INCEPTION-Lite commit
`5a9cd1c9b252ffd1797575267c4528fa5b8e81c5`,
GitHub Actions: https://github.com/Rouhin1997/PULSEJET_beta/actions/runs/39423881379

The combined image layers Python/SymPy and this repository's code over the
previously compiled `PULSEJET` CUDA 12.3.2 and PulsarX image, so it is larger
than the separate minimal CPU-only image. You only need the one SIF if you
want both injection and PulsarX truth folding. **The original master/main
branches and production Docker Hub `pulsejet:latest` tag are unchanged.**

On the HPC's SingularityCE 3.11:

    singularity pull inception_lite.sif \
      docker://rouhin1997/pulsejet:inception-lite-5a9cd1c9b252ffd1797575267c4528fa5b8e81c5

Then, assuming your JSON file and SIGPROC .fil are on writable host scratch
and the same absolute paths are visible within Singularity:

    singularity exec --bind /scratch:/scratch inception_lite.sif \
      python3 -m inception_lite --config /scratch/injections.json \
      --dry-run

    singularity exec --bind /scratch:/scratch inception_lite.sif \
      python3 -m inception_lite --config /scratch/injections.json \
      --perform-folds yes

For injection and CPU PulsarX folding, `--nv` is not needed. It is needed
when you subsequently execute the CUDA PULSEJET search on the cluster.

Real GitHub Docker image integration testing injected a source with
`P0=0.02 s, A0=400 m/s^2, A1=-3 m/s^3` and verified that real
`psrfold_fil` created a non-empty folded `.ar` archive and PNG.
PulsarX's native PlotX pipeline produced no PNG in that headless CI
container; INCEPTION-Lite now renders a small PNG **from the actual
PulsarX folded PSRFITS archive**. It does not invent the profile from
injection parameters. Files are under
`case_0001/folds/<pulsar-id>/`.

**Boundaries:** this is a successful software/integration test, not yet
a scientific recovery or sensitivity benchmark on the actual A40 cluster.
The plots should be examined, and PULSEJET GPU recovery still needs
independent injection-to-search validation.


## PulsarX archive plot correction (2026-10-27)

The PulsarX `psrfold_fil` archive writer reapplies channel-dependent
dispersion delays when exporting the frequency-resolved PSRFITS `.ar`.
Plotting scripts **must undo those phase shifts before averaging channels**.
Both the headless PNG fallback (`inception_lite/fold_plots.py`) and the
three-panel scientific diagnostic (`scripts/plot_pulsarx_diagnostic.py`)
now share a correction using `SUBINT.DM`, per-subintegration `PERIOD`,
`DAT_FREQ`, `DAT_SCL`, `DAT_OFFS`, and `DAT_WTS`.
Do not apply it twice or to an archive already dedispersed by other tools.

On existing M62 W3 injections, corrected archive peak/MAD improved from
2.58 to 31.40 (I11, stationary), 2.78 to 12.80 (I12, acceleration),
2.23 to 33.98 (I13, jerk), and 2.30 to 23.66 (I14, mixed).
These are **profile-contrast ratios, not calibrated detection S/N**.
PRESTO independently recovered stationary I11 at about 60.8 sigma.
This validates those *truth folds*, not PULSEJET candidate recovery.

Already generated PNGs do not change automatically; rerun the diagnostic
plotter on existing `.ar` archives, or regenerate the fallback images.
Container images are immutable and require rebuilding to incorporate this fix.


### M62 20-source truth-fold visual quality-control batch

To regenerate the *existing* five PulsarX truth-fold archives in each of
W1/W2/W3/W4 (20 total) as corrected three-panel diagnostics and produce
a machine-readable QC table, run from a **host source checkout** on the
cluster (not from the older immutable `inception_lite.sif` source):

```bash
REPO=/home/rnag/PULSEJET_validation/inception_injection_test_tool
ROOT=/home/rnag/PULSEJET_validation/stage6_m62/m62_15m_injections
SIF=/home/rnag/PULSEJET_validation/images/inception_lite.sif

singularity exec --cleanenv "$SIF" env PYTHONNOUSERSITE=1 python3 \
  "$REPO/scripts/batch_m62_pulsarx_diagnostics.py" \
  --root "$ROOT" \
  --out "$ROOT/pulsarx_diagnostics_corrected"
```

Output: 20 corrected `.png` diagnostics and
`pulsarx_diagnostics_corrected/m62_20_truth_fold_qc.csv`. The CSV records
archive DM versus injected DM and raw/corrected profile peak/MAD as a
**descriptive, non-calibrated** contrast. The command checks that it found
four five-source manifests and exactly one `.ar` per source, and returns a
nonzero exit status on missing/failed plot outputs. It does not regenerate
`.fil` or `.ar` files; the previous plots in `folds/` are unchanged.
Examine all PNGs rather than declaring a scientific pass from the peak/MAD
metric alone. Low injected S/N or DM=0 may show little or no contrast gain.

PRESTO `prepfold` writes a different `.pfd` format and is not read by
this PulsarX-specific PSRFITS handler. A future PRESTO plotter must parse
its own native products; this fix cannot affect PRESTO folding.


### Completed M62 twenty-source visual validation checkpoint (2026-10-27)

The user ran the corrected M62 batch on all four completed observational
injection windows and shared `m62_20_truth_fold_qc.csv` plus the full
twenty three-panel PNG set. **All 20 PulsarX known-parameter truth archives
were processed** (5 each W1–W4); archive DM matches injected DM in each,
and every archive has 89×32×128 subintegrations/channels/phase bins with
positive channel weights. After channel alignment, median integrated
profile peak/MAD is 25.128 (range 7.103–58.391). Every source with
nonzero additional DM improved; the DM=0 W1/I01 control stayed at 52.484.
Even the faint injected target-S/N cases W2/I10=10, W4/I18=12 and
W4/I20=15 show visually identifiable folded profile peaks, though weaker
time/frequency waterfalls. Visual inspection found no obvious systematic
pulse-phase drift across the windows. Peak/MAD is **not** calibrated
detection S/N or completeness. A separate PRESTO fixed-parameter control
for W3/I11 yielded `~60.8 sigma` in its `.bestprof`.

This passes a **truth-fold plotting/consistency gate**, not a blind PULSEJET
candidate-recovery test. The saved PulsarX `.ar`, original injected
`.fil`, injection physics and existing images were not modified by this
correction. New plotted files:
`$HOME/PULSEJET_validation/stage6_m62/m62_15m_injections/pulsarx_diagnostics_corrected/`.
The companion PULSEJET development log and cluster handoff describe the
next controlled GPU searches and a known still-unfixed `PeakFinder`
frequency-bound bug in the older pinned PULSEJET image.
