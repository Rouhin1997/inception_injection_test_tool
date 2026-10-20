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
