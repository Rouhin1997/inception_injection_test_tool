# INCEPTION-Lite / PULSEJET — 2026-10-20 cross-project restart checkpoint

**The user is pausing until tomorrow and will resume in a fresh chat, beginning with Stage 1: Singularity container verification on their A40 cluster.**

**Authoritative detailed roadmap and PULSEJET scientific backlog:**
https://github.com/Rouhin1997/PULSEJET_beta/blob/dev/polynomial-search/tests/PHASE_B_CLUSTER_HANDOFF.md

This short note preserves this repository's role in the plan.

## Repository policy

- Work on \`dev/pulsejet-injections\` exclusively, preserve \`main\` upstream.
- Original INCEPTION injection engine \`injector/\` remains unchanged.
- Lightweight Python runner in \`inception_lite/\`, no Nextflow or MPIfR cluster configuration needed.

## Published validated image

Use the **immutable** version:

    docker://rouhin1997/pulsejet:inception-lite-5a9cd1c9b252ffd1797575267c4528fa5b8e81c5

Image source commit: \`5a9cd1c9b252ffd1797575267c4528fa5b8e81c5\`
Build+actual PulsarX fold test: https://github.com/Rouhin1997/PULSEJET_beta/actions/runs/39423881379

Original PULSEJET-only image, for baseline independent tests:

    docker://rouhin1997/pulsejet:dev-polynomial-513e80d1ace21868fd3d3d37734bcd79782abce3

Both tags live on \`rouhin1997/pulsejet\` Docker Hub repository. The INCEPTION-Lite image *inherits* PULSEJET/PulsarX and adds injector/runner, so it can perform injection + truth folding in one SIF. No cluster SIF runtime tests performed yet.

## First step tomorrow

    singularity pull inception_lite.sif \
      docker://rouhin1997/pulsejet:inception-lite-5a9cd1c9b252ffd1797575267c4528fa5b8e81c5

    singularity exec inception_lite.sif python3 -m inception_lite --help
    singularity exec inception_lite.sif psrfold_fil --help

Check host/SIF bind mounts and storage before running actual data. For PULSEJET GPU search use \`singularity exec --nv\` on an allocated A40 compute node.

## Capabilities already verified in CI

- Mode \`together\`: multiple pulsars in one input SIGPROC .fil and one output .fil.
- Modes \`separate\`, \`grid\`, \`random\`: distinct output filterbanks, reproducible seeds/manifests.
- CLI: \`python3 -m inception_lite --config run.json --dry-run\` preview only; \`--perform-folds yes|no\` overrides JSON folding setting.
- Only .fil injection currently supported; no direct .fits/PSRFITS injection.
- Physical injection uses first-order Doppler acceleration A0 (m/s²), jerk A1 (m/s³); use \`frame: topo\`, \`ACCEPOCH: 0.5\`, no explicit PEPOCH for first controlled runs. Same physical epoch and file interval must be matched against PULSEJET polynomial search.
- CI: single and multiple real filterbank injections, re-seeding reproducibility, known-frequency coherence, extreme jerk sign check all passed.
- Integrated image: **real** \`psrfold_fil\` produced a folded \`.ar\` for a test jerk pulsar; because native PlotX PNG was absent in headless CI, new \`inception_lite/fold_plots.py\` rendered a PNG from actual folded archive; CI verified both. Do not mistake this for PULSEJET search recovery.

## Scientific workflow to resume after cluster checks

1. Stage 1: Pull and verify both SIFs, A40 runtime access, library/CPU compatibility and scratch paths.
2. Stage 2: Independent CPU vs actual GPU resampler parity; create the GPU harness, not built yet.
3. Stage 3: Simple single-injection truth-fold, PULSEJET search, XML-to-PulsarX conversion and recovered-parameter fold; compare a=j=0, both signs, combined a/j.
4. Stage 4: Pilot parameter sweep and matched acceleration-only baseline; full sensitivity/mismatch survey later.
5. Stage 5: Multi-source source association with \`--candidate_diagnostics\`, only revise distillation on evidence.
6. Stage 6: Large-scale scientifically robust sensitivity, performance and survey validation.

**Time realism:** 2–3 days at 5–6 hours/day (10–18 h) may suffice for stages 1–3 and perhaps a tiny Stage 4 pilot if cluster runs cooperate; cannot guarantee complete stages 4–6 within that window.

Start next chat with: "Resume Stage 1 from the PULSEJET 2026-10-20 cluster handoff."
