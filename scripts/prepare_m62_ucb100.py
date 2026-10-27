#!/usr/bin/env python3
"""Deterministic, 20x5 full-orbit ultra-compact M62 injection campaign.

Creates configuration JSONs and physical-truth CSV only. Does NOT inject.
Design: 90 circular binaries (25-75 min), 10 isolated controls. The sample
is an intentionally jerk-stratified stress test, not a population prior.
Uses existing 4 x ~15 minute observational windows, reused 5 times each.
For binary injections, do NOT also set polynomial A0/A1.
"""
import argparse
import csv
import hashlib
import json
import math
import random
from pathlib import Path
from inception_lite.runner import plan_cases, read_filterbank_header

G = 6.67430e-11
MSUN = 1.98847e30
C = 299792458.0
SEED = 20261027
BINS = [
    ("isolated", 10, None, None),
    ("j_000_025", 15, 0.0, 0.25),
    ("j_025_060", 20, 0.25, 0.60),
    ("j_060_120", 25, 0.60, 1.20),
    ("j_120_210", 30, 1.20, 2.10),
]
FIELDS = [
    "id", "batch", "window", "motion_class", "f_spin_hz", "p_spin_s",
    "dm_add", "snr_nominal", "duty_cycle", "pb_min", "mp_msun",
    "mc_msun", "inc_deg", "phase_mid_rad", "x_lt_s", "t0_mjd",
    "acc_mid_mps2", "jerk_mid_mps3", "amax_mps2", "jmax_mps3",
    "source_seed",
]


def orbital_params(pb_min, mc, inc_deg, phase_mid):
    pb_sec = pb_min * 60.0
    omega = 2.0 * math.pi / pb_sec
    mtotal = 1.4 + mc
    semimajor = (G * mtotal * MSUN / omega**2) ** (1.0 / 3.0)
    x = (semimajor / C) * mc / mtotal * math.sin(math.radians(inc_deg))
    amax = C * x * omega**2
    jmax = amax * omega
    # T0 selected below so phase_mid = Omega * (t_mid - T0).
    return x, -amax * math.sin(phase_mid), -jmax * math.cos(phase_mid), amax, jmax


def build_population(seed=SEED):
    rng = random.Random(seed)
    by_class = []
    for kind, count, jlo, jhi in BINS:
        accepted = []
        for number in range(count):
            if kind == "isolated":
                accepted.append(dict(motion_class=kind, pb_min=None,
                    mp_msun=None, mc_msun=None, inc_deg=None,
                    phase_mid_rad=None, x_lt_s=None,
                    acc_mid_mps2=0.0, jerk_mid_mps3=0.0,
                    amax_mps2=0.0, jmax_mps3=0.0))
                continue
            for attempt in range(100000):
                pb_min = rng.uniform(25.0, 75.0)
                mc = rng.uniform(0.03, 0.20)
                inc = math.degrees(math.acos(rng.uniform(0, 1)))
                phase = rng.uniform(0, 2 * math.pi)
                x, acc, jerk, amax, jmax = orbital_params(pb_min, mc, inc, phase)
                if (jlo <= abs(jerk) < jhi and amax <= 550.0
                        and jmax <= 2.5):
                    accepted.append(dict(motion_class=kind, pb_min=pb_min,
                        mp_msun=1.4, mc_msun=mc, inc_deg=inc,
                        phase_mid_rad=phase, x_lt_s=x,
                        acc_mid_mps2=acc, jerk_mid_mps3=jerk,
                        amax_mps2=amax, jmax_mps3=jmax))
                    break
            else:
                raise RuntimeError("Could not populate constrained class " + kind)
        by_class.append(accepted)
    # Round-robin classes across batches, shuffling deterministically within
    # each stratum to prevent all extreme systems accumulating on one node.
    ordered = []
    while any(by_class):
        for group in by_class:
            if group:
                ordered.append(group.pop(0))
    assert len(ordered) == 100
    return ordered


def generate(root, destination, seed=SEED):
    if destination.exists():
        raise FileExistsError("Refusing to overwrite existing campaign: " + str(destination))
    backgrounds = []
    for i in range(1, 5):
        path = root / "stage6_m62" / "m62_15m_windows" / ("W%d_M62_background.fil" % i)
        if not path.is_file():
            raise FileNotFoundError("Missing validated M62 background: " + str(path))
        hdr = read_filterbank_header(path)
        if hdr["nbits"] != 8 or hdr["actual_nsamples"] < 4194304:
            raise RuntimeError("Unexpected M62 input format/length: " + str(path))
        backgrounds.append((path, hdr))
    rng = random.Random(seed ^ 0xA5A51234)
    pop = build_population(seed)
    rows = []
    configs = []
    for index in range(100):
        batch = index // 5 + 1
        slot = index % 5
        window = (batch - 1) % 4 + 1
        background, hdr = backgrounds[window - 1]
        item = dict(pop[index])
        # Well-separated frequency strata avoid harmonics of other injected
        # sources falling directly on their fundamental spin frequencies.
        lo, hi = [(231.,255.),(267.,292.),(305.,330.),(342.,367.),(380.,403.)][slot]
        f = rng.uniform(lo, hi)
        dm = rng.uniform(7.2, 8.8)
        snr = rng.uniform(15.0, 85.0)
        duty = rng.uniform(.03, .08)
        name = "UCB_B%02d_S%02d" % (batch, slot + 1)
        record = dict(id=name, batch=batch, window=window,
                      f_spin_hz=f, p_spin_s=1/f, dm_add=dm,
                      snr_nominal=snr, duty_cycle=duty,
                      source_seed=seed + index + 1, **item)
        if item["pb_min"] is not None:
            duration = hdr["actual_nsamples"] * hdr["tsamp"]
            omega = 2*math.pi/(item["pb_min"]*60.0)
            # INCEPTION orbital phase uses t relative to T0; store T0 at
            # the physical source's chosen orbital phase at the FULL
            # input-filterbank midpoint.
            t0 = hdr["tstart"] + (duration/2 - item["phase_mid_rad"]/omega)/86400.0
            record["t0_mjd"] = t0
        else:
            record["t0_mjd"] = None
        rows.append(record)
    destination.mkdir(parents=True)  # only after all source files validated
    try:
        for batch in range(1, 21):
            window = (batch-1) % 4 + 1
            background, _ = backgrounds[window-1]
            pulses = []
            for item in rows[(batch-1)*5:batch*5]:
                p = dict(ID=item["id"],P0=item["p_spin_s"],DM=item["dm_add"],
                         SNR=item["snr_nominal"],duty_cycle=item["duty_cycle"],
                         seed=item["source_seed"])
                if item["pb_min"] is not None:
                    p.update(binary_period=item["pb_min"]/60.0,  # HOURS
                        M1=1.4,M2=item["mc_msun"],inc=item["inc_deg"],
                        T0=item["t0_mjd"],ecc=0.0,AoP=0.0)
                pulses.append(p)
            cfg = dict(id="ucb100_batch_%02d" % batch,
                seed=seed+1000*batch,mode="together",
                filterbank=str(background),
                output_dir=str(destination/"injections"/("batch_%02d" % batch)),
                offline=True,ncpu=48,stats_samples=200000,gulp_size_GB=0.02,
                psr_global=dict(frame="topo",ACCEPOCH=0.5,DM_ref="inf",
                                create_parfile="par"),
                pulsars=pulses,
                # Full orbital models are not supported by the current
                # automatic PulsarX polynomial truth-fold wrapper.
                fold=dict(enabled=False))
            cases=plan_cases(cfg)
            assert len(cases)==1 and len(cases[0][1])==5
            conf=destination/"configs"/("batch_%02d.json" % batch)
            conf.parent.mkdir(parents=True,exist_ok=True)
            with conf.open("x") as fp:
                json.dump(cfg,fp,indent=2)
                fp.write("\n")
            configs.append(conf)
        with (destination/"ucb100_truth.csv").open("x",newline="") as fp:
            writer=csv.DictWriter(fp,fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        manifest = {
            "experiment": "M62 UCB100 20x5 full orbital injections",
            "seed": seed, "total_sources":len(rows),"batches":20,
            "binary_sources":sum(x["pb_min"] is not None for x in rows),
            "isolated_controls":sum(x["pb_min"] is None for x in rows),
            "source_files": [{
                "path":str(p),"size_bytes":p.stat().st_size,
                "nsamples":h["actual_nsamples"],"tsamp_s":h["tsamp"],
                "tstart_mjd":h["tstart"],"nbits":h["nbits"],
                "nchans":h["nchans"]} for p,h in backgrounds],
            "config_sha256": {p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in configs},
            "search_envelope": {"acc_mps2":[-550,550],
                                "jerk_mps3":[-2.5,2.5],
                                "dm_add":[7.2,8.8]},
            "notes":["Full Keplerian circular-orbit injection for 90 sources, 10 stationary controls.",
                     "Physical midpoint a/j from same orbit; no polynomial A0/A1 is injected with binary.",
                     "25-75 min binary period; 0.03-0.20 Msun companion; population is jerk-stratified, not an unbiased astrophysical prior.",
                     "Automatic truth folding disabled: current wrapper cannot fold full binary orbits.",
                     "Injected DM is ADDITIONAL synthetic DM, not original M62 total DM.",
                     "Identical backgrounds reused five times; samples are not statistically independent noise observations."]}
        with (destination/"campaign_manifest.json").open("x") as fp:
            json.dump(manifest,fp,indent=2)
            fp.write("\n")
        print("Prepared",len(configs),"batch configs, 100 pulsars, 90 full binary orbits.")
        print("Plan:",destination)
        print("Source truth:",destination/"ucb100_truth.csv")
        print("Manifest:",destination/"campaign_manifest.json")
        print("Ranges: |a_mid|<=%.3f m/s^2, |j_mid|<=%.6f m/s^3, max j amplitude=%.6f" % (
            max(abs(x["acc_mid_mps2"]) for x in rows),
            max(abs(x["jerk_mid_mps3"]) for x in rows),
            max(x["jmax_mps3"] for x in rows)))
        print("Estimated total injected .fil bytes (20 outputs):",sum(
            backgrounds[(b-1)%4][0].stat().st_size for b in range(1,21)))
        return rows
    except Exception:
        print("WARNING: partial plan may exist in",destination,
              "; inspect manually, do NOT overwrite.",flush=True)
        raise


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root",type=Path,default=Path("/home/rnag/PULSEJET_validation"))
    ap.add_argument("--out",type=Path)
    args=ap.parse_args()
    root=args.root.resolve()
    dest=(args.out or root/"stage10_ucb100").resolve()
    generate(root,dest)


if __name__=="__main__":
    main()
