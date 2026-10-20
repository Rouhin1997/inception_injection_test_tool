"""Standalone, Nextflow-free INCEPTION injection controller.

The original injector is executed unchanged as a subprocess. This module
controls reproducible configurations, files, and optional PulsarX folding.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import itertools
import importlib.metadata
import json
import math
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
INJECTOR = REPO / "injector" / "SCRIPT_inject_pulsars.py"
SAFE_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_-]{0,63}$")


def positive_number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be numeric")
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def check_pulsars(pulsars):
    if not isinstance(pulsars, list) or not pulsars:
        raise ValueError("pulsars must be a nonempty list")
    ids = set()
    for pulsar in pulsars:
        if not isinstance(pulsar, dict):
            raise ValueError("Each pulsar must be an object")
        identifier = pulsar.get("ID", "")
        if not isinstance(identifier, str) or not SAFE_ID.fullmatch(identifier):
            raise ValueError(f"Unsafe or missing pulsar ID: {identifier!r}")
        if identifier in ids:
            raise ValueError(f"Duplicate pulsar ID: {identifier}")
        ids.add(identifier)
        if ("P0" not in pulsar) and ("F0" not in pulsar):
            raise ValueError(f"{identifier}: provide P0 or F0")
        positive_number(pulsar.get("P0", pulsar.get("F0")), identifier + " period/frequency")
        if "SNR" not in pulsar:
            raise ValueError(f"{identifier}: SNR required")
        positive_number(pulsar["SNR"], identifier + " SNR")
        for param in ("A0", "A1", "DM"):
            if param in pulsar:
                v = pulsar[param]
                if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
                    raise ValueError(f"{identifier}: {param} must be finite numeric")
        if pulsar.get("DM", 0) < 0:
            raise ValueError(f"{identifier}: negative dispersion measure")


def draw(spec, rng):
    if not isinstance(spec, dict):
        raise ValueError("Random parameter specifications must be objects")
    dist = spec.get("distribution")
    if dist == "choice":
        values = spec.get("values")
        if not isinstance(values, list) or not values:
            raise ValueError("choice requires nonempty values")
        return rng.choice(values)
    low, high = spec.get("min"), spec.get("max")
    if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x)
           for x in (low, high)) or low > high:
        raise ValueError("random min/max must be finite, numeric, and ordered")
    if dist == "uniform":
        return rng.uniform(low, high)
    if dist == "loguniform":
        if low <= 0:
            raise ValueError("loguniform bounds must be positive")
        return math.exp(rng.uniform(math.log(low), math.log(high)))
    raise ValueError(f"Unsupported random distribution: {dist}")


def plan_cases(config):
    """Return list of (case_id, individual pulsar dictionaries).

    modes: together, separate, grid, random.
    'together': inject all listed pulsars into the SAME filterbank.
    'separate': inject each listed pulsar into a separate filterbank.
    'grid' / 'random': one baseline pulsar per generated filterbank.
    """
    mode = config.get("mode", "together")
    pulsars = config.get("pulsars")
    check_pulsars(pulsars)
    seed = config.get("seed")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("An explicit nonnegative integer seed is required")
    base_id = config.get("id", "inception_lite")
    if not isinstance(base_id, str) or not SAFE_ID.fullmatch(base_id):
        raise ValueError("id must begin with a letter and contain only letters, digits, '_' or '-'")
    cases = []
    if mode == "together":
        cases = [("case_0001", copy.deepcopy(pulsars))]
    elif mode == "separate":
        cases = [(f"case_{i:04}", [copy.deepcopy(p)]) for i, p in enumerate(pulsars, 1)]
    elif mode in ("grid", "random"):
        if len(pulsars) != 1:
            raise ValueError(f"{mode} mode requires precisely one baseline pulsar")
        base = pulsars[0]
        if mode == "grid":
            grid = config.get("grid")
            if not isinstance(grid, dict) or not grid:
                raise ValueError("grid mode requires a nonempty grid object")
            if any(not isinstance(v, list) or not v for v in grid.values()):
                raise ValueError("each grid parameter must have a nonempty list")
            if "ID" in grid or "seed" in grid:
                raise ValueError("grid cannot override ID or seed")
            for n, combination in enumerate(itertools.product(*grid.values()), 1):
                p = {**copy.deepcopy(base), **dict(zip(grid, combination))}
                p["ID"] = f"{base['ID']}_{n:04}"
                cases.append((f"case_{n:04}", [p]))
        else:
            randomized = config.get("random")
            if not isinstance(randomized, dict):
                raise ValueError("random mode requires a random object")
            count = randomized.get("count")
            if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
                raise ValueError("random count must be a positive integer")
            specs = randomized.get("parameters")
            if not isinstance(specs, dict) or not specs or "ID" in specs or "seed" in specs:
                raise ValueError("random.parameters must be nonempty and cannot override ID/seed")
            rng = random.Random(seed)
            for n in range(1, count + 1):
                p = {**copy.deepcopy(base), **{k: draw(v, rng) for k, v in specs.items()}}
                p["ID"] = f"{base['ID']}_{n:04}"
                cases.append((f"case_{n:04}", [p]))
    else:
        raise ValueError(f"Unsupported mode: {mode}")
    if len(cases) > config.get("max_cases", 1000):
        raise ValueError("Too many planned cases; raise max_cases explicitly if intended")
    for _, batch in cases:
        check_pulsars(batch)
        for pulsar in batch:
            global_args = config.get("psr_global", {})
            if "frame" not in pulsar and "frame" not in global_args:
                raise ValueError("Set frame explicitly to topo or bary; upstream defaults are not appropriate for blind validation")
            a = pulsar.get("A0", global_args.get("A0", 0))
            j = pulsar.get("A1", global_args.get("A1", 0))
            if (a or j) and ("PEPOCH" in pulsar or "PEPOCH" in global_args):
                raise ValueError("Do not set PEPOCH with A0/A1 in the validated mode: original spin_ref and ACCEPOCH are both offsets")
    return cases


def case_seed(global_seed, name):
    # Stable across platforms and independent of Python's randomized hash().
    digest = hashlib.sha256(f"{global_seed}:{name}".encode()).digest()
    return int.from_bytes(digest[:4], "big")


def read_filterbank_header(path):
    """Read only the fields needed for safe fold epoch selection."""
    import struct
    fields = {}
    ints = {"machine_id", "telescope_id", "data_type", "nchans", "nbits", "nifs",
            "scan_number", "barycentric", "pulsarcentric", "nbeams", "ibeam", "nsamples"}
    floats = {"tstart", "tsamp", "fch1", "foff", "refdm", "az_start", "za_start",
              "src_raj", "src_dej"}
    strings = {"source_name", "rawdatafile"}
    def s(handle):
        raw = handle.read(4)
        if len(raw) != 4:
            raise ValueError("Unexpected EOF in filterbank header")
        n = struct.unpack("<i", raw)[0]
        if not 0 <= n <= 80:
            raise ValueError("Invalid SIGPROC header string length")
        return handle.read(n).decode("utf-8")
    with path.open("rb") as handle:
        if s(handle) != "HEADER_START":
            raise ValueError("Missing SIGPROC HEADER_START")
        while True:
            tag = s(handle)
            if tag == "HEADER_END":
                offset = handle.tell()
                break
            if tag in ints:
                fields[tag] = struct.unpack("<i", handle.read(4))[0]
            elif tag in floats:
                fields[tag] = struct.unpack("<d", handle.read(8))[0]
            elif tag in strings:
                fields[tag] = s(handle)
            elif tag == "signed":
                fields[tag] = struct.unpack("<b", handle.read(1))[0]
            else:
                raise ValueError(f"Unknown SIGPROC header field: {tag}")
    nbits, nchans = fields.get("nbits"), fields.get("nchans")
    if not isinstance(nbits, int) or not isinstance(nchans, int) or nbits <= 0 or nchans <= 0:
        raise ValueError("Invalid channel/bit count in filterbank")
    if not math.isfinite(fields.get("tsamp", math.nan)) or fields["tsamp"] <= 0:
        raise ValueError("Invalid sample interval")
    data_bits = (path.stat().st_size - offset) * 8
    bits_per_sample = nchans * nbits * fields.get("nifs", 1)
    if data_bits <= 0 or data_bits % bits_per_sample != 0:
        raise ValueError("Malformed filterbank data length")
    fields["actual_nsamples"] = data_bits // bits_per_sample
    fields["data_offset"] = offset
    return fields


def fold_truth(run_dir, output_fil, pulses, cfg, header):
    fold = cfg.get("fold", {})
    if not fold.get("enabled", False):
        return []
    if any(p.get("frame", cfg.get("psr_global", {}).get("frame", "bary")) != "topo" for p in pulses):
        raise ValueError("First-phase automatic folding supports topocentric injections only")
    if any(p.get("ACCEPOCH", cfg.get("psr_global", {}).get("ACCEPOCH", 0.5)) != 0.5 for p in pulses):
        raise ValueError("First-phase folding requires midpoint ACCEPOCH=0.5")
    template = Path(fold.get("template", ""))
    if not template.is_file():
        raise ValueError("fold.template must name an existing PulsarX template")
    binary = fold.get("binary", "psrfold_fil")
    if not shutil.which(binary):
        raise ValueError(f"Cannot find folding executable: {binary}")
    mid_epoch = header["tstart"] + header["actual_nsamples"] * header["tsamp"] / (2 * 86400)
    results = []
    for p in pulses:
        # Use a folding truth file consistent with INCEPTION's *injected*
        # first-order Doppler phase. Keep its own historical .candfile intact
        # for cross-checks: it includes a slightly different a^2/c^2 term.
        if p.get("binary_period") or p.get("pint_polycos") or any(
                k in p for k in ("F1", "F2", "P1", "P2", "presto_z", "presto_w")):
            raise ValueError("Automatic truth folding is currently limited to simple AX-only topocentric injections")
        if "PEPOCH" in p or "PEPOCH" in cfg.get("psr_global", {}):
            raise ValueError("Explicit PEPOCH needs a verified epoch transformation; omit it for midpoint folding")
        f0 = 1.0 / float(p["P0"]) if "P0" in p else float(p["F0"])
        a = float(p.get("A0", cfg.get("psr_global", {}).get("A0", 0.0)))
        jerk = float(p.get("A1", cfg.get("psr_global", {}).get("A1", 0.0)))
        dm = float(p.get("DM", cfg.get("psr_global", {}).get("DM", 0.0)))
        snr = float(p["SNR"])
        c = 299792458.0
        f1 = -a * f0 / c
        f2 = -jerk * f0 / c
        fold_dir = run_dir / "folds" / p["ID"]
        fold_dir.mkdir(parents=True, exist_ok=False)
        cand = fold_dir / f"{p['ID']}_phase_truth.cands"
        cand.write_text(
            "#id DM acc F0 F1 F2 S/N" + chr(10) +
            f"0 {dm:.16g} {a:.16g} {f0:.16g} {f1:.16g} {f2:.16g} {snr:.16g}" + chr(10),
            encoding="utf-8",
        )
        args = [binary, "--template", str(template.resolve()),
                "--candfile", str(cand.resolve()), "--pepoch", f"{mid_epoch:.13f}",
                "--nosearch", "-b", str(fold.get("nbin", 64)),
                "-o", str((fold_dir / p["ID"]).resolve()),
                "-f", str(output_fil.resolve())]
        with (fold_dir / "fold.log").open("w") as output:
            subprocess.run(args, check=True, stdout=output, stderr=subprocess.STDOUT)
        products = [str(f.relative_to(run_dir)) for f in fold_dir.rglob("*") if f.is_file()]
        results.append({"pulsar": p["ID"], "command": args, "products": products})
    return results


def run(config_path, *, dry_run=False, selected_case=None):
    cfg = json.loads(Path(config_path).read_text(encoding="utf-8"))
    cases = plan_cases(cfg)
    if selected_case is not None:
        cases = [c for c in cases if c[0] == selected_case]
        if not cases:
            raise ValueError(f"Unknown case {selected_case!r}")
    base = Path(cfg["filterbank"]).expanduser().resolve()
    out_root = Path(cfg["output_dir"]).expanduser().resolve()
    for name, pulsars in cases:
        print(f"[INCEPTION-LITE] {name}: " +
              ", ".join(f"{p['ID']} P0={p.get('P0','-')} DM={p.get('DM',0)} A0={p.get('A0',0)} A1={p.get('A1',0)}" for p in pulsars))
    if dry_run:
        print(f"[INCEPTION-LITE] dry run only: {len(cases)} planned filterbank(s)")
        return cases
    if not base.is_file():
        raise ValueError(f"Input filterbank not found: {base}")
    header = read_filterbank_header(base)
    if cfg.get("fold", {}).get("enabled", False):
        fold = cfg["fold"]
        if not Path(fold.get("template", "")).is_file():
            raise ValueError("Cannot start: fold.template is missing")
        if shutil.which(fold.get("binary", "psrfold_fil")) is None:
            raise ValueError("Cannot start: PulsarX psrfold_fil is unavailable")
        # Validate before doing expensive injections.
        if any(p.get("frame", cfg.get("psr_global", {}).get("frame", "bary")) != "topo"
               for _, ps in cases for p in ps):
            raise ValueError("Automatic folding currently supports topocentric injections only")
    for name, pulsars in cases:
        dest = out_root / name
        if dest.exists():
            raise FileExistsError(f"Refusing to overwrite existing run: {dest}")
    out_root.mkdir(parents=True, exist_ok=True)
    records = []
    for name, pulsars in cases:
        dest = out_root / name
        dest.mkdir()
        resolved = {
            "psr_global": {
                **copy.deepcopy(cfg.get("psr_global", {})),
                "injection_id": f"{cfg.get('id','inception_lite')}_{name}",
                "global_seed": case_seed(cfg["seed"], name),
            },
            "pulsars": pulsars,
        }
        (dest / "resolved_injection.json").write_text(
            json.dumps(resolved, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        cmd = [sys.executable, str(INJECTOR), "--signal",
               str((dest / "resolved_injection.json").resolve()),
               "--fb", str(base), "--output", str(dest), "--ncpu",
               str(cfg.get("ncpu", 1)), "--gulp_size_GB",
               str(cfg.get("gulp_size_GB", 0.01)), "--stats_samples",
               str(cfg.get("stats_samples", 100000))]
        if cfg.get("offline", True):
            cmd.append("--offline")
        if cfg.get("cache"):
            cmd += ["--cache", str(Path(cfg["cache"]).expanduser().resolve())]
        if cfg.get("ephem"):
            cmd += ["--ephem", str(cfg["ephem"])]
        if cfg.get("generate_fb") is not None:
            gen = cfg["generate_fb"]
            if (not isinstance(gen, list) or len(gen) != 2 or
                    not all(isinstance(x, (int, float)) for x in gen) or
                    not math.isfinite(gen[0]) or gen[1] <= 0):
                raise ValueError("generate_fb must be [mean, positive_std]")
            cmd += ["--generate_fb", f"{gen[0]},{gen[1]}"]
        with (dest / "inject.log").open("w", encoding="utf-8") as log:
            subprocess.run(cmd, check=True, cwd=REPO, stdout=log, stderr=subprocess.STDOUT)
        output_fil = dest / f"{base.stem}_{resolved['psr_global']['injection_id']}.fil"
        if not output_fil.is_file():
            raise RuntimeError(f"Injector did not create expected output: {output_fil}")
        actual = read_filterbank_header(output_fil)
        for key in ("tstart", "tsamp", "nchans", "nbits", "actual_nsamples"):
            if actual[key] != header[key]:
                raise RuntimeError(f"Output changed filterbank {key}: {actual[key]} vs {header[key]}")
        report = list(dest.glob("report_*.json"))
        if len(report) != 1:
            raise RuntimeError(f"Expected exactly one injection report, got {len(report)}")
        folding = fold_truth(dest, output_fil, pulsars, cfg, header)
        record = {
            "status": "success", "case": name, "input_fil": str(base),
            "output_fil": str(output_fil), "report": str(report[0]),
            "injection_id": resolved["psr_global"]["injection_id"],
            "seed": resolved["psr_global"]["global_seed"],
            "pulsars": pulsars, "header": header,
            "inject_command": cmd, "folding": folding,
            "provenance": {
                "python": sys.version.split()[0],
                "injector_script_sha256": hashlib.sha256(INJECTOR.read_bytes()).hexdigest(),
                "pulsar_model_sha256": hashlib.sha256(
                    (REPO / "injector" / "pulsar_model.py").read_bytes()).hexdigest(),
                "software": {
                    package: importlib.metadata.version(package)
                    for package in ("numpy", "scipy", "astropy", "sympy", "pandas")
                },
            },
        }
        (dest / "run_manifest.json").write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        records.append(record)
        print(f"[INCEPTION-LITE] success: {output_fil}")
    return records


def main():
    p = argparse.ArgumentParser(description="Nextflow-free INCEPTION injection runner")
    p.add_argument("--config", required=True, help="Portable run JSON configuration")
    p.add_argument("--dry-run", action="store_true", help="Plan only; no file changes or dependencies")
    p.add_argument("--case", help="Only run one case_0001 etc (useful for SLURM arrays)")
    args = p.parse_args()
    try:
        run(args.config, dry_run=args.dry_run, selected_case=args.case)
    except (ValueError, KeyError, FileNotFoundError, FileExistsError, subprocess.CalledProcessError) as exc:
        p.exit(2, f"[INCEPTION-LITE] ERROR: {exc}\n")


if __name__ == "__main__":
    main()
