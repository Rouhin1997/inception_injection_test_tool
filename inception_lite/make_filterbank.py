"""Small deterministic 8-bit SIGPROC filterbank template for injector tests.

This creates noise-only data. The real signal is injected by INCEPTION itself.
It does not replace independently validated physical injection models.
"""
from __future__ import annotations
import argparse
import random
import struct
from pathlib import Path


def make_filterbank(path, *, samples=32768, nchans=8, tsamp=0.001,
                    tstart=60000.0, fch1=1400.0, foff=-1.0,
                    mean=128.0, std=6.0, seed=771):
    if samples <= 0 or nchans < 2 or tsamp <= 0 or std <= 0 or seed < 0:
        raise ValueError("samples>0, nchans>=2, tsamp>0, std>0 and seed>=0 required")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite {path}")
    rng = random.Random(seed)
    def field(f, value):
        b = value.encode("ascii")
        f.write(struct.pack("<i", len(b)) + b)
    with path.open("wb") as f:
        field(f, "HEADER_START")
        for key, value in [
            ("source_name", "INCEPTION_LITE_TEST"),
            ("telescope_id", 64),  # MeerKAT: supported by INCEPTION's site table
            ("machine_id", 1), ("data_type", 1),
            ("tstart", tstart), ("tsamp", tsamp), ("fch1", fch1), ("foff", foff),
            ("src_raj", 123456.0), ("src_dej", -300000.0),
            ("nchans", nchans), ("nbits", 8), ("nifs", 1),
        ]:
            field(f, key)
            if isinstance(value, str):
                field(f, value)
            elif isinstance(value, int):
                f.write(struct.pack("<i", value))
            else:
                f.write(struct.pack("<d", float(value)))
        field(f, "HEADER_END")
        for _ in range(samples):
            f.write(bytes(max(0, min(255, round(rng.gauss(mean, std))))
                          for _ in range(nchans)))
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True)
    p.add_argument("--samples", type=int, default=32768)
    p.add_argument("--nchans", type=int, default=8)
    p.add_argument("--tsamp", type=float, default=0.001)
    p.add_argument("--tstart", type=float, default=60000)
    p.add_argument("--seed", type=int, default=771)
    p.add_argument("--mean", type=float, default=128)
    p.add_argument("--std", type=float, default=6)
    p.add_argument("--fch1", type=float, default=1400)
    p.add_argument("--foff", type=float, default=-1.0)
    args = p.parse_args()
    path = make_filterbank(args.output, samples=args.samples, nchans=args.nchans,
                           tsamp=args.tsamp, tstart=args.tstart, seed=args.seed,
                           mean=args.mean, std=args.std,
                           fch1=args.fch1, foff=args.foff)
    print("Created noise-only SIGPROC template:", path)


if __name__ == "__main__":
    main()
