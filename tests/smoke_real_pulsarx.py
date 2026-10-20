"""Real PulsarX truth-fold smoke, run INSIDE the combined injection+fold image.

This test is NOT part of the lightweight Python-only GitHub CI. It runs only
where the actual PulsarX native binary and fold template have been installed.
"""
from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path

from inception_lite.make_filterbank import make_filterbank
from inception_lite.runner import run


class RealPulsarX(unittest.TestCase):
    def test_inject_then_fold_to_real_png_and_archive(self):
        binary = shutil.which("psrfold_fil")
        template = Path("/usr/local/share/pulsarx/meerkat_fold.template")
        self.assertIsNotNone(binary, "PulsarX psrfold_fil must be installed")
        self.assertTrue(template.is_file(), "PulsarX fold template must be available")
        with tempfile.TemporaryDirectory() as scratch:
            base = Path(scratch)
            fil = make_filterbank(
                base / "noise.fil", samples=16384,
                nchans=8, tsamp=0.001, seed=20261021
            )
            cfg = {
                "id": "real_fold",
                "seed": 20261021,
                "mode": "together",
                "filterbank": str(fil),
                "output_dir": str(base / "results"),
                "ncpu": 1,
                "offline": True,
                "gulp_size_GB": 0.001,
                "stats_samples": 3000,
                "psr_global": {
                    "frame": "topo",
                    "ACCEPOCH": 0.5,
                    "create_parfile": "pulsarx",
                },
                "pulsars": [
                    {
                        "ID": "TEST_AJ",
                        "P0": 0.02,
                        "DM": 0,
                        "SNR": 100,
                        "A0": 400,
                        "A1": -3,
                        "duty_cycle": 0.1,
                    }
                ],
                "fold": {
                    "enabled": True,
                    "binary": binary,
                    "template": str(template),
                    "nbin": 32,
                    "nsubband": 8,
                },
            }
            cfg_path = base / "plan.json"
            cfg_path.write_text(json.dumps(cfg))
            try:
                results = run(cfg_path)
            except Exception:
                for f in (base / "results").rglob("*.log"):
                    print("\n[REAL FOLD FAILURE] ", f)
                    print(f.read_text(errors="replace")[-8000:])
                raise
            self.assertEqual(len(results), 1)
            work = Path(results[0]["output_fil"]).parent
            self.assertTrue((work / "TEST_AJ.candfile").is_file())
            self.assertTrue(results[0]["folding"], "No real folding results")
            fold_dir = work / "folds" / "TEST_AJ"
            pngs = list(fold_dir.glob("*.png"))
            archives = list(fold_dir.glob("*.ar"))
            print("[REAL PULSARX FOLD] figures:", [x.name for x in pngs])
            print("[REAL PULSARX FOLD] archives:", [x.name for x in archives])
            self.assertGreater(len(pngs), 0, "PulsarX did not generate diagnostic PNG")
            self.assertGreater(len(archives), 0, "PulsarX did not generate a folded archive")
            self.assertGreater(pngs[0].stat().st_size, 1000)
            self.assertGreater(archives[0].stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
