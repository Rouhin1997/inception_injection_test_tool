"""End-to-end CPU injection smoke: real INCEPTION code, real small .fil files.

No GPU or Nextflow. Tests that the output has intact metadata, its data
change from noise to pulses, its zero-motion pulse is coherent at known P0,
and seed-reproduced runs give byte-identical injected filterbanks.
"""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from inception_lite.make_filterbank import make_filterbank
from inception_lite.runner import read_filterbank_header, run


def data(path):
    header = read_filterbank_header(path)
    raw = np.fromfile(path, dtype=np.uint8, offset=header["data_offset"])
    return raw.reshape(header["actual_nsamples"], header["nchans"])


class EndToEnd(unittest.TestCase):
    def test_two_real_injections_and_reproducibility(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base = make_filterbank(root / "noise.fil", samples=8192, nchans=8,
                                   tsamp=0.001, seed=303)
            cfg = {
                "id": "smoke", "seed": 20261020,
                "mode": "separate", "filterbank": str(base),
                "output_dir": str(root / "run_a"), "ncpu": 1,
                "offline": True, "gulp_size_GB": 0.001,
                "stats_samples": 3000,
                "psr_global": {"frame": "topo", "ACCEPOCH": 0.5,
                               "create_parfile": "pulsarx"},
                "pulsars": [
                    {"ID": "STEADY", "P0": 0.02, "SNR": 90,
                     "DM": 0, "A0": 0, "A1": 0, "duty_cycle": 0.1},
                    {"ID": "MOVING", "P0": 0.025, "SNR": 90,
                     "DM": 0, "A0": 2000, "A1": -30, "duty_cycle": 0.1},
                ],
                "fold": {"enabled": False},
            }
            plan = root / "config.json"
            plan.write_text(json.dumps(cfg))
            results_a = run(plan)
            self.assertEqual(len(results_a), 2)
            baseline = data(base).astype(float)
            steady = data(Path(results_a[0]["output_fil"])).astype(float)
            motion = data(Path(results_a[1]["output_fil"])).astype(float)
            self.assertEqual(steady.shape, baseline.shape)
            self.assertGreater(np.mean(steady != baseline), 0.01)
            self.assertGreater(np.mean(motion != baseline), 0.01)
            self.assertNotEqual(hashlib.sha256(steady.tobytes()).digest(),
                                hashlib.sha256(motion.tobytes()).digest())
            for result in results_a:
                manifest = json.loads((Path(result["output_fil"]).parent / "run_manifest.json").read_text())
                self.assertEqual(manifest["status"], "success")
                self.assertTrue(list(Path(result["output_fil"]).parent.glob("*.candfile")))
                self.assertTrue(Path(manifest["report"]).exists())

            # Fold the injected zero-motion signal at its known 20 ms period;
            # a correctly injected signal must have a visible phase-dependent
            # excess above noise when summed over channels.
            mean_signal = np.mean(steady-baseline, axis=1)
            phase = np.arange(mean_signal.size) % 20
            profile = np.array([mean_signal[phase == k].mean() for k in range(20)])
            print("[INJECTION] coherent profile amplitude", float(np.ptp(profile)))
            self.assertGreater(np.ptp(profile), 0.5)
            self.assertGreater(float(profile.max()), 0.5)

            cfg["output_dir"] = str(root / "run_b")
            plan.write_text(json.dumps(cfg))
            results_b = run(plan)
            for a, b in zip(results_a, results_b):
                self.assertEqual(Path(a["output_fil"]).read_bytes(),
                                 Path(b["output_fil"]).read_bytes())
                self.assertEqual(a["seed"], b["seed"])


if __name__ == "__main__":
    unittest.main()
