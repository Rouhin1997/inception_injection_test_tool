"""Deliberately extreme jerk-injection falsification test (no PULSEJET involved).

This is not an astrophysical population model: high jerk creates measurable
phase drift in a tiny CI filterbank. Physical phase is computed independently
from the first-order INCEPTION acceleration-jerk convention.
"""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from inception_lite.make_filterbank import make_filterbank
from inception_lite.runner import read_filterbank_header, run


def samples(path):
    h = read_filterbank_header(path)
    return np.fromfile(path, dtype=np.uint8, offset=h["data_offset"]).reshape(
        h["actual_nsamples"], h["nchans"]).astype(float)


class StressJerk(unittest.TestCase):
    def test_injected_waveform_matches_jerk_sign_and_midpoint(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base = make_filterbank(root / "noise.fil", samples=8192,
                                   nchans=8, tsamp=.001, seed=407)
            cfg = {
                "id": "stress", "seed": 1234, "mode": "together",
                "filterbank": str(base), "output_dir": str(root/"outputs"),
                "ncpu": 1, "offline": True, "gulp_size_GB": .001,
                "stats_samples": 4000,
                "psr_global": {"frame": "topo", "ACCEPOCH": 0.5,
                               "create_parfile": "pulsarx"},
                "pulsars": [{"ID": "JERK", "P0": .025, "DM": 0,
                             "SNR": 140, "A0": 0, "A1": -1000000,
                             "duty_cycle": .08}],
            }
            config_path = root/"config.json"
            config_path.write_text(json.dumps(cfg))
            output = Path(run(config_path)[0]["output_fil"])
            signal = (samples(output)-samples(base)).mean(axis=1)
            t = np.arange(len(signal))*.001
            tau = t - len(t)*.001/2
            f0 = 40.
            c = 299792458.
            phase_true = f0*(tau - cfg["pulsars"][0]["A1"]*tau**3/(6*c))
            phase_wrong = f0*(tau + cfg["pulsars"][0]["A1"]*tau**3/(6*c))
            def first_harmonic(phase):
                return abs(np.mean(signal*np.exp(-2j*np.pi*phase)))
            recovered = first_harmonic(phase_true)
            uncorrected = first_harmonic(f0*t)
            wrong = first_harmonic(phase_wrong)
            print(f"[JERK VALIDATION] true={recovered:.4g} uncorrected={uncorrected:.4g} reverse_sign={wrong:.4g}")
            self.assertGreater(recovered, 1.7*uncorrected)
            self.assertGreater(recovered, 1.7*wrong)


if __name__ == "__main__":
    unittest.main()
