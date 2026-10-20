"""End-to-end multiple-pulsar, generated-noise, multi-CPU injection smoke."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from inception_lite.make_filterbank import make_filterbank
from inception_lite.runner import read_filterbank_header, run


class MultiplePulsars(unittest.TestCase):
    def test_two_pulsars_one_filterbank_parallel_workers(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            infile = make_filterbank(root/"noise.fil", samples=8192, nchans=8,
                                     tsamp=.001, seed=941)
            cfg = {
                "id": "multiple", "seed": 424242, "mode": "together",
                "filterbank": str(infile), "output_dir": str(root/"results"),
                "ncpu": 2, "gulp_size_GB": .001, "stats_samples": 1000,
                "offline": True, "generate_fb": [128, 6],
                "psr_global": {"frame": "topo", "ACCEPOCH": .5,
                               "create_parfile": "pulsarx"},
                "pulsars": [
                    {"ID": "PSR_ONE", "P0": .02, "DM": 0, "SNR": 60, "A0": 0, "A1": 0},
                    {"ID": "PSR_TWO", "P0": .03, "DM": 0, "SNR": 60, "A0": 100, "A1": -2},
                ],
            }
            cfgpath = root/"plan.json"
            cfgpath.write_text(json.dumps(cfg))
            result = run(cfgpath)
            self.assertEqual(len(result), 1)
            output = Path(result[0]["output_fil"])
            header = read_filterbank_header(output)
            self.assertEqual(header["actual_nsamples"], 8192)
            files = {p.name for p in output.parent.glob("*.candfile")}
            self.assertEqual(files, {"PSR_ONE.candfile", "PSR_TWO.candfile"})
            values = np.fromfile(output, dtype=np.uint8, offset=header["data_offset"])
            self.assertGreater(float(np.std(values)), 2)
            self.assertLess(float(np.std(values)), 40)


if __name__ == "__main__":
    unittest.main()
