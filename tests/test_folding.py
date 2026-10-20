"""Test optional PulsarX handoff without needing a real pulsar installation."""
import json
import math
import tempfile
import unittest
from pathlib import Path

from inception_lite.make_filterbank import make_filterbank
from inception_lite.runner import fold_truth, read_filterbank_header


class TruthFolding(unittest.TestCase):
    def test_fold_handoff_keeps_jerk_and_midpoint_epoch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fil = make_filterbank(root / "input.fil", samples=1024, nchans=8, tsamp=.001)
            header = read_filterbank_header(fil)
            template = root / "fold.template"
            template.write_text("# stub fold template")
            binary = root / "fake_psrfold_fil"
            binary.write_text(
                "#!/usr/bin/env python3\n"
                "import pathlib,sys\n"
                "args=sys.argv\n"
                "c=pathlib.Path(args[args.index('--candfile')+1])\n"
                "assert len(c.read_text().splitlines())==2\n"
                "p=pathlib.Path(args[args.index('-o')+1])\n"
                "p.with_suffix('.png').write_bytes(b'fake fold')\n"
                "p.with_suffix('.ar').write_bytes(b'fake archive')\n"
            )
            binary.chmod(0o755)
            cfg = {"psr_global": {"frame": "topo", "ACCEPOCH": .5},
                   "fold": {"enabled": True, "template": str(template),
                            "binary": str(binary), "nbin": 64}}
            pulsar = {"ID": "TEST", "P0": .02, "DM": 33, "SNR": 40,
                      "A0": 300, "A1": -4}
            results = fold_truth(root, fil, [pulsar], cfg, header)
            self.assertEqual(len(results), 1)
            products = list((root / "folds" / "TEST").glob("*"))
            self.assertTrue(any(x.suffix == ".png" for x in products))
            self.assertTrue(any(x.suffix == ".ar" for x in products))
            lines = (root / "folds" / "TEST" / "TEST_phase_truth.cands").read_text().splitlines()
            self.assertEqual(len(lines), 2)
            fields = lines[1].split()
            c = 299792458.
            self.assertAlmostEqual(float(fields[3]), 50.)
            self.assertAlmostEqual(float(fields[4]), -300*50/c, places=13)
            self.assertAlmostEqual(float(fields[5]), 4*50/c, places=13)
            self.assertIn("--nosearch", results[0]["command"])
            self.assertIn("--pepoch", results[0]["command"])
            self.assertIn("--candfile", results[0]["command"])

    def test_reject_barycentric_shortcuts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fil = make_filterbank(root / "input.fil", samples=128)
            template = root / "template"
            template.write_text("x")
            cfg = {"fold": {"enabled": True, "template": str(template)},
                   "psr_global": {"frame": "bary"}}
            pulsar = {"ID": "TEST", "P0": .02, "SNR": 20}
            with self.assertRaisesRegex(ValueError, "topocentric"):
                fold_truth(root, fil, [pulsar], cfg, read_filterbank_header(fil))


if __name__ == "__main__":
    unittest.main()
