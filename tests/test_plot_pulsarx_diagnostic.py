"""Unit tests for injection-truth annotation on PulsarX diagnostic plots.

These tests ensure source metadata comes from the manifest, not guesses from
archive names or arbitrary labels; FITS plot content is tested on real data.
"""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "plot_pulsarx_diagnostic.py"
spec = importlib.util.spec_from_file_location("plot_pulsarx_diagnostic", SCRIPT)
plotter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plotter)


class TestManifestTruth(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.case = Path(self.temp.name) / "case_0001"
        self.fold = self.case / "folds" / "M62_W1_I01"
        self.fold.mkdir(parents=True)
        self.archive = self.fold / "M62_W1_I01_dummy.ar"
        self.archive.touch()
        self.manifest = self.case / "run_manifest.json"
        self.manifest.write_text(json.dumps({
            "status": "success",
            "pulsars": [
                {"ID": "M62_W1_I01", "P0": 1 / 247, "DM": 0,
                 "A0": 0, "A1": 0, "SNR": 100, "duty_cycle": 0.05},
                {"ID": "M62_W1_I02", "P0": 1 / 263, "DM": 4,
                 "A0": 40, "A1": 0, "SNR": 60, "duty_cycle": 0.05}
            ],
        }))
        (self.case / "resolved_injection.json").write_text(json.dumps({
            "psr_global": {"frame": "topo", "ACCEPOCH": 0.5}
        }))

    def test_truth_metadata(self):
        t = plotter.read_injection_truth(self.archive, "M62_W1_I01")
        self.assertEqual(t["id"], "M62_W1_I01")
        self.assertAlmostEqual(t["frequency"], 247)
        self.assertAlmostEqual(t["period"], 1 / 247)
        self.assertEqual(t["dm"], 0)
        self.assertEqual(t["a"], 0)
        self.assertEqual(t["jerk"], 0)
        self.assertEqual(t["snr"], 100)
        self.assertEqual(t["duty"], 0.05)
        self.assertEqual(t["frame"], "topo")
        self.assertEqual(t["accep"], 0.5)

    def test_refuse_other_sources_truth(self):
        with self.assertRaises(ValueError):
            plotter.read_injection_truth(self.archive, "DOES_NOT_EXIST")

    def test_missing_manifest_is_not_fabricated(self):
        self.manifest.unlink()
        self.assertIsNone(plotter.read_injection_truth(self.archive, "M62_W1_I01"))
        with self.assertRaises(FileNotFoundError):
            plotter.read_injection_truth(self.archive, "M62_W1_I01", self.manifest)


if __name__ == "__main__":
    unittest.main()
