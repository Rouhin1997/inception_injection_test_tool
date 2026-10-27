"""Static and planning-regression tests for the M62 UCB100 campaign.

No actual full-observation injections; uses placeholder test files and a
mocked SIGPROC header reader to validate deterministic planning quickly.
"""
import csv
import importlib.util
import json
import math
import tempfile
import unittest
from collections import Counter
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "prepare_m62_ucb100.py"
SPEC = importlib.util.spec_from_file_location("prepare_m62_ucb100", SCRIPT)
PLAN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PLAN)


class UCB100PlanTests(unittest.TestCase):
    def test_physical_midpoint_motion_and_strata(self):
        pop = PLAN.build_population()
        self.assertEqual(len(pop), 100)
        counts = Counter(s["motion_class"] for s in pop)
        self.assertEqual(counts, {
            "isolated": 10,
            "j_000_025": 15,
            "j_025_060": 20,
            "j_060_120": 25,
            "j_120_210": 30})
        self.assertEqual(pop, PLAN.build_population(), "Seed must reproduce catalogue")
        for p in pop:
            if p["pb_min"] is None:
                self.assertEqual(p["acc_mid_mps2"], 0)
                self.assertEqual(p["jerk_mid_mps3"], 0)
            else:
                self.assertGreaterEqual(p["pb_min"], 25)
                self.assertLessEqual(p["pb_min"], 75)
                self.assertLessEqual(p["amax_mps2"], 550)
                self.assertLessEqual(p["jmax_mps3"], 2.5)
                x, a, j, am, jm = PLAN.orbital_params(
                    p["pb_min"], p["mc_msun"], p["inc_deg"], p["phase_mid_rad"])
                self.assertAlmostEqual(a, p["acc_mid_mps2"])
                self.assertAlmostEqual(j, p["jerk_mid_mps3"])
                self.assertAlmostEqual(x, p["x_lt_s"])
        for i in range(20):
            self.assertEqual(len(pop[i*5:(i+1)*5]), 5)

    def test_generate_20_safe_configs_with_true_orbital_t0(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            base = root / "stage6_m62" / "m62_15m_windows"
            base.mkdir(parents=True)
            for i in range(1,5):
                (base / ("W%d_M62_background.fil" % i)).write_bytes(b"placeholder")
            header = dict(nbits=8,nchans=256,actual_nsamples=14835712,
                          tsamp=0.00006024,tstart=60200.12345)
            # Reproduce the actual frontend2 failure sequence: both Slurm
            # pilots exited before planning and left job logs in an otherwise
            # empty campaign directory. Keep those logs, then create config.
            dest=root/"stage10_ucb100"
            (dest/"logs").mkdir(parents=True)
            old_job_log=dest/"logs"/"pilot_node1_66570_1.out"
            old_job_log.write_text("Missing injection SIF/config\\n")
            with patch.object(PLAN, "read_filterbank_header", return_value=header):
                rows=PLAN.generate(root,dest)
            self.assertEqual(old_job_log.read_text(),"Missing injection SIF/config\\n")
            self.assertEqual(len(rows),100)
            self.assertEqual(len(list((dest/"configs").glob("batch_*.json"))),20)
            with (dest/"ucb100_truth.csv").open() as fh:
                csvrows=list(csv.DictReader(fh))
            self.assertEqual(len(csvrows),100)
            for i in range(1,21):
                config=json.loads((dest/"configs"/("batch_%02d.json" % i)).read_text())
                self.assertEqual(config["mode"],"together")
                self.assertEqual(len(config["pulsars"]),5)
                self.assertFalse(config["fold"]["enabled"])
                self.assertEqual(config["psr_global"]["frame"],"topo")
                for pulsar in config["pulsars"]:
                    self.assertNotIn("A0",pulsar)
                    self.assertNotIn("A1",pulsar)
                    if "binary_period" in pulsar:
                        self.assertGreater(pulsar["T0"],50000)
                        self.assertEqual(pulsar["ecc"],0)
            for row in rows:
                if row["pb_min"] is None:
                    self.assertIsNone(row["t0_mjd"])
                else:
                    phase=2*math.pi/(row["pb_min"]*60.0)
                    model_phase=phase*((header["tstart"]-
                                      row["t0_mjd"])*86400 +
                                      header["actual_nsamples"]*header["tsamp"]/2)
                    self.assertAlmostEqual(model_phase,row["phase_mid_rad"],places=4)
            with patch.object(PLAN,"read_filterbank_header",return_value=header):
                with self.assertRaises(FileExistsError):
                    PLAN.generate(root,dest)
            for i in range(1,5):
                self.assertEqual(
                  (base / ("W%d_M62_background.fil" % i)).read_bytes(),
                  b"placeholder")
            # Only logs/ may pre-exist; unrelated files or a partial old
            # catalogue must still cause a hard failure.
            unsafe=root/"stage10_wrong_partial"
            unsafe.mkdir()
            (unsafe/"configs").mkdir()
            with patch.object(PLAN,"read_filterbank_header",return_value=header):
                with self.assertRaises(FileExistsError):
                    PLAN.generate(root,unsafe)



if __name__=="__main__":
    unittest.main()
