import copy
import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from inception_lite.make_filterbank import make_filterbank
from inception_lite.runner import (
    apply_folding_override, case_seed, check_pulsars, draw, main, plan_cases,
    read_filterbank_header, run,
)


def config():
    return {
        "id": "aj_test",
        "seed": 20261020,
        "filterbank": "/not/required/for/dry-run.fil",
        "output_dir": "/not/required/for/dry-run",
        "psr_global": {"frame": "topo", "ACCEPOCH": 0.5},
        "pulsars": [{"ID": "PSR_TEST", "P0": 0.02, "DM": 40,
                     "SNR": 30, "A0": 20, "A1": 1.2}],
    }


class Planning(unittest.TestCase):
    def test_multiple_same_filterbank(self):
        x = config()
        x["pulsars"].append({"ID": "PSR_2", "F0": 80, "DM": 9, "SNR": 20})
        planned = plan_cases(x)
        self.assertEqual(len(planned), 1)
        self.assertEqual(len(planned[0][1]), 2)

    def test_separate(self):
        x = config()
        x["pulsars"].append({"ID": "PSR_2", "P0": .04, "SNR": 20})
        x["mode"] = "separate"
        self.assertEqual(len(plan_cases(x)), 2)

    def test_parameter_grid(self):
        x = config()
        x["mode"] = "grid"
        x["grid"] = {"A0": [-10, 0, 10], "A1": [-1, 1]}
        cases = plan_cases(x)
        self.assertEqual(len(cases), 6)
        self.assertEqual([p[0]["A0"] for _, p in cases], [-10,-10,0,0,10,10])
        self.assertEqual(len({p[0]["ID"] for _, p in cases}), 6)

    def test_seeded_random(self):
        x = config()
        x["mode"] = "random"
        x["random"] = {
            "count": 7,
            "parameters": {"P0": {"distribution": "loguniform", "min": .005, "max": .3},
                           "A1": {"distribution": "choice", "values": [-3, 0, 3]},
                           "DM": {"distribution": "uniform", "min": 0, "max": 30}},
        }
        self.assertEqual(plan_cases(x), plan_cases(x))
        vals = [p[0] for _, p in plan_cases(x)]
        self.assertTrue(all(.005 <= p["P0"] <= .3 for p in vals))
        self.assertTrue(all(p["A1"] in (-3, 0, 3) for p in vals))
        self.assertNotEqual(case_seed(x["seed"], "case_0001"),
                            case_seed(x["seed"], "case_0002"))

    def test_untrusted_output_ids_rejected(self):
        x = config()
        x["pulsars"][0]["ID"] = "../unsafe"
        with self.assertRaisesRegex(ValueError, "Unsafe"):
            plan_cases(x)

    def test_nonfinite_motion_rejected(self):
        x = config()
        x["pulsars"][0]["A1"] = float("nan")
        with self.assertRaisesRegex(ValueError, "finite"):
            plan_cases(x)

    def test_no_implicit_unseeded_runs(self):
        x = config()
        del x["seed"]
        with self.assertRaises(ValueError):
            plan_cases(x)

    def test_prevent_double_midpoint_epoch_shift(self):
        x = config()
        x["psr_global"]["PEPOCH"] = 0.5
        with self.assertRaisesRegex(ValueError, "PEPOCH"):
            plan_cases(x)

    def test_require_explicit_topocentric_or_barycentric_frame(self):
        x = config()
        del x["psr_global"]["frame"]
        with self.assertRaisesRegex(ValueError, "frame explicitly"):
            plan_cases(x)

    def test_grid_case_limit(self):
        x = config()
        x["mode"] = "grid"
        x["max_cases"] = 1
        x["grid"] = {"A0": [1, 2, 3]}
        with self.assertRaisesRegex(ValueError, "Too many"):
            plan_cases(x)

    def test_dry_run_no_files_needed(self):
        x = config()
        x["mode"] = "grid"
        x["grid"] = {"A0": [-1, 0, 1]}
        with tempfile.TemporaryDirectory() as tmp:
            x["output_dir"] = str(Path(tmp) / "out")
            cfg = Path(tmp) / "plan.json"
            cfg.write_text(json.dumps(x))
            result = run(cfg, dry_run=True)
            self.assertEqual(len(result), 3)
            self.assertFalse(Path(x["output_dir"]).exists())

    def test_folding_cli_precedence_and_config_preservation(self):
        x = config()
        x["fold"] = {"enabled": False, "template": "my_fold.template", "nbin": 64}
        yes = apply_folding_override(copy.deepcopy(x), True)
        self.assertTrue(yes["fold"]["enabled"])
        self.assertEqual(yes["fold"]["template"], "my_fold.template")
        self.assertFalse(x["fold"]["enabled"])
        x["fold"]["enabled"] = True
        no = apply_folding_override(copy.deepcopy(x), False)
        self.assertFalse(no["fold"]["enabled"])
        default = apply_folding_override(copy.deepcopy(x), None)
        self.assertTrue(default["fold"]["enabled"])
        with self.assertRaises(ValueError):
            apply_folding_override(x, "yes")

    def test_cli_yes_no_flag_and_underscore_alias(self):
        for arg, value, expected in (
            ("--perform-folds", "yes", True),
            ("--perform-folds", "no", False),
            ("--perform_folds", "yes", True),
        ):
            with self.subTest(arg=arg, value=value):
                with mock.patch("sys.argv", ["inception-lite", "--config", "dummy.json", arg, value]):
                    with mock.patch("inception_lite.runner.run") as mocked:
                        main()
                        self.assertEqual(mocked.call_args.kwargs["perform_folds"], expected)
        with mock.patch("sys.argv", ["inception-lite", "--config", "dummy.json"]):
            with mock.patch("inception_lite.runner.run") as mocked:
                main()
                self.assertIsNone(mocked.call_args.kwargs["perform_folds"])

    def test_invalid_folding_cli_choice_is_rejected(self):
        with mock.patch("sys.argv", ["inception-lite", "--config", "dummy.json",
                                     "--perform-folds", "maybe"]):
            with self.assertRaises(SystemExit):
                main()

    def test_filterbank_header_and_reproducibility(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp) / "a.fil", Path(tmp) / "b.fil"
            make_filterbank(a, samples=256, nchans=4, seed=333)
            make_filterbank(b, samples=256, nchans=4, seed=333)
            self.assertEqual(a.read_bytes(), b.read_bytes())
            h = read_filterbank_header(a)
            self.assertEqual(h["actual_nsamples"], 256)
            self.assertEqual(h["nchans"], 4)
            self.assertEqual(h["telescope_id"], 64)

    def test_output_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            t = Path(tmp)
            base = make_filterbank(t / "in.fil", samples=32)
            dest = t / "out" / "case_0001"
            dest.mkdir(parents=True)
            x = config()
            x["filterbank"] = str(base)
            x["output_dir"] = str(t / "out")
            c = t / "plan.json"
            c.write_text(json.dumps(x))
            with self.assertRaises(FileExistsError):
                run(c)


class Physics(unittest.TestCase):
    def test_phase_and_frequency_match_first_order_doppler(self):
        from injector.pulsar_model import PulsarModel
        from scipy.constants import c
        p = PulsarModel.__new__(PulsarModel)
        p.FX_list = [40.0]
        p.AX_list = [400.0, -3.0]
        p.spin_ref = 0.0
        p.accepoch = 500.0
        p.ephemeris_freq = 0.0
        p.get_spin_functions({"phase_offset": 0.0})
        for t in (0., 100., 500., 1000.):
            tau = t - 500
            phase = 40 * (tau - 400 * tau**2/(2*c) + 3*tau**3/(6*c))
            f = 40 * (1 - (400 * tau - 1.5*tau**2)/c)
            self.assertAlmostEqual(p.phase_func(t), phase, places=8)
            self.assertAlmostEqual(p.spin_func(0.0), 40., places=8)
            h = 0.005
            df = (p.phase_func(t+h)-p.phase_func(t-h))/(2*h)
            self.assertAlmostEqual(df, f, delta=1e-6)
        self.assertAlmostEqual(p.FX_doppler[1], -40 * 400 / c, delta=1e-12)
        self.assertAlmostEqual(p.FX_doppler[2], 40 * 3 / c, delta=1e-12)

    def test_midpoint_model_versus_arrival_time_inverse(self):
        # Independent inverse of PULSEJET's midpoint arrival equation.
        # Difference begins at second order in v/c, not a sign or epoch error.
        c = 299792458.0
        f0 = 40.
        a, jerk = 400., -3.
        for tau in (-500., -250., 0., 200., 500.):
            exact_emission = tau
            for _ in range(8):
                delay = (a * exact_emission**2 / 2 +
                         jerk * exact_emission**3 / 6) / c
                derivative = 1 + (a*exact_emission + jerk*exact_emission**2/2)/c
                exact_emission -= (exact_emission+delay-tau)/derivative
            linear_doppler = tau - (a*tau**2/2+jerk*tau**3/6)/c
            # Relative frequency and time offsets differ by O((v/c)^2).
            self.assertLess(abs(exact_emission-linear_doppler), 0.03)
            self.assertLess(abs(f0*(exact_emission-linear_doppler)), 1.1)


if __name__ == "__main__":
    unittest.main()
