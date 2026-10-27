"""Regression tests: PulsarX channel dispersion shifts must be undone before plotting."""
import unittest

import numpy as np

from inception_lite.fold_plots import restore_pulsarx_dispersion


class TestPulsarXArchiveAlignment(unittest.TestCase):
    @staticmethod
    def dispersed_cube(dm, periods):
        nt = len(periods)
        nc, nb = 32, 128
        freq = np.tile(np.linspace(552.4, 1079.4, nc), (nt, 1))
        intrinsic = np.exp(-0.5 * ((np.arange(nb) - 63) / 2.2) ** 2)
        delay = 4148.808 * dm * (freq**-2 - np.max(freq, axis=1, keepdims=True)**-2)
        shifts = np.rint(delay * nb / periods[:, None]).astype(int)
        # Emulate upstream PulsarX archivewriter.h reintroducing channel delays.
        raw = np.array([[np.roll(intrinsic, shifts[t, c])
                         for c in range(nc)] for t in range(nt)])
        return raw, freq, intrinsic

    def test_stationary_profile_recovers(self):
        periods = np.full(3, 1 / 241)
        raw, freq, pulse = self.dispersed_cube(4, periods)
        restored = restore_pulsarx_dispersion(raw, freq, periods, 4)
        np.testing.assert_allclose(restored, np.broadcast_to(pulse, restored.shape), atol=1e-12)
        self.assertLess(np.max(raw.mean(axis=(0, 1))), np.max(pulse) * 0.5)

    def test_evolving_folding_period_recovers(self):
        periods = np.array([1 / 200, 1 / 280, 1 / 340])
        raw, freq, pulse = self.dispersed_cube(12, periods)
        restored = restore_pulsarx_dispersion(raw, freq, periods, 12)
        np.testing.assert_allclose(restored, np.broadcast_to(pulse, restored.shape), atol=1e-12)

    def test_zero_dm_preserves_profile(self):
        periods = np.array([0.005, 0.004])
        raw, freq, _ = self.dispersed_cube(0, periods)
        np.testing.assert_array_equal(
            restore_pulsarx_dispersion(raw, freq, periods, 0), raw)

    def test_invalid_period_is_rejected(self):
        periods = np.full(2, 0.005)
        raw, freq, _ = self.dispersed_cube(4, periods)
        with self.assertRaisesRegex(ValueError, "Invalid SUBINT"):
            restore_pulsarx_dispersion(raw, freq, np.array([0, 0.005]), 4)


if __name__ == "__main__":
    unittest.main()
