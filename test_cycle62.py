# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Offline tests for cycle62.py (PREREG-2). Made with AI agents working in MV2's name."""
import unittest

import numpy as np

import cycle
import cycle62 as c6


def window_mids():
    mids = np.array([(b + t) / 2.0 for _, b, t in cycle.STAGES_GTS2020])
    lo, hi = c6.WINDOW62
    return np.sort(mids[(mids >= lo) & (mids <= hi)])


class TestRangeThrough(unittest.TestCase):
    def test_basic(self):
        sets = [{"a"}, {"a", "b"}, {"b"}, {"c"}]
        self.assertEqual(c6.range_through(sets), [1, 2, 1, 1])

    def test_gap_filled(self):
        self.assertEqual(c6.range_through([{"a"}, set(), {"a"}]), [1, 1, 1])

    def test_short_lived(self):
        mids = [100.0, 80.0, 40.0, 10.0]
        sets = [{"a", "c"}, {"a"}, {"b"}, {"b", "c"}]
        self.assertEqual(c6.range_through(sets), [2, 2, 2, 2])
        self.assertEqual(c6.range_through(sets, mids, 45.0), [1, 1, 1, 1])
        self.assertEqual(c6.range_through(sets, mids, 10.0), [0, 0, 0, 0])


class TestDetrend(unittest.TestCase):
    def test_cubic_removed(self):
        t = window_mids()
        y = 2.0 + 0.1 * t - 1e-3 * t ** 2 + 1e-6 * t ** 3
        self.assertLess(float(np.max(np.abs(c6.detrend_poly(t, y, 3)))), 1e-6)

    def test_rows_match_single(self):
        t = window_mids()
        rng = np.random.default_rng(1)
        Y = rng.standard_normal((3, len(t)))
        R = c6.detrend_rows_poly(t, Y, 3)
        for k in range(3):
            np.testing.assert_allclose(R[k], c6.detrend_poly(t, Y[k], 3), atol=1e-8)


class TestMatchedFilter(unittest.TestCase):
    def setUp(self):
        self.t = window_mids()
        rng = np.random.default_rng(7)
        self.trend = 500.0 + 2.0 * self.t - 0.002 * self.t ** 2
        self.noise = 0.3 * rng.standard_normal(len(self.t))

    def test_in_phase_signal_detected(self):
        y = self.trend + 3.0 * cycle.template(self.t, 62.0, 27.8) + self.noise
        res = c6.matched_filter_poly(self.t, y, 62.0, 27.8, 2000, np.random.default_rng(3))
        self.assertGreater(res["r"], 0.8)
        self.assertLess(res["p_ar1"], 0.05)
        self.assertLess(res["p_phase"], 0.1)
        self.assertLess(res["phase_offset_myr"], 3.0)

    def test_antiphase_signal_rejected(self):
        y = self.trend + 3.0 * cycle.template(self.t, 62.0, 27.8 + 31.0) + self.noise
        res = c6.matched_filter_poly(self.t, y, 62.0, 27.8, 2000, np.random.default_rng(3))
        self.assertLess(res["r"], -0.8)
        self.assertGreater(res["p_ar1"], 0.5)

    def test_periodogram_finds_62(self):
        y = self.trend + 3.0 * cycle.template(self.t, 62.0, 27.8) + self.noise
        pg = c6.periodogram_poly(self.t, y, 500, np.random.default_rng(5))
        self.assertTrue(58.0 <= pg["global_peak_period_myr"] <= 66.0)
        self.assertTrue(pg["global_peak_in_band"])


class TestDeclining(unittest.TestCase):
    def test_known_ages(self):
        got = [bool(x) for x in c6.declining(np.array([20.0, 40.0, 380.0, 399.0, 368.0, 27.8]))]
        self.assertEqual(got, [True, False, True, True, False, True])

    def test_mb_events(self):
        rows = c6.mb_events_declining()
        self.assertEqual(len(rows), 19)
        self.assertEqual(sum(r["declining_2020"] for r in rows), 16)
        self.assertEqual(sum(r["declining_2012"] for r in rows), 15)

    def test_big_five_tops_declining(self):
        tops = {n: t for n, _, t in cycle.STAGES_GTS2020}
        self.assertTrue(all(bool(c6.declining(tops[b])) for b in cycle.BIG_FIVE_BINS))

    def test_phase_contrast_detects(self):
        t = np.array(sorted(tp for _, _, tp in cycle.STAGES_GTS2020
                            if cycle.WINDOW[0] <= tp <= cycle.WINDOW[1]))
        rng = np.random.default_rng(11)
        y = c6.declining(t).astype(float) + 0.01 * rng.standard_normal(len(t))
        res = c6.phase_contrast(t, y, 2000, np.random.default_rng(2))
        self.assertGreater(res["diff_normal_scores"], 1.0)
        self.assertLess(res["p_random_phase"], 0.05)
        self.assertTrue(0 < res["n_declining"] < res["n"])


if __name__ == "__main__":
    unittest.main()
