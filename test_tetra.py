# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Offline tests for tetra.py (PREREG-3). Run: python -m unittest -v test_tetra"""
import math
import unittest

import numpy as np

import cycle
import tetra

ST = cycle.STAGES_GTS2020
NAMES = [s[0] for s in ST]
TOPS = np.array([s[2] for s in ST], float)


class EventsAndPhase(unittest.TestCase):
    def test_event_bins(self):
        got = {ST[i][0] for i in tetra.event_bins(ST)}
        self.assertEqual(got, {"Tortonian", "Priabonian", "Maastrichtian", "Tithonian",
                               "Rhaetian", "Norian", "Carnian", "Changhsingian",
                               "Capitanian"})

    def test_event_bin_prefers_nearest_top(self):
        # 251.9 Ma is 0.002 from the Changhsingian top and 1.1 from the Induan top.
        self.assertEqual([ST[i][0] for i in tetra.event_bins(ST, [251.9])], ["Changhsingian"])
        # 215 Ma is far from any top: the containing stage is used.
        self.assertEqual([ST[i][0] for i in tetra.event_bins(ST, [215.0])], ["Norian"])

    def test_part_a_offsets(self):
        d = cycle.circ_dist([a for _, a in tetra.TETRAPOD_EVENTS], 27.5, 7.0)
        self.assertEqual(int(np.sum(d <= 3.0)), 6)
        self.assertAlmostEqual(float(d[2]), 4.0, places=6)    # K-Pg, 66 Ma
        self.assertAlmostEqual(float(d[5]), 12.0, places=6)   # mid-Norian, 215 Ma

    def test_circular_phase(self):
        t0, r = tetra.circular_phase([7.0, 34.5], 27.5)
        self.assertAlmostEqual(t0, 7.0, places=9)
        self.assertAlmostEqual(r, 1.0, places=9)
        t0, r = tetra.circular_phase(tetra.PULSES_RCZ2021, 27.5)
        self.assertTrue(6.5 < t0 < 8.5, t0)
        self.assertTrue(0.25 < r < 0.45, r)

    def test_queries(self):
        # Amendment 1: PBDB rejects envtype=terr; 'terrestrial' = any terrestrial env.
        for url in tetra.QUERIES.values():
            self.assertIn("envtype=terrestrial&", url)
            self.assertNotIn("envtype=terr&", url)
            self.assertIn("idreso=genus", url)
            self.assertIn("show=class", url)
        self.assertIn("Tetrapoda%5EMammalia", tetra.QUERIES["tetra_nonmammal"])


class Analyze(unittest.TestCase):
    def synthetic(self, seed=1, amp=1.0, noise=0.2):
        rng = np.random.default_rng(seed)
        mu = amp * cycle.template(TOPS, 27.5, 7.0) + noise * rng.standard_normal(len(ST))
        mu[NAMES.index("Maastrichtian")] = 2.0
        return mu

    def test_window_has_57_bins(self):
        res, _ = tetra.analyze(ST, self.synthetic(), 200, np.random.default_rng(3),
                               full=False)
        self.assertEqual(res["n_bins"], 57)
        self.assertEqual(res["window_bins"][0], "Gzhelian")
        self.assertEqual(res["window_bins"][-1], "Piacenzian")

    def test_planted_cycle_is_found(self):
        res, (f, p) = tetra.analyze(ST, self.synthetic(), 400, np.random.default_rng(7))
        self.assertTrue(res["sanity"]["passed"])
        self.assertGreater(res["H_rcz"]["r"], 0.7)
        self.assertLess(res["H_rcz"]["p_ar1"], 0.05)
        self.assertLess(res["H_rcz"]["phase_offset_myr"], 2.0)
        self.assertGreater(res["leave_events_out"]["r"], 0.7)
        self.assertEqual(len(res["leave_events_out"]["dropped_bins"]), 9)
        self.assertEqual(res["leave_events_out"]["n"], 48)
        self.assertEqual(res["decision"], "SUPPORTED_BEYOND_EVENTS")
        self.assertEqual(len(f), len(p))
        self.assertTrue(25.0 < res["periodogram"]["global_peak_period_myr"] < 30.0)

    def test_antiphase_is_not_supported(self):
        rng = np.random.default_rng(5)
        mu = -cycle.template(TOPS, 27.5, 7.0) + 0.2 * rng.standard_normal(len(ST))
        mu[NAMES.index("Maastrichtian")] = 2.0
        res, _ = tetra.analyze(ST, mu, 300, np.random.default_rng(11), full=False)
        self.assertLess(res["H_rcz"]["r"], 0.0)
        self.assertEqual(res["decision"], "NOT_SUPPORTED")

    def test_sanity_gate(self):
        mu = self.synthetic()
        mu[NAMES.index("Maastrichtian")] = -5.0
        res, _ = tetra.analyze(ST, mu, 200, np.random.default_rng(2), full=False)
        self.assertFalse(res["sanity"]["checks"]["maastrichtian_top_quartile"])
        self.assertEqual(res["decision"], "DATA_LIMITED")

    def test_empty_series(self):
        res, (f, p) = tetra.analyze(ST, [math.nan] * len(ST), 100, np.random.default_rng(1))
        self.assertEqual(res["n_bins"], 0)
        self.assertEqual(res["decision"], "DATA_LIMITED")
        self.assertIsNone(f)

    def test_noise_runs(self):
        rng = np.random.default_rng(9)
        mu = rng.standard_normal(len(ST))
        res, _ = tetra.analyze(ST, mu, 200, np.random.default_rng(4))
        for key in ("H_rcz", "leave_events_out", "H_mb2014", "H_rcz_ranks",
                    "duration_artifact"):
            self.assertTrue(0.0 < res[key]["p_ar1"] <= 1.0)
            self.assertTrue(-1.0 <= res[key]["r"] <= 1.0)
        self.assertIn(res["decision"], ("DATA_LIMITED", "NOT_SUPPORTED",
                                        "SUPPORTED_ONLY_BY_EVENTS",
                                        "SUPPORTED_BEYOND_EVENTS"))


if __name__ == "__main__":
    unittest.main()
