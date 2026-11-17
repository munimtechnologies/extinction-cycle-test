# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Offline tests for cycle.py (no network). Run: python -m unittest -v test_cycle"""
import csv
import io
import math
import unittest

import numpy as np

import cycle


class TestStages(unittest.TestCase):
    def test_contiguous(self):
        st = cycle.STAGES_GTS2020
        self.assertEqual(len(st), 96)
        for (n1, b1, t1), (n2, b2, t2) in zip(st, st[1:]):
            self.assertGreater(b1, t1, n1)
            self.assertAlmostEqual(t1, b2, places=6, msg=n1 + "/" + n2)
        self.assertEqual(st[-1], ("Quaternary", 2.58, 0.0))
        self.assertIn("Pridoli", [s[0] for s in st])

    def test_big_five_names_exist(self):
        names = {s[0] for s in cycle.STAGES_GTS2020}
        for b in cycle.BIG_FIVE_BINS:
            self.assertIn(b, names)


class TestAssign(unittest.TestCase):
    def setUp(self):
        self.a = cycle.BinAssigner()
        self.names = [s[0] for s in cycle.STAGES_GTS2020]

    def name(self, mx, mn):
        i = self.a(mx, mn)
        return None if i is None else self.names[i]

    def test_exact_stage(self):
        self.assertEqual(self.name(72.2, 66.0), "Maastrichtian")
        self.assertEqual(self.name(66.0, 61.66), "Danian")
        self.assertEqual(self.name(538.8, 529.0), "Fortunian")
        self.assertEqual(self.name(1.8, 0.774), "Quaternary")
        self.assertEqual(self.name(3.6, 2.58), "Piacenzian")

    def test_reversed_and_tolerance(self):
        self.assertEqual(self.name(66.0, 72.2), "Maastrichtian")
        self.assertEqual(self.name(72.25, 66.0), "Maastrichtian")

    def test_spanning(self):
        self.assertIsNone(self.name(72.2, 61.66))
        self.assertIsNone(self.name(600.0, 540.0))
        self.assertIsNone(self.name(None, 3.0))


class TestRates(unittest.TestCase):
    def setUp(self):
        self.sets = [{"A", "B", "C"}, {"A", "B"}, {"A", "C"}, {"A"}, {"A"}]

    def test_three_timer(self):
        mu, ps, c = cycle.three_timer(self.sets)
        self.assertAlmostEqual(ps[1], 0.5)
        self.assertAlmostEqual(ps[2], 1.0)
        self.assertAlmostEqual(ps[3], 1.0)
        self.assertTrue(math.isnan(ps[0]) and math.isnan(ps[4]))
        self.assertAlmostEqual(mu[1], math.log(2.0))
        self.assertAlmostEqual(mu[2], 0.0)
        self.assertTrue(math.isnan(mu[3]))
        self.assertEqual(c["two_timers"][1], 2)
        self.assertEqual(c["three_timers"][1], 1)
        self.assertEqual(c["part_timers"][1], 1)

    def test_foote(self):
        q, prop, c = cycle.foote_rates(self.sets)
        self.assertEqual(c["Nbt"][1], 2)
        self.assertEqual(c["NbL"][1], 1)
        self.assertAlmostEqual(q[1], -math.log(2 / 3))
        self.assertAlmostEqual(q[2], math.log(2.0))
        self.assertAlmostEqual(prop[1], 1 / 3)
        self.assertTrue(math.isnan(q[0]))

    def test_accumulator_from_csv(self):
        text = ("occurrence_no,collection_no,accepted_name,accepted_rank,max_ma,min_ma,"
                "phylum,genus\n"
                "1,10,Spirifer striatus,species,346.7,330.3,Brachiopoda,Spirifer\n"
                "2,11,Spirifer,genus,346.7,330.3,Brachiopoda,\n"
                "3,12,Ammonites,genus,200,150,Mollusca,Ammonites\n"
                "4,13,Foo bar,species,72.2,66,Mollusca,NO_GENUS_SPECIFIED\n"
                "5,14,,genus,72.2,66,Mollusca,\n")
        acc = cycle.Accumulator(96)
        acc.add_rows(csv.DictReader(io.StringIO(text)), cycle.BinAssigner())
        names = [s[0] for s in cycle.STAGES_GTS2020]
        vis = names.index("Visean")
        maa = names.index("Maastrichtian")
        self.assertEqual(acc.sets[vis], {"Brachiopoda|Spirifer"})
        self.assertEqual(acc.n_occ[vis], 2)
        self.assertEqual(acc.n_coll[vis], 2)
        self.assertEqual(acc.sets[maa], {"Mollusca|Foo"})
        self.assertEqual(acc.stats, {"rows": 5, "assigned": 3, "no_genus": 1, "no_bin": 1})


class TestStats(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(1)
        self.t = np.cumsum(rng.uniform(2, 8, 95))
        self.t = self.t[self.t < 470]
        self.rng = rng

    def test_circ_dist(self):
        P, t0 = 27.16, 9.5
        self.assertAlmostEqual(float(cycle.circ_dist(t0, P, t0)), 0.0)
        self.assertAlmostEqual(float(cycle.circ_dist(t0 + 3 * P + 2, P, t0)), 2.0)
        self.assertAlmostEqual(float(cycle.circ_dist(t0 - 1, P, t0)), 1.0)
        self.assertAlmostEqual(float(cycle.circ_dist(t0 + P / 2, P, t0)), P / 2)

    def test_detrend(self):
        y = 2 + 0.1 * self.t
        self.assertLess(np.max(np.abs(cycle.detrend(self.t, y))), 1e-9)
        Y = np.vstack([y, 3 * y - 1])
        self.assertLess(np.max(np.abs(cycle.detrend_rows(self.t, Y))), 1e-8)

    def test_matched_filter_detects_injected_cycle(self):
        y = cycle.template(self.t, 27.16, 9.5) + 0.5 * self.rng.standard_normal(len(self.t))
        y = y + 0.002 * self.t
        res = cycle.matched_filter_test(self.t, y, 27.16, 9.5, 2000, self.rng)
        self.assertGreater(res["r"], 0.5)
        self.assertLess(res["p_ar1"], 0.01)
        self.assertLess(res["p_phase"], 0.2)
        self.assertLess(res["phase_offset_myr"], 3.0)
        anti = cycle.matched_filter_test(self.t, -y, 27.16, 9.5, 500, self.rng)
        self.assertLess(anti["r"], -0.5)
        self.assertGreater(anti["p_ar1"], 0.9)

    def test_noise_not_significant_on_average(self):
        ps = []
        for _ in range(20):
            y = self.rng.standard_normal(len(self.t))
            ps.append(cycle.matched_filter_test(self.t, y, 27.16, 9.5, 300, self.rng)["p_ar1"])
        self.assertGreater(np.mean(ps), 0.3)

    def test_fit_phase(self):
        y = 2.0 * cycle.template(self.t, 27.0, 11.8)
        amp, tmax = cycle.fit_phase(self.t, y, 27.0)
        self.assertAlmostEqual(amp, 2.0, places=6)
        self.assertAlmostEqual(tmax, 11.8, places=6)

    def test_periodogram(self):
        y = np.sin(2 * np.pi * self.t / 27.0) + 0.3 * self.rng.standard_normal(len(self.t))
        res = cycle.periodogram_test(self.t, y, 1000, self.rng)
        self.assertTrue(26.0 <= res["global_peak_period_myr"] <= 28.5)
        self.assertTrue(res["global_peak_in_band"])
        self.assertLess(res["p_band_ar1"], 0.01)
        self.assertEqual(len(res["freqs"]), len(res["power"]))

    def test_ls_matches_direct_formula(self):
        y = self.rng.standard_normal(len(self.t))
        f = 0.037
        w = 2 * np.pi * f
        tau = math.atan2(np.sum(np.sin(2 * w * self.t)), np.sum(np.cos(2 * w * self.t))) / (2 * w)
        yc = y - y.mean()
        c = np.cos(w * (self.t - tau))
        s = np.sin(w * (self.t - tau))
        direct = ((yc @ c) ** 2 / (c @ c) + (yc @ s) ** 2 / (s @ s)) / (2 * y.var(ddof=1))
        got = cycle.ls_power(y, cycle.ls_matrices(self.t, np.array([f])))[0, 0]
        self.assertAlmostEqual(direct, got, places=9)

    def test_binom_and_events(self):
        self.assertEqual(cycle.binom_sf(0, 5, 0.3), 1.0)
        self.assertAlmostEqual(cycle.binom_sf(5, 5, 0.3), 0.3 ** 5)
        self.assertAlmostEqual(cycle.binom_sf(1, 3, 0.5), 0.875)
        y = cycle.template(self.t, 27.16, 9.5)
        ev = cycle.event_test(self.t, y, 27.16, 9.5)
        self.assertGreater(ev["n_peaks"], 5)
        self.assertGreaterEqual(ev["hits"], 0.6 * ev["n_peaks"])

    def test_peaks(self):
        # 0.75-quantile of these nine values is 2, so the peak of height 1 is dropped
        self.assertEqual(cycle.peaks([0, 5, 0, 1, 0, 9, 0, 2, 0], q=0.75), [1, 5, 7])

    def test_normal_scores(self):
        z = cycle.normal_scores([3.0, 1.0, 2.0, 10.0])
        self.assertTrue(z[1] < z[2] < z[0] < z[3])
        self.assertAlmostEqual(float(np.sum(z)), 0.0, places=9)

    def test_redated_hits(self):
        rows = cycle.redated_hits()
        self.assertEqual(len(rows), 19)
        self.assertEqual(sum(r["hit_2012"] for r in rows), 10)

    def test_jitter(self):
        y = cycle.template(self.t, 27.16, 9.5)
        j = cycle.jitter_r(self.t, y, 27.16, 9.5, 50, self.rng)
        self.assertGreater(j["r_p05"], 0.8)


if __name__ == "__main__":
    unittest.main()
