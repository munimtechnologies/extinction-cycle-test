# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Offline tests for grade_cycle2.py (standard library only).
Made with AI agents working in MV2's name."""
import json
import os
import tempfile
import unittest

import grade_cycle2 as gc

NAN = float("nan")


def _mf(r=0.1, p=0.2, off=20.0, n=90):
    return {"n": n, "r": r, "lag1": 0.8, "p_ar1": p, "p_phase": 0.3,
            "fit_amplitude": 10.0, "fit_tmax_ma": 27.8 + off, "phase_offset_myr": off,
            "detrend_deg": 3}


def make(r=0.1, p=0.2, off=20.0, passed=True, in_band=False, p_band=0.5,
         r_short=0.1, p_short=0.3, r_coll=0.1, p_coll=0.3, p_h7=0.3,
         decision=None, n=90, n_window=90):
    res = {"prereg": "PREREG-2.md (extinction-cycle-test)", "status": "ok",
           "sanity": {"checks": {"induan_lt_changhsingian": passed,
                                 "danian_lt_maastrichtian": True,
                                 "rhuddanian_lt_katian": True},
                      "passed": passed},
           "window_bins": ["bin%d" % i for i in range(n_window)],
           "H6_range_through": _mf(r, p, off, n),
           "H6_periodogram": {"band": [58.0, 66.0], "band_max_power": 3.0,
                              "p_band_ar1": p_band,
                              "global_peak_period_myr": 61.0 if in_band else 140.0,
                              "global_peak_in_band": in_band, "top_peaks": []},
           "H6_jitter": {"r_p05": 0.0, "r_p50": 0.1, "r_p95": 0.2, "frac_r_pos": 0.8},
           "H6_quadratic_detrend": _mf(0.05, 0.4, off, n),
           "short_lived": _mf(r_short, p_short, off, n),
           "sampled_in_bin": _mf(0.0, 0.5, off, n),
           "log_collections": _mf(r_coll, p_coll, off, n),
           "short_lived_periodogram": {"band": [58.0, 66.0], "p_band_ar1": 0.5,
                                       "global_peak_period_myr": 100.0,
                                       "global_peak_in_band": False, "top_peaks": []},
           "H7_mu3t_declining": {"n": 83, "n_declining": 40, "diff_normal_scores": 0.1,
                                 "p_random_phase": p_h7}}
    if decision is None:
        decision = ("SANITY_FAILED" if not passed
                    else ("SUPPORTED" if (r > 0 and p < 0.05) else "NOT_SUPPORTED"))
    res["decision_H6"] = decision
    return res


class TestGrade(unittest.TestCase):
    def rows(self, res):
        return {r["id"]: r["holds"] for r in gc.grade(res)["rows"]}

    def test_default(self):
        gr = gc.grade(make())
        self.assertEqual(gr["status"], "GRADED")
        self.assertEqual(gr["decision"], "NOT_SUPPORTED")
        self.assertEqual(gr["headline"], "NOT_SUPPORTED")
        self.assertEqual(gr["failure_mode"], "PERIOD_FAILED")
        r = self.rows(make())
        self.assertTrue(r["Q1"])
        self.assertTrue(r["Q2"])
        self.assertFalse(r["Q3"])
        self.assertFalse(r["Q4"])     # offset 20 > 8
        self.assertFalse(r["Q5"])
        self.assertFalse(r["Q6"])
        self.assertFalse(r["Q7"])
        self.assertFalse(r["Q8"])
        self.assertFalse(r["Q9"])
        self.assertEqual(gr["brier"]["n"], 9)

    def test_supported(self):
        gr = gc.grade(make(r=0.4, p=0.01))
        self.assertEqual(gr["decision"], "SUPPORTED")
        self.assertEqual(gr["headline"], "PREDICTS_DIVERSITY")
        self.assertIsNone(gr["failure_mode"])

    def test_supported_with_sampling(self):
        gr = gc.grade(make(r=0.4, p=0.01, r_coll=0.3, p_coll=0.01))
        self.assertEqual(gr["headline"], "PREDICTS_DIVERSITY_AND_SAMPLING")

    def test_phase_failed(self):
        gr = gc.grade(make(in_band=True, off=20.0))
        self.assertEqual(gr["failure_mode"], "PHASE_FAILED")

    def test_weak(self):
        gr = gc.grade(make(p_band=0.01, off=3.0))
        self.assertEqual(gr["failure_mode"], "WEAK")
        r = self.rows(make(p_band=0.01, off=3.0))
        self.assertTrue(r["Q4"])
        self.assertTrue(r["Q6"])

    def test_negative_r_small_p(self):
        gr = gc.grade(make(r=-0.4, p=0.01))
        self.assertEqual(gr["decision"], "NOT_SUPPORTED")
        self.assertFalse(self.rows(make(r=-0.4, p=0.01))["Q2"])

    def test_secondary_rows(self):
        r = self.rows(make(r_short=0.3, p_short=0.01, r_coll=0.3, p_coll=0.02, p_h7=0.01))
        self.assertTrue(r["Q7"])
        self.assertTrue(r["Q8"])
        self.assertTrue(r["Q9"])

    def test_nan_rows_false(self):
        r = self.rows(make(p=NAN, p_band=NAN, p_short=NAN, p_coll=NAN, p_h7=NAN, off=NAN))
        for k in ("Q3", "Q4", "Q6", "Q7", "Q8", "Q9"):
            self.assertFalse(r[k], k)

    def test_sanity_failed(self):
        gr = gc.grade(make(passed=False, r=0.4, p=0.01))
        self.assertEqual(gr["decision"], "SANITY_FAILED")
        self.assertEqual(gr["headline"], "SANITY_FAILED")
        self.assertEqual(gr["status"], "GRADED_SANITY_FAILED")
        self.assertEqual(gr["brier"]["n"], 1)
        self.assertIsNone(gr["failure_mode"])

    def test_check_code(self):
        gr = gc.grade(make(decision="SUPPORTED"))
        self.assertEqual(gr["status"], "CHECK_CODE")

    def test_check_data(self):
        gr = gc.grade(make(n=90, n_window=89))
        self.assertEqual(gr["status"], "CHECK_DATA")
        self.assertFalse(gr["window_check"]["passed"])

    def test_job_failed(self):
        gr = gc.grade({"status": "data_incomplete"})
        self.assertEqual(gr["status"], "JOB_FAILED")
        self.assertIn("Nothing is graded", gc.render_md(gr, {}))

    def test_brier_value(self):
        # default: Q1 T, Q2 T, Q3-Q9 F
        gr = gc.grade(make())
        expect = ((0.15 ** 2) + (0.30 ** 2) + (0.40 ** 2) + (0.45 ** 2) + (0.30 ** 2)
                  + (0.30 ** 2) + (0.40 ** 2) + (0.30 ** 2) + (0.40 ** 2)) / 9.0
        self.assertAlmostEqual(gr["brier"]["brier"], expect, places=9)

    def test_main(self):
        res = make(r=0.4, p=0.01)
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "results62.json")
            with open(p, "w", encoding="utf-8") as f:
                f.write(json.dumps(res))
            self.assertEqual(gc.main(["--results", p, "--out-dir", d]), 0)
            with open(os.path.join(d, "graded2.json"), encoding="utf-8") as f:
                gr = json.load(f)
            self.assertEqual(gr["headline"], "PREDICTS_DIVERSITY")
            with open(os.path.join(d, "GRADED-2.md"), encoding="utf-8") as f:
                md = f.read()
            self.assertIn("| Q3 |", md)
            self.assertIn("PREDICTS_DIVERSITY", md)


if __name__ == "__main__":
    unittest.main()
