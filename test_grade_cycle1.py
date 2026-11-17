# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Offline tests for grade_cycle1.py (standard library only).
Made with AI agents working in MV2's name."""
import json
import os
import tempfile
import unittest

import grade_cycle1 as gc

NAN = float("nan")


def _blk(r=0.1, p=0.2, off=5.0, p_binom=0.3):
    return {"MB2014": {"matched_filter": {"n": 85, "r": r, "p_ar1": p, "p_phase": 0.3,
                                          "fit_tmax_ma": 9.5 + off, "phase_offset_myr": off},
                       "events": {"n_peaks": 12, "hits": 3, "p0": 0.22, "p_binom": p_binom}}}


def make(r=0.1, p=0.2, passed=True, p_wo=0.4, p_ns=0.3, in_band=False, p_band=0.5,
         off=5.0, p_binom=0.3, decision=None, tops=None):
    mu = _blk(r, p, off, p_binom)
    mu["periodogram"] = {"band": [26.0, 28.5], "band_max_power": 3.0, "p_band_ar1": p_band,
                         "global_peak_period_myr": 62.0 if not in_band else 27.3,
                         "global_peak_in_band": in_band, "top_peaks": []}
    mu["jitter_MB2014"] = {"r_p05": 0.0, "r_p50": 0.1, "r_p95": 0.2, "frac_r_pos": 0.9}
    tops = dict(gc.BIG_FIVE_TOPS_MA) if tops is None else tops
    res = {"status": "ok",
           "sanity": {"top10_mu3t_bins": ["Changhsingian"], "big_five_in_top10": 5,
                      "passed": passed, "n_bins_in_window": 85},
           "mu3t": mu,
           "mu3t_bins": [{"bin": k, "top_ma": v, "mu3t": 1.0} for k, v in tops.items()],
           "mu3t_normal_scores": _blk(0.05, p_ns),
           "mu3t_without_big_five": _blk(0.02, p_wo)}
    if decision is None:
        decision = ("SANITY_FAILED" if not passed
                    else ("SUPPORTED" if (r > 0 and p < 0.05) else "NOT_SUPPORTED"))
    res["decision_H1"] = decision
    return res


class TestEphemeris(unittest.TestCase):
    def test_offsets_of_big_five(self):
        expect = {"Maastrichtian": 2.18, "Rhaetian": 1.78, "Changhsingian": 2.038,
                  "Frasnian": 9.57, "Hirnantian": 0.96, "Katian": 1.14}
        for k, v in expect.items():
            self.assertAlmostEqual(gc.offset(gc.BIG_FIVE_TOPS_MA[k]), v, places=6, msg=k)

    def test_template_signs(self):
        self.assertGreater(gc.template(66.0), 0.85)
        self.assertLess(gc.template(372.15), -0.55)
        self.assertAlmostEqual(gc.template(9.5), 1.0, places=12)


class TestGrade(unittest.TestCase):
    def rows(self, res):
        return {r["id"]: r["holds"] for r in gc.grade(res)["rows"]}

    def test_default(self):
        gr = gc.grade(make())
        self.assertEqual(gr["status"], "GRADED")
        self.assertEqual(gr["decision"], "NOT_SUPPORTED")
        self.assertEqual(gr["headline"], "NOT_SUPPORTED")
        r = self.rows(make())
        self.assertTrue(r["P1"])
        self.assertTrue(r["P2"])
        self.assertFalse(r["P3"])
        self.assertTrue(r["P4"])
        self.assertTrue(r["P5"])
        self.assertTrue(r["P6"])
        self.assertTrue(r["P7"])
        self.assertFalse(r["P8"])     # offset 5 > 3
        self.assertTrue(r["P9"])

    def test_supported_beyond(self):
        gr = gc.grade(make(r=0.3, p=0.01, p_wo=0.02))
        self.assertEqual(gr["decision"], "SUPPORTED")
        self.assertEqual(gr["headline"], "PREDICTS_BEYOND_MASS_EXTINCTIONS")

    def test_supported_rests_on_big_five(self):
        gr = gc.grade(make(r=0.3, p=0.01, p_wo=0.3))
        self.assertEqual(gr["headline"], "RESTS_ON_MASS_EXTINCTIONS")

    def test_negative_r_small_p_not_supported(self):
        gr = gc.grade(make(r=-0.3, p=0.01))
        self.assertEqual(gr["decision"], "NOT_SUPPORTED")
        self.assertFalse(self.rows(make(r=-0.3, p=0.01))["P2"])

    def test_nan_p_rows_false(self):
        r = self.rows(make(p_wo=NAN, p_ns=NAN, p_band=NAN, p_binom=NAN))
        for k in ("P4", "P5", "P7", "P9"):
            self.assertFalse(r[k], k)

    def test_in_band_and_phase(self):
        r = self.rows(make(in_band=True, off=2.5))
        self.assertFalse(r["P6"])
        self.assertTrue(r["P8"])

    def test_sanity_failed(self):
        gr = gc.grade(make(passed=False, r=0.3, p=0.01))
        self.assertEqual(gr["decision"], "SANITY_FAILED")
        self.assertEqual(gr["status"], "GRADED_SANITY_FAILED")
        self.assertEqual(gr["brier"]["n"], 1)

    def test_check_code(self):
        gr = gc.grade(make(decision="SUPPORTED"))
        self.assertEqual(gr["status"], "CHECK_CODE")

    def test_check_data(self):
        tops = dict(gc.BIG_FIVE_TOPS_MA)
        tops["Frasnian"] = 371.1
        del tops["Rhaetian"]
        gr = gc.grade(make(tops=tops))
        self.assertEqual(gr["status"], "CHECK_DATA")
        self.assertEqual(gr["tops_check"]["missing_bins"], ["Rhaetian"])
        self.assertEqual([m["bin"] for m in gr["tops_check"]["mismatches"]], ["Frasnian"])

    def test_job_failed(self):
        gr = gc.grade({"status": "data_incomplete: missing Mollusca"})
        self.assertEqual(gr["status"], "JOB_FAILED")
        self.assertIn("Nothing is graded", gc.render_md(gr, {}))

    def test_main(self):
        res = make(r=0.3, p=0.01, p_wo=0.3)
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "results.json")
            with open(p, "w", encoding="utf-8") as f:
                f.write(json.dumps(res))
            self.assertEqual(gc.main(["--results", p, "--out-dir", d]), 0)
            with open(os.path.join(d, "graded1.json"), encoding="utf-8") as f:
                gr = json.load(f)
            self.assertEqual(gr["headline"], "RESTS_ON_MASS_EXTINCTIONS")
            with open(os.path.join(d, "GRADED-1.md"), encoding="utf-8") as f:
                md = f.read()
            self.assertIn("| P3 |", md)
            self.assertIn("RESTS_ON_MASS_EXTINCTIONS", md)


if __name__ == "__main__":
    unittest.main()
