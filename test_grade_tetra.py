# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Offline tests for grade_tetra.py (standard library only).
Made with AI agents working in MV2's name."""
import json
import os
import tempfile
import unittest

import grade_tetra as gt

NAN = float("nan")
EV = list(gt.EVENT_STAGES)


def _mf(r=0.1, p=0.2, pph=0.3, n=50, off=10.0):
    return {"n": n, "r": r, "lag1": 0.4, "p_ar1": p, "p_phase": pph,
            "fit_amplitude": 0.1, "fit_tmax_ma": 7.0 + off, "phase_offset_myr": off}


def make(r=0.1, p=0.2, pph=0.3, r2=0.05, p2=0.4, passed=True, in_band=False,
         p_band=0.5, r_mb=0.1, p_dur=0.5, r_nm=0.1, p_binom=0.5, decision=None,
         n_window=50, n_h1=None, n_h2=None, dropped=None, nonmammal=True):
    names = EV + ["s%d" % i for i in range(n_window - len(EV))]
    dropped = EV if dropped is None else dropped
    n_h1 = n_window if n_h1 is None else n_h1
    n_h2 = (n_window - len(dropped)) if n_h2 is None else n_h2
    lo = _mf(r2, p2, 0.4, n_h2)
    lo["dropped_bins"] = list(dropped)
    prim = {"window": [2.5, 300.0], "n_bins": n_window, "window_bins": names,
            "sanity": {"checks": {"maastrichtian_top_quartile": passed,
                                  "enough_bins": True}, "passed": passed},
            "H_rcz": _mf(r, p, pph, n_h1),
            "leave_events_out": lo,
            "H_mb2014": _mf(r_mb, 0.3, 0.3, n_h1),
            "H_rcz_ranks": _mf(r, 0.3, 0.3, n_h1),
            "H_rcz_jitter": {"r_p05": 0.0, "r_p50": 0.1, "r_p95": 0.2, "frac_r_pos": 0.8},
            "periodogram": {"band": [26.0, 29.0], "band_max_power": 2.0,
                            "p_band_ar1": p_band,
                            "global_peak_period_myr": 27.4 if in_band else 41.0,
                            "global_peak_in_band": in_band, "top_peaks": []},
            "event_test": {"n_peaks": 9, "hits": 2, "p0": 6 / 27.5, "p_binom": p_binom,
                           "peak_ages_ma": [], "peak_offsets_myr": []},
            "duration_artifact": _mf(0.0, p_dur, 0.5, n_h1)}
    if decision is None:
        if not passed:
            decision = "DATA_LIMITED"
        elif r > 0 and p < 0.05:
            decision = ("SUPPORTED_BEYOND_EVENTS" if (r2 > 0 and p2 < 0.05)
                        else "SUPPORTED_ONLY_BY_EVENTS")
        else:
            decision = "NOT_SUPPORTED"
    prim["decision"] = decision
    res = {"prereg": "PREREG-3.md + PREREG-3-AMENDMENT-1.md (extinction-cycle-test)",
           "status": "ok", "primary": prim, "decision": decision,
           "window_260": {"decision": "NOT_SUPPORTED", "H_rcz": _mf(0.05, 0.4)},
           "foote_q": {"decision": "NOT_SUPPORTED", "H_rcz": _mf(0.05, 0.4)}}
    if nonmammal:
        res["nonmammal"] = {"decision": "NOT_SUPPORTED", "H_rcz": _mf(r_nm, 0.4)}
    return res


class TestGrade(unittest.TestCase):
    def rows(self, res):
        return {r["id"]: r["holds"] for r in gt.grade(res)["rows"]}

    def test_default(self):
        gr = gt.grade(make())
        self.assertEqual(gr["status"], "GRADED")
        self.assertEqual(gr["decision"], "NOT_SUPPORTED")
        self.assertEqual(gr["headline"], "NOT_SUPPORTED")
        r = self.rows(make())
        expect = {"T1": True, "T2": True, "T3": False, "T4": True, "T5": False,
                  "T6": True, "T7": True, "T8": True, "T9": True, "T10": True,
                  "T11": True}
        self.assertEqual(r, expect)
        self.assertEqual(gr["brier"]["n"], 11)

    def test_brier_value(self):
        # default: all rows TRUE except T3 and T5
        gr = gt.grade(make())
        sq = 0.0
        for k, p, _ in gt.PRIORS:
            hit = 0.0 if k in ("T3", "T5") else 1.0
            sq += (p - hit) ** 2
        self.assertAlmostEqual(gr["brier"]["brier"], sq / 11.0, places=9)

    def test_only_by_events(self):
        gr = gt.grade(make(r=0.3, p=0.01))
        self.assertEqual(gr["decision"], "SUPPORTED_ONLY_BY_EVENTS")
        self.assertEqual(gr["headline"], "RESTS_ON_FITTED_EVENTS")
        self.assertEqual(gr["status"], "GRADED")
        self.assertTrue(self.rows(make(r=0.3, p=0.01))["T3"])

    def test_beyond_events(self):
        gr = gt.grade(make(r=0.3, p=0.01, r2=0.25, p2=0.02))
        self.assertEqual(gr["decision"], "SUPPORTED_BEYOND_EVENTS")
        self.assertEqual(gr["headline"], "INDEPENDENT_SUPPORT")
        self.assertFalse(self.rows(make(r=0.3, p=0.01, r2=0.25, p2=0.02))["T4"])

    def test_beyond_events_duration_artifact(self):
        gr = gt.grade(make(r=0.3, p=0.01, r2=0.25, p2=0.02, p_dur=0.01))
        self.assertEqual(gr["decision"], "SUPPORTED_BEYOND_EVENTS")
        self.assertEqual(gr["headline"], "SUPPORT_BUT_DURATION_ARTIFACT")

    def test_negative_r_small_p(self):
        gr = gt.grade(make(r=-0.3, p=0.01))
        self.assertEqual(gr["decision"], "NOT_SUPPORTED")
        r = self.rows(make(r=-0.3, p=0.01))
        self.assertFalse(r["T2"])
        self.assertFalse(r["T3"])
        self.assertFalse(r["T10"])   # non-mammal r 0.1 has the other sign

    def test_periodogram_rows(self):
        r = self.rows(make(in_band=True, p_band=0.01))
        self.assertFalse(r["T6"])
        self.assertFalse(r["T7"])

    def test_phase_and_events(self):
        r = self.rows(make(pph=0.01, p_binom=0.01, r_mb=-0.2))
        self.assertTrue(r["T5"])
        self.assertFalse(r["T11"])
        self.assertFalse(r["T8"])

    def test_nan_never_satisfies(self):
        r = self.rows(make(p=NAN, pph=NAN, p2=NAN, p_band=NAN, p_dur=NAN, p_binom=NAN,
                           r_nm=NAN))
        for k in ("T3", "T4", "T5", "T7", "T9", "T10", "T11"):
            self.assertFalse(r[k], k)

    def test_nonmammal_missing(self):
        r = self.rows(make(nonmammal=False))
        self.assertIsNone(r["T10"])
        gr = gt.grade(make(nonmammal=False))
        self.assertEqual(gr["brier"]["n"], 10)

    def test_nonmammal_too_short(self):
        res = make()
        res["nonmammal"] = {"decision": "DATA_LIMITED", "n_bins": 5}
        self.assertFalse(self.rows(res)["T10"])

    def test_data_limited(self):
        gr = gt.grade(make(passed=False, r=0.3, p=0.01))
        self.assertEqual(gr["decision"], "DATA_LIMITED")
        self.assertEqual(gr["headline"], "DATA_LIMITED")
        self.assertEqual(gr["status"], "GRADED_DATA_LIMITED")
        self.assertEqual(gr["brier"]["n"], 1)
        self.assertIn("nothing is claimed", gt.render_md(gr, make(passed=False)))

    def test_data_limited_short_series(self):
        # tetra.analyze stops before any test when fewer than 12 bins exist
        res = {"status": "ok", "decision": "DATA_LIMITED",
               "primary": {"window": [2.5, 300.0], "n_bins": 8,
                           "window_bins": ["a", "b", "c", "d", "e", "f", "g", "h"],
                           "sanity": {"checks": {"maastrichtian_top_quartile": False,
                                                 "enough_bins": False}, "passed": False},
                           "decision": "DATA_LIMITED"}}
        gr = gt.grade(res)
        self.assertEqual(gr["status"], "GRADED_DATA_LIMITED")
        self.assertEqual(gr["decision"], "DATA_LIMITED")
        r = {x["id"]: x["holds"] for x in gr["rows"]}
        self.assertFalse(r["T1"])
        self.assertIsNone(r["T2"])
        self.assertTrue(gr["checks"]["data_ok"])
        self.assertIn("Graded result", gt.render_md(gr, res))

    def test_check_code_decision(self):
        gr = gt.grade(make(decision="SUPPORTED_ONLY_BY_EVENTS"))
        self.assertEqual(gr["status"], "CHECK_CODE")
        self.assertEqual(gr["decision"], "NOT_SUPPORTED")

    def test_check_code_events(self):
        gr = gt.grade(make(dropped=EV[:8]))
        self.assertFalse(gr["checks"]["events_ok"])
        self.assertEqual(gr["status"], "CHECK_CODE")

    def test_check_data(self):
        gr = gt.grade(make(n_h1=49))
        self.assertEqual(gr["status"], "CHECK_DATA")
        gr = gt.grade(make(n_h2=42))
        self.assertEqual(gr["status"], "CHECK_DATA")

    def test_job_failed(self):
        gr = gt.grade({"status": "data_incomplete"})
        self.assertEqual(gr["status"], "JOB_FAILED")
        self.assertIn("Nothing is graded", gt.render_md(gr, {}))

    def test_main(self):
        res = make(r=0.3, p=0.01)
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "results_tetra.json")
            with open(p, "w", encoding="utf-8") as f:
                f.write(json.dumps(res))
            self.assertEqual(gt.main(["--results", p, "--out-dir", d]), 0)
            with open(os.path.join(d, "graded3.json"), encoding="utf-8") as f:
                gr = json.load(f)
            self.assertEqual(gr["headline"], "RESTS_ON_FITTED_EVENTS")
            with open(os.path.join(d, "GRADED-3.md"), encoding="utf-8") as f:
                md = f.read()
            self.assertIn("| T4 |", md)
            self.assertIn("RESTS_ON_FITTED_EVENTS", md)


if __name__ == "__main__":
    unittest.main()
