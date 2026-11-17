# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""test_graded_results.py: runs the graders that were published before the results
(grade_cycle2.py, grade_tetra.py) on the published job outputs in results/, and checks
that they give exactly the labels and row grades reported in RESULTS-2.md and
RESULTS-3.md (which were graded by hand from the same rules).

    python -m unittest -v test_graded_results

Standard library only. Made with AI agents working in MV2's name.
"""
import json
import os
import unittest

import grade_cycle2
import grade_tetra

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(*parts):
    with open(os.path.join(HERE, *parts), encoding="utf-8") as f:
        return json.load(f)


class TestResult2(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = _load("results", "prereg2", "results62.json")
        cls.gr = grade_cycle2.grade(cls.res)

    def test_labels(self):
        self.assertEqual(self.gr["status"], "GRADED")
        self.assertEqual(self.gr["decision"], "SUPPORTED")
        self.assertEqual(self.gr["headline"], "PREDICTS_DIVERSITY_AND_SAMPLING")
        self.assertIsNone(self.gr["failure_mode"])
        self.assertTrue(self.gr["decision_agrees_with_job"])

    def test_every_row_holds(self):
        self.assertEqual([r["id"] for r in self.gr["rows"]],
                         ["Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9"])
        self.assertEqual([r["holds"] for r in self.gr["rows"]], [True] * 9)

    def test_brier(self):
        b = self.gr["brier"]
        self.assertEqual(b["n"], 9)
        self.assertAlmostEqual(b["brier"], 2.965 / 9, places=6)

    def test_window(self):
        self.assertEqual(self.gr["window_check"]["n_h6"], 91)
        self.assertEqual(self.gr["window_check"]["n_window_bins"], 91)

    def test_markdown_renders(self):
        md = grade_cycle2.render_md(self.gr, self.res)
        self.assertIn("PREDICTS_DIVERSITY_AND_SAMPLING", md)


class TestResult3(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = _load("results", "prereg3", "results_tetra.json")
        cls.gr = grade_tetra.grade(cls.res)

    def test_labels(self):
        self.assertEqual(self.gr["status"], "GRADED_DATA_LIMITED")
        self.assertEqual(self.gr["decision"], "DATA_LIMITED")
        self.assertEqual(self.gr["headline"], "DATA_LIMITED")
        self.assertTrue(self.gr["decision_agrees_with_job"])

    def test_rows(self):
        got = {r["id"]: r["holds"] for r in self.gr["rows"]}
        want = {"T1": False, "T2": True, "T3": False, "T4": True, "T5": False,
                "T6": True, "T7": True, "T8": True, "T9": True, "T10": True,
                "T11": True}
        self.assertEqual(got, want)

    def test_brier_covers_t1_only(self):
        b = self.gr["brier"]
        self.assertEqual(b["n"], 1)
        self.assertAlmostEqual(b["brier"], 0.49, places=6)

    def test_checks(self):
        ck = self.gr["checks"]
        self.assertEqual((ck["n_bins"], ck["n_h1"], ck["n_h2"]), (38, 38, 30))
        self.assertTrue(ck["data_ok"])
        self.assertTrue(ck["events_ok"])
        self.assertEqual(ck["event_stages_missing_from_window"], ["Changhsingian"])
        self.assertEqual(len(ck["dropped_bins"]), 8)

    def test_markdown_renders(self):
        md = grade_tetra.render_md(self.gr, self.res)
        self.assertIn("GRADED_DATA_LIMITED", md)


if __name__ == "__main__":
    unittest.main()
