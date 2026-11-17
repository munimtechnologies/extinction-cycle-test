# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""grade_cycle2.py: the mechanical reading of PREREG-2.md (the out-of-sample
test of the 62-Myr marine diversity cycle, ephemeris maxima at 27.8 + 62 n Ma).

Written 2026-11-17 (about 19:40 UTC), BEFORE the job extinction-cycle-62-v1
had produced any output. See PRE-RESULT-NOTE-2.md for what had been seen.

It turns the job's results62.json into the graded table (Q1-Q9), the
registered H6 decision and the registered headline. No new statistic is
computed. Reading choices fixed now:

  * Q2/Q3/Q7/Q8 read the "r" and "p_ar1" of H6_range_through, short_lived and
    log_collections. "r > 0 and p_ar1 < 0.05" needs both numbers finite.
  * Q4 is TRUE if H6_range_through.phase_offset_myr <= 8.0 (finite).
  * Q5 is TRUE if H6_periodogram.global_peak_in_band is True (band 58-66 Myr).
  * Q6 is TRUE if H6_periodogram.p_band_ar1 < 0.05 (finite).
  * Q9 is TRUE if H7_mu3t_declining.p_random_phase < 0.05 (finite).
  * If the sanity gate S2 fails, rows are shown but not graded; the Brier
    score then covers Q1 only.

Headline (PREREG-2 "Predictions"):
  SANITY_FAILED                       Q1 fails
  PREDICTS_DIVERSITY                  Q3 holds, Q8 fails
  PREDICTS_DIVERSITY_AND_SAMPLING     Q3 holds and Q8 holds (the same rhythm
                                      is in the number of collections)
  NOT_SUPPORTED                       Q3 fails
For NOT_SUPPORTED, "failure_mode" says whether the period or only the phase
failed: PERIOD_FAILED if neither Q5 nor Q6 holds; PHASE_FAILED if Q5 or Q6
holds but Q4 fails; WEAK if Q5 or Q6 holds and Q4 holds.

Checks (never change the label): the decision is compared with the job's
decision_H6 (CHECK_CODE on disagreement), and H6_range_through.n with the
number of window_bins (CHECK_DATA on a mismatch).

    python grade_cycle2.py --results results62.json --out-dir .

Outputs graded2.json and GRADED-2.md. Standard library only.
Made with AI agents working in MV2's name.
"""
import argparse
import json
import math
import os

ALPHA = 0.05
PHASE_TOL = 8.0

PRIORS = [
    ("Q1", 0.85, "Sanity gate S2 passes (richness drops across the end-Permian, end-Cretaceous and end-Ordovician)"),
    ("Q2", 0.70, "H6 r > 0 (range-through richness vs cos(2 pi (t - 27.8)/62))"),
    ("Q3", 0.40, "H6 SUPPORTED: r > 0 and p_AR1 < 0.05"),
    ("Q4", 0.45, "The fitted 62-Myr phase is within +/-8 Myr of 27.8 Ma"),
    ("Q5", 0.30, "The global periodogram peak (20-200 Myr) lies in 58-66 Myr"),
    ("Q6", 0.30, "The 58-66 Myr band p_AR1 < 0.05"),
    ("Q7", 0.40, "Short-lived genera: r > 0 and p_AR1 < 0.05"),
    ("Q8", 0.30, "log(collections): r > 0 and p_AR1 < 0.05 (the cycle is also in sampling)"),
    ("Q9", 0.40, "H7 (extinction higher in declining half-cycles) p_random_phase < 0.05"),
]


def _f(x):
    if isinstance(x, bool):
        return float("nan")
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def _finite(x):
    v = _f(x)
    return not (math.isnan(v) or math.isinf(v))


def _get(d, *keys):
    for k in keys:
        if not isinstance(d, dict):
            return None
        d = d.get(k)
    return d


def _sig_pos(block):
    """r > 0 and p_ar1 < ALPHA, both finite."""
    if not isinstance(block, dict):
        return False
    r, p = block.get("r"), block.get("p_ar1")
    return bool(_finite(r) and _f(r) > 0 and _finite(p) and _f(p) < ALPHA)


def _p_below(p):
    return bool(_finite(p) and _f(p) < ALPHA)


def grade_rows(res):
    g = {k: None for k, _, _ in PRIORS}
    san = res.get("sanity")
    if isinstance(san, dict) and "passed" in san:
        g["Q1"] = bool(san["passed"])
    h6 = res.get("H6_range_through")
    if isinstance(h6, dict):
        r = h6.get("r")
        g["Q2"] = bool(_finite(r) and _f(r) > 0)
        g["Q3"] = _sig_pos(h6)
        off = h6.get("phase_offset_myr")
        g["Q4"] = bool(_finite(off) and _f(off) <= PHASE_TOL)
    pg = res.get("H6_periodogram")
    if isinstance(pg, dict):
        if "global_peak_in_band" in pg:
            g["Q5"] = pg["global_peak_in_band"] is True
        g["Q6"] = _p_below(pg.get("p_band_ar1"))
    if isinstance(res.get("short_lived"), dict):
        g["Q7"] = _sig_pos(res["short_lived"])
    if isinstance(res.get("log_collections"), dict):
        g["Q8"] = _sig_pos(res["log_collections"])
    if isinstance(res.get("H7_mu3t_declining"), dict):
        g["Q9"] = _p_below(res["H7_mu3t_declining"].get("p_random_phase"))
    return g


def decision_from(g):
    if not g.get("Q1"):
        return "SANITY_FAILED"
    return "SUPPORTED" if g.get("Q3") else "NOT_SUPPORTED"


def headline_from(g):
    if not g.get("Q1"):
        return "SANITY_FAILED"
    if g.get("Q3"):
        return "PREDICTS_DIVERSITY_AND_SAMPLING" if g.get("Q8") else "PREDICTS_DIVERSITY"
    return "NOT_SUPPORTED"


def failure_mode(g):
    if not g.get("Q1") or g.get("Q3"):
        return None
    period = bool(g.get("Q5") or g.get("Q6"))
    if not period:
        return "PERIOD_FAILED"
    return "WEAK" if g.get("Q4") else "PHASE_FAILED"


def brier(g, only=None):
    sq, n = 0.0, 0
    for k, p, _ in PRIORS:
        if g.get(k) is None or (only is not None and k not in only):
            continue
        sq += (p - (1.0 if g[k] else 0.0)) ** 2
        n += 1
    return {"n": n, "brier": (sq / n) if n else float("nan"), "brier_coinflip": 0.25}


def check_window(res):
    wb = res.get("window_bins")
    n = _get(res, "H6_range_through", "n")
    nb = len(wb) if isinstance(wb, list) else None
    ok = nb is not None and isinstance(n, int) and not isinstance(n, bool) and n == nb
    return {"n_h6": n, "n_window_bins": nb, "passed": bool(ok)}


def grade(res):
    out = {"source": "grade_cycle2.py (published before results)", "rows": [],
           "job_status": res.get("status") if isinstance(res, dict) else None}
    if not isinstance(res, dict) or res.get("status") != "ok" or "sanity" not in res:
        out["status"] = "JOB_FAILED"
        out["decision"] = out["headline"] = out["failure_mode"] = None
        return out
    g = grade_rows(res)
    out["decision"] = decision_from(g)
    out["headline"] = headline_from(g)
    out["failure_mode"] = failure_mode(g)
    out["decision_from_job"] = res.get("decision_H6")
    out["decision_agrees_with_job"] = out["decision"] == out["decision_from_job"]
    for k, p, text in PRIORS:
        out["rows"].append({"id": k, "prior": p, "prediction": text, "holds": g[k]})
    out["brier"] = brier(g) if g.get("Q1") else brier(g, only=("Q1",))
    out["window_check"] = check_window(res)
    status = "GRADED" if g.get("Q1") else "GRADED_SANITY_FAILED"
    if not out["decision_agrees_with_job"]:
        status = "CHECK_CODE"
    if not out["window_check"]["passed"]:
        status = "CHECK_DATA"
    out["status"] = status
    return out


def _fmt(x, nd=3):
    if x is None:
        return "n/a"
    if isinstance(x, bool):
        return "yes" if x else "no"
    v = _f(x)
    return f"{v:.{nd}f}" if _finite(v) else "n/a"


def render_md(gr, res):
    L = ["# Graded result: extinction-cycle-test, PREREG-2 (62-Myr diversity ephemeris)", "",
         "Produced mechanically by `grade_cycle2.py`, which was published before the result",
         "existed. Made with AI agents working in MV2's name.", "",
         f"- Status: **{gr.get('status')}**",
         f"- H6 decision: **{gr.get('decision')}**",
         f"- Registered headline: **{gr.get('headline')}**"]
    if gr.get("status") == "JOB_FAILED":
        L += ["", f"The job did not finish (status `{gr.get('job_status')}`). Nothing is graded."]
        return "\n".join(L) + "\n"
    if gr.get("failure_mode"):
        L.append(f"- What failed: **{gr['failure_mode']}**")
    L += [f"- Decision written by the job: {gr.get('decision_from_job')} "
          f"(agrees: {_fmt(gr.get('decision_agrees_with_job'))})", ""]
    if gr.get("status") == "GRADED_SANITY_FAILED":
        L += ["The sanity gate failed, so no hypothesis is graded; rows Q2-Q9 are shown for",
              "the record only.", ""]
    L += ["| # | Prediction | Prior | Holds? |", "|---|---|---|---|"]
    for r in gr["rows"]:
        h = "n/a" if r["holds"] is None else ("YES" if r["holds"] else "no")
        L.append(f"| {r['id']} | {r['prediction']} | {int(round(100 * r['prior']))}% | {h} |")
    b = gr["brier"]
    L += ["", f"Brier score of the registered priors: {_fmt(b['brier'])} over {b['n']} rows "
          f"(a 50% forecaster scores 0.250).", "", "## Key numbers", ""]
    for key, name in (("H6_range_through", "range-through richness (cubic detrend)"),
                      ("H6_quadratic_detrend", "range-through richness (quadratic detrend)"),
                      ("short_lived", "short-lived genera"),
                      ("sampled_in_bin", "sampled-in-bin richness"),
                      ("log_collections", "log(collections)")):
        mf = res.get(key) or {}
        L.append(f"- {name}: r {_fmt(mf.get('r'))}, n {mf.get('n', 'n/a')}, "
                 f"p_AR1 {_fmt(mf.get('p_ar1'), 4)}, p_phase {_fmt(mf.get('p_phase'), 4)}, "
                 f"fitted max {_fmt(mf.get('fit_tmax_ma'), 2)} Ma (offset "
                 f"{_fmt(mf.get('phase_offset_myr'), 2)} Myr)")
    for key, name in (("H6_periodogram", "Periodogram (richness)"),
                      ("short_lived_periodogram", "Periodogram (short-lived)")):
        pg = res.get(key) or {}
        L.append(f"- {name}: global peak {_fmt(pg.get('global_peak_period_myr'), 2)} Myr "
                 f"(in 58-66 band: {_fmt(pg.get('global_peak_in_band'))}); band-max p_AR1 "
                 f"{_fmt(pg.get('p_band_ar1'), 4)}")
    jt = res.get("H6_jitter") or {}
    L.append(f"- Timescale jitter: r 5-50-95% {_fmt(jt.get('r_p05'))} / {_fmt(jt.get('r_p50'))}"
             f" / {_fmt(jt.get('r_p95'))}; fraction r > 0 {_fmt(jt.get('frac_r_pos'))}")
    h7 = res.get("H7_mu3t_declining") or {}
    L.append(f"- H7: n {h7.get('n', 'n/a')}, declining {h7.get('n_declining', 'n/a')}, "
             f"difference in mean normal score {_fmt(h7.get('diff_normal_scores'))}, "
             f"p_random_phase {_fmt(h7.get('p_random_phase'), 4)}")
    san = res.get("sanity") or {}
    checks = san.get("checks") or {}
    L.append("- Sanity S2: " + (", ".join(f"{k} {_fmt(v)}" for k, v in sorted(checks.items()))
                                or "n/a"))
    wc = gr.get("window_check") or {}
    L.append(f"- Window check: H6 n {wc.get('n_h6')}, window bins {wc.get('n_window_bins')}")
    return "\n".join(L) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args(argv)
    with open(args.results, encoding="utf-8") as f:
        res = json.load(f)
    gr = grade(res)
    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, "graded2.json"), "w", encoding="utf-8") as f:
        json.dump(gr, f, indent=1, default=str)
    with open(os.path.join(args.out_dir, "GRADED-2.md"), "w", encoding="utf-8") as f:
        f.write(render_md(gr, res))
    print(json.dumps({"status": gr["status"], "decision": gr["decision"],
                      "headline": gr["headline"], "failure_mode": gr["failure_mode"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
