# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""grade_cycle1.py: the mechanical reading of PREREG.md (PREREG-1, the
out-of-sample test of Melott & Bambach's 27.16-Myr extinction ephemeris).

Written and published 2026-11-17 (about 18:45 UTC), BEFORE the job
extinction-cycle-v1 had produced any result and before any extinction rate
had been computed.

It turns the job's results.json into the graded prediction table (P1-P9),
the registered H1 decision and the registered headline, row by row. No new
statistic is computed. Reading choices fixed now (PRE-RESULT-NOTE-1.md):

  * P3 is graded on the numbers alone (r > 0 and p_ar1 < 0.05 for mu3t at
    MB2014). The H1 decision additionally needs the sanity gate, as in job.py.
  * "p >= 0.05" rows (P4, P5, P7, P9) are TRUE only if the p value is finite.
  * P6 is TRUE if the job reports global_peak_in_band == False.
  * P8 is TRUE if phase_offset_myr <= 3.0 (finite).
  * If the sanity gate fails, the rows are still computed and shown, but
    PREREG.md says no hypothesis is graded: the Brier score then covers P1 only.

Headline (PREREG "Predictions"): PREDICTS_BEYOND_MASS_EXTINCTIONS if P3 holds
and P4 fails; RESTS_ON_MASS_EXTINCTIONS if P3 and P4 both hold;
NOT_SUPPORTED if P3 fails; SANITY_FAILED if P1 fails.

Checks: the label is compared with the job's decision_H1 (CHECK_CODE on
disagreement), and the GTS2020 top ages of the six Big Five bins in
results.json's mu3t_bins are compared with the registered stage table
(CHECK_DATA on a mismatch or a missing bin). Neither changes the label.

    python grade_cycle1.py --results results.json --out-dir .

Outputs graded1.json and GRADED-1.md. Standard library only.
Made with AI agents working in MV2's name.
"""
import argparse
import json
import math
import os

ALPHA = 0.05
PERIOD, T0 = 27.16, 9.5
HIT = 3.0

PRIORS = [
    ("P1", 0.85, "Sanity gate S1 passes (>= 4 of the Big Five among the 10 highest mu3t bins)"),
    ("P2", 0.75, "H1 r > 0 at the MB2014 ephemeris (27.16 Myr, 9.5 Ma)"),
    ("P3", 0.40, "H1 SUPPORTED: r > 0 and p_AR1 < 0.05"),
    ("P4", 0.75, "With the six Big Five bins removed, p_AR1 >= 0.05"),
    ("P5", 0.65, "On normal scores (ranks), p_AR1 >= 0.05"),
    ("P6", 0.75, "The global periodogram peak (10-100 Myr) is NOT in 26-28.5 Myr"),
    ("P7", 0.70, "The band-max p_AR1 (H2) >= 0.05"),
    ("P8", 0.45, "The fitted phase at 27.16 Myr is within +/-3 Myr of 9.5 Ma"),
    ("P9", 0.60, "H3 binomial p >= 0.05"),
]

# GTS2020 top ages (Ma) of the Big Five bins, from the registered stage table
# (cycle.STAGES_GTS2020), and the value of the MB2014 template there. Public
# inputs only; see PRE-RESULT-NOTE-1.md.
BIG_FIVE_TOPS_MA = {"Katian": 445.2, "Hirnantian": 443.1, "Frasnian": 372.15,
                    "Changhsingian": 251.902, "Rhaetian": 201.4, "Maastrichtian": 66.0}
TOP_TOL = 0.001


def template(t, period=PERIOD, t0=T0):
    return math.cos(2 * math.pi * (t - t0) / period)


def offset(t, period=PERIOD, t0=T0):
    """Distance (Myr) from age t to the nearest predicted maximum."""
    return abs(((t - t0 + period / 2.0) % period) - period / 2.0)


def _f(x):
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


def _p_at_least(p):
    return bool(_finite(p) and _f(p) >= ALPHA)


def grade_rows(res):
    g = {k: None for k, _, _ in PRIORS}
    san = res.get("sanity")
    if isinstance(san, dict) and "passed" in san:
        g["P1"] = bool(san["passed"])
    mf = _get(res, "mu3t", "MB2014", "matched_filter")
    if not isinstance(mf, dict):
        return g
    r, p = _f(mf.get("r")), _f(mf.get("p_ar1"))
    g["P2"] = bool(_finite(r) and r > 0)
    g["P3"] = bool(_finite(r) and r > 0 and _finite(p) and p < ALPHA)
    g["P4"] = _p_at_least(_get(res, "mu3t_without_big_five", "MB2014", "matched_filter", "p_ar1"))
    g["P5"] = _p_at_least(_get(res, "mu3t_normal_scores", "MB2014", "matched_filter", "p_ar1"))
    pg = _get(res, "mu3t", "periodogram") or {}
    if "global_peak_in_band" in pg:
        g["P6"] = not bool(pg["global_peak_in_band"])
    g["P7"] = _p_at_least(pg.get("p_band_ar1"))
    off = _f(mf.get("phase_offset_myr"))
    g["P8"] = bool(_finite(off) and off <= HIT)
    g["P9"] = _p_at_least(_get(res, "mu3t", "MB2014", "events", "p_binom"))
    return g


def decision_from(g):
    if not g.get("P1"):
        return "SANITY_FAILED"
    return "SUPPORTED" if g.get("P3") else "NOT_SUPPORTED"


def headline_from(g):
    if not g.get("P1"):
        return "SANITY_FAILED"
    if not g.get("P3"):
        return "NOT_SUPPORTED"
    return "RESTS_ON_MASS_EXTINCTIONS" if g.get("P4") else "PREDICTS_BEYOND_MASS_EXTINCTIONS"


def brier(g, only=None):
    sq, n = 0.0, 0
    for k, p, _ in PRIORS:
        if g.get(k) is None or (only is not None and k not in only):
            continue
        sq += (p - (1.0 if g[k] else 0.0)) ** 2
        n += 1
    return {"n": n, "brier": (sq / n) if n else float("nan"), "brier_coinflip": 0.25}


def check_tops(res):
    bins = res.get("mu3t_bins")
    seen = {}
    if isinstance(bins, list):
        for b in bins:
            if isinstance(b, dict) and b.get("bin") in BIG_FIVE_TOPS_MA:
                seen[b["bin"]] = _f(b.get("top_ma"))
    bad = [{"bin": k, "job": v, "registered": BIG_FIVE_TOPS_MA[k]}
           for k, v in sorted(seen.items())
           if not (_finite(v) and abs(v - BIG_FIVE_TOPS_MA[k]) <= TOP_TOL)]
    missing = sorted(set(BIG_FIVE_TOPS_MA) - set(seen))
    return {"mismatches": bad, "missing_bins": missing, "passed": not bad and not missing}


def grade(res):
    out = {"source": "grade_cycle1.py (published before results)", "rows": [],
           "job_status": res.get("status") if isinstance(res, dict) else None}
    if not isinstance(res, dict) or res.get("status") != "ok" or "sanity" not in res:
        out["status"] = "JOB_FAILED"
        out["decision"] = out["headline"] = None
        return out
    g = grade_rows(res)
    out["decision"] = decision_from(g)
    out["headline"] = headline_from(g)
    out["decision_from_job"] = res.get("decision_H1")
    out["decision_agrees_with_job"] = out["decision"] == out["decision_from_job"]
    for k, p, text in PRIORS:
        out["rows"].append({"id": k, "prior": p, "prediction": text, "holds": g[k]})
    out["brier"] = brier(g) if g.get("P1") else brier(g, only=("P1",))
    out["tops_check"] = check_tops(res)
    status = "GRADED" if g.get("P1") else "GRADED_SANITY_FAILED"
    if not out["decision_agrees_with_job"]:
        status = "CHECK_CODE"
    if not out["tops_check"]["passed"]:
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
    L = ["# Graded result: extinction-cycle-test, PREREG-1 (27.16-Myr ephemeris)", "",
         "Produced mechanically by `grade_cycle1.py`, which was published before the result",
         "existed. Made with AI agents working in MV2's name.", "",
         f"- Status: **{gr.get('status')}**",
         f"- H1 decision: **{gr.get('decision')}**",
         f"- Registered headline: **{gr.get('headline')}**"]
    if gr.get("status") == "JOB_FAILED":
        L += ["", f"The job did not finish (status `{gr.get('job_status')}`). Nothing is graded."]
        return "\n".join(L) + "\n"
    L += [f"- Decision written by the job: {gr.get('decision_from_job')} "
          f"(agrees: {_fmt(gr.get('decision_agrees_with_job'))})", ""]
    if gr.get("status") == "GRADED_SANITY_FAILED":
        L += ["The sanity gate failed, so (PREREG.md) no hypothesis is graded; rows P2-P9 are",
              "shown for the record only.", ""]
    L += ["| # | Prediction | Prior | Holds? |", "|---|---|---|---|"]
    for r in gr["rows"]:
        h = "n/a" if r["holds"] is None else ("YES" if r["holds"] else "no")
        L.append(f"| {r['id']} | {r['prediction']} | {int(round(100 * r['prior']))}% | {h} |")
    b = gr["brier"]
    L += ["", f"Brier score of the registered priors: {_fmt(b['brier'])} over {b['n']} rows "
          f"(a 50% forecaster scores 0.250).", "", "## Key numbers", ""]
    for blk, name in (("mu3t", "mu3t"), ("mu3t_without_big_five", "mu3t without Big Five"),
                      ("mu3t_normal_scores", "mu3t normal scores"), ("foote_q", "Foote q"),
                      ("prop_ext", "proportion last seen")):
        for eph in ("MB2014", "MB2017", "RC2015"):
            mf = _get(res, blk, eph, "matched_filter") or {}
            ev = _get(res, blk, eph, "events") or {}
            L.append(f"- {name}, {eph}: r {_fmt(mf.get('r'))}, n {mf.get('n', 'n/a')}, "
                     f"p_AR1 {_fmt(mf.get('p_ar1'), 4)}, p_phase {_fmt(mf.get('p_phase'), 4)}, "
                     f"fitted max {_fmt(mf.get('fit_tmax_ma'), 2)} Ma (offset "
                     f"{_fmt(mf.get('phase_offset_myr'), 2)} Myr); events {ev.get('hits', 'n/a')}"
                     f"/{ev.get('n_peaks', 'n/a')}, p_binom {_fmt(ev.get('p_binom'), 4)}")
    pg = _get(res, "mu3t", "periodogram") or {}
    L.append(f"- Periodogram (mu3t): global peak {_fmt(pg.get('global_peak_period_myr'), 2)} Myr "
             f"(in 26-28.5 band: {_fmt(pg.get('global_peak_in_band'))}); band-max p_AR1 "
             f"{_fmt(pg.get('p_band_ar1'), 4)}")
    jt = _get(res, "mu3t", "jitter_MB2014") or {}
    L.append(f"- Timescale jitter (H5): r 5-50-95% {_fmt(jt.get('r_p05'))} / "
             f"{_fmt(jt.get('r_p50'))} / {_fmt(jt.get('r_p95'))}; fraction r > 0 "
             f"{_fmt(jt.get('frac_r_pos'))}")
    san = res.get("sanity") or {}
    L.append(f"- Sanity: Big Five events in top 10: {san.get('big_five_in_top10', 'n/a')}; "
             f"top 10 bins: {', '.join(san.get('top10_mu3t_bins') or [])}")
    tc = gr.get("tops_check") or {}
    L += ["", "## Big Five top-age check", "",
          f"mismatches: {len(tc.get('mismatches', []))}; missing: "
          f"{', '.join(tc.get('missing_bins', [])) or 'none'}"]
    for m in tc.get("mismatches", []):
        L.append(f"- {m['bin']}: job {_fmt(m['job'], 3)}, registered {m['registered']}")
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
    with open(os.path.join(args.out_dir, "graded1.json"), "w", encoding="utf-8") as f:
        json.dump(gr, f, indent=1, default=str)
    with open(os.path.join(args.out_dir, "GRADED-1.md"), "w", encoding="utf-8") as f:
        f.write(render_md(gr, res))
    print(json.dumps({"status": gr["status"], "decision": gr["decision"],
                      "headline": gr["headline"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
