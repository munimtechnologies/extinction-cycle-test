# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""grade_tetra.py: the mechanical reading of PREREG-3.md (+ PREREG-3-AMENDMENT-1.md),
the out-of-sample test of the claimed 27.5-Myr cycle in land-tetrapod extinction
(Rampino, Caldeira & Zhu 2021; ephemeris maxima at 7.0 + 27.5 n Ma).

Written 2026-11-17 (about 20:20 UTC), BEFORE the job extinction-cycle-tetra-v2
had produced any output and before any PBDB tetrapod occurrence data had been
analysed. See PRE-RESULT-NOTE-3.md for what had been seen.

It turns the job's results_tetra.json into the graded table (T1-T11), the
registered decision and the registered headline. No new statistic is
computed. Reading choices fixed now (keys are those written by tetra.analyze):

  * T1  primary.sanity.passed is True.
  * T2  primary.H_rcz.r > 0 (finite).
  * T3  primary.H_rcz: r > 0 and p_ar1 < 0.05 (both finite).
  * T4  primary.leave_events_out.p_ar1 >= 0.05 (finite; read literally as
        registered, "H2 p_ar1 >= 0.05").
  * T5  primary.H_rcz.p_phase < 0.05 (finite).
  * T6  primary.periodogram.global_peak_in_band is False (band 26-29 Myr,
        periodogram 10-100 Myr).
  * T7  primary.periodogram.p_band_ar1 >= 0.05 (finite).
  * T8  primary.H_mb2014.r > 0 (finite).
  * T9  primary.duration_artifact.p_ar1 >= 0.05 (finite).
  * T10 nonmammal.H_rcz.r and primary.H_rcz.r are both finite, non-zero and of
        the same sign. FALSE if the non-mammal series was too short to test.
  * T11 primary.event_test.p_binom >= 0.05 (finite).
  * A row whose block is missing from the results is n/a (None), not graded.
  * NaN never satisfies an inequality: a NaN p-value makes "p < 0.05" FALSE
    and "p >= 0.05" FALSE as well.
  * If T1 fails, rows are shown but not graded; the Brier score covers T1 only.

Decision (PREREG-3 "Decision rule"):
  DATA_LIMITED               T1 fails
  SUPPORTED_BEYOND_EVENTS    H1 and H2 both have r > 0 and p_ar1 < 0.05
  SUPPORTED_ONLY_BY_EVENTS   H1 holds, H2 does not
  NOT_SUPPORTED              otherwise

Headline (PREREG-3: "independently supports" only for SUPPORTED_BEYOND_EVENTS
and only if H7 shows no duration artifact):
  DATA_LIMITED
  INDEPENDENT_SUPPORT              SUPPORTED_BEYOND_EVENTS and T9 holds
  SUPPORT_BUT_DURATION_ARTIFACT    SUPPORTED_BEYOND_EVENTS and T9 fails
  RESTS_ON_FITTED_EVENTS           SUPPORTED_ONLY_BY_EVENTS
  NOT_SUPPORTED                    NOT_SUPPORTED

Checks (never change the label, only the status):
  CHECK_CODE  the decision differs from primary.decision written by the job,
              or a registered event stage is in the window but was not dropped
              by the leave-events-out test.
  CHECK_DATA  primary.H_rcz.n differs from primary.n_bins, or
              leave_events_out.n differs from n_bins minus the dropped stages.

    python grade_tetra.py --results results_tetra.json --out-dir .

Outputs graded3.json and GRADED-3.md. Standard library only.
Made with AI agents working in MV2's name.
"""
import argparse
import json
import math
import os

ALPHA = 0.05

# The nine stages that hold RCZ2021's tetrapod events (PREREG-3, H2).
EVENT_STAGES = ("Tortonian", "Priabonian", "Maastrichtian", "Tithonian", "Rhaetian",
                "Norian", "Carnian", "Changhsingian", "Capitanian")

PRIORS = [
    ("T1", 0.70, "Sanity gate passes (>= 40 window stages with a rate; Maastrichtian in the top quartile)"),
    ("T2", 0.65, "H1 r > 0 (detrended mu vs cos(2 pi (t - 7.0)/27.5))"),
    ("T3", 0.25, "H1 r > 0 and p_AR1 < 0.05"),
    ("T4", 0.85, "H2 (leave the nine event stages out) p_AR1 >= 0.05"),
    ("T5", 0.20, "H1 p_phase < 0.05 (the published phase beats random phases)"),
    ("T6", 0.80, "The global periodogram peak (10-100 Myr) is outside 26-29 Myr"),
    ("T7", 0.75, "H4 26-29 Myr band p_AR1 >= 0.05"),
    ("T8", 0.60, "H6 marine ephemeris (9.5 + 27.16 n Ma) r > 0"),
    ("T9", 0.80, "H7 stage-duration artifact p_AR1 >= 0.05"),
    ("T10", 0.65, "Non-mammal series has the same sign of r as H1"),
    ("T11", 0.70, "H5 event test p_binom >= 0.05"),
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


def _lt(x, a=ALPHA):
    return bool(_finite(x) and _f(x) < a)


def _ge(x, a=ALPHA):
    return bool(_finite(x) and _f(x) >= a)


def _pos(x):
    return bool(_finite(x) and _f(x) > 0)


def _sig_pos(block):
    """r > 0 and p_ar1 < ALPHA, both finite."""
    if not isinstance(block, dict):
        return False
    return _pos(block.get("r")) and _lt(block.get("p_ar1"))


def _block(d, key):
    v = d.get(key) if isinstance(d, dict) else None
    return v if isinstance(v, dict) else None


def _sign(x):
    if not _finite(x) or _f(x) == 0:
        return 0
    return 1 if _f(x) > 0 else -1


def grade_rows(res):
    g = {k: None for k, _, _ in PRIORS}
    prim = _block(res, "primary") or {}
    san = _block(prim, "sanity")
    if san is not None and "passed" in san:
        g["T1"] = san["passed"] is True
    h1 = _block(prim, "H_rcz")
    if h1 is not None:
        g["T2"] = _pos(h1.get("r"))
        g["T3"] = _sig_pos(h1)
        g["T5"] = _lt(h1.get("p_phase"))
    h2 = _block(prim, "leave_events_out")
    if h2 is not None:
        g["T4"] = _ge(h2.get("p_ar1"))
    pg = _block(prim, "periodogram")
    if pg is not None:
        if "global_peak_in_band" in pg:
            g["T6"] = pg["global_peak_in_band"] is False
        g["T7"] = _ge(pg.get("p_band_ar1"))
    h6 = _block(prim, "H_mb2014")
    if h6 is not None:
        g["T8"] = _pos(h6.get("r"))
    h7 = _block(prim, "duration_artifact")
    if h7 is not None:
        g["T9"] = _ge(h7.get("p_ar1"))
    nm = _block(res, "nonmammal")
    if h1 is not None and nm is not None:
        s1 = _sign(h1.get("r"))
        s2 = _sign((_block(nm, "H_rcz") or {}).get("r"))
        g["T10"] = bool(s1 != 0 and s1 == s2)
    ev = _block(prim, "event_test")
    if ev is not None:
        g["T11"] = _ge(ev.get("p_binom"))
    return g


def decision_from(res, g):
    if not g.get("T1"):
        return "DATA_LIMITED"
    prim = _block(res, "primary") or {}
    if g.get("T3"):
        return ("SUPPORTED_BEYOND_EVENTS" if _sig_pos(_block(prim, "leave_events_out"))
                else "SUPPORTED_ONLY_BY_EVENTS")
    return "NOT_SUPPORTED"


def headline_from(decision, g):
    if decision == "SUPPORTED_BEYOND_EVENTS":
        return "INDEPENDENT_SUPPORT" if g.get("T9") else "SUPPORT_BUT_DURATION_ARTIFACT"
    if decision == "SUPPORTED_ONLY_BY_EVENTS":
        return "RESTS_ON_FITTED_EVENTS"
    return decision


def brier(g, only=None):
    sq, n = 0.0, 0
    for k, p, _ in PRIORS:
        if g.get(k) is None or (only is not None and k not in only):
            continue
        sq += (p - (1.0 if g[k] else 0.0)) ** 2
        n += 1
    return {"n": n, "brier": (sq / n) if n else float("nan"), "brier_coinflip": 0.25}


def _int(x):
    return x if isinstance(x, int) and not isinstance(x, bool) else None


def checks(res):
    prim = _block(res, "primary") or {}
    wb = prim.get("window_bins") if isinstance(prim.get("window_bins"), list) else []
    lo = _block(prim, "leave_events_out") or {}
    dropped = lo.get("dropped_bins") if isinstance(lo.get("dropped_bins"), list) else []
    n_bins = _int(prim.get("n_bins"))
    n_h1 = _int((_block(prim, "H_rcz") or {}).get("n"))
    n_h2 = _int(lo.get("n"))
    out = {"n_bins": n_bins, "n_h1": n_h1, "n_h2": n_h2, "dropped_bins": dropped,
           "event_stages_in_window": [s for s in EVENT_STAGES if s in wb],
           "event_stages_missing_from_window": [s for s in EVENT_STAGES if s not in wb]}
    has_tests = n_h1 is not None
    out["data_ok"] = (not has_tests) or bool(
        n_bins is not None and n_h1 == n_bins
        and (n_h2 is None or n_h2 == n_bins - len(dropped)))
    out["events_ok"] = (not lo) or all(s in dropped for s in out["event_stages_in_window"])
    return out


def grade(res):
    out = {"source": "grade_tetra.py (published before results)", "rows": [],
           "job_status": res.get("status") if isinstance(res, dict) else None}
    if (not isinstance(res, dict) or res.get("status") != "ok"
            or not isinstance(_block(_block(res, "primary") or {}, "sanity"), dict)):
        out["status"] = "JOB_FAILED"
        out["decision"] = out["headline"] = None
        return out
    g = grade_rows(res)
    out["decision"] = decision_from(res, g)
    out["headline"] = headline_from(out["decision"], g)
    out["decision_from_job"] = (_block(res, "primary") or {}).get("decision")
    out["decision_agrees_with_job"] = out["decision"] == out["decision_from_job"]
    for k, p, text in PRIORS:
        out["rows"].append({"id": k, "prior": p, "prediction": text, "holds": g[k]})
    out["brier"] = brier(g) if g.get("T1") else brier(g, only=("T1",))
    out["checks"] = checks(res)
    status = "GRADED" if g.get("T1") else "GRADED_DATA_LIMITED"
    if not out["decision_agrees_with_job"] or not out["checks"]["events_ok"]:
        status = "CHECK_CODE"
    if not out["checks"]["data_ok"]:
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


def _mf_line(name, mf):
    mf = mf or {}
    return (f"- {name}: r {_fmt(mf.get('r'))}, n {mf.get('n', 'n/a')}, "
            f"p_AR1 {_fmt(mf.get('p_ar1'), 4)}, p_phase {_fmt(mf.get('p_phase'), 4)}, "
            f"fitted max {_fmt(mf.get('fit_tmax_ma'), 2)} Ma (offset "
            f"{_fmt(mf.get('phase_offset_myr'), 2)} Myr)")


def render_md(gr, res):
    L = ["# Graded result: extinction-cycle-test, PREREG-3 (land tetrapods, 27.5 Myr)", "",
         "Produced mechanically by `grade_tetra.py`, which was published before the result",
         "existed. Made with AI agents working in MV2's name.", "",
         f"- Status: **{gr.get('status')}**",
         f"- Decision: **{gr.get('decision')}**",
         f"- Registered headline: **{gr.get('headline')}**"]
    if gr.get("status") == "JOB_FAILED":
        L += ["", f"The job did not finish (status `{gr.get('job_status')}`). Nothing is graded."]
        return "\n".join(L) + "\n"
    L += [f"- Decision written by the job: {gr.get('decision_from_job')} "
          f"(agrees: {_fmt(gr.get('decision_agrees_with_job'))})", ""]
    if gr.get("status") == "GRADED_DATA_LIMITED":
        L += ["The sanity gate T1 failed, so nothing is claimed; rows T2-T11 are shown for",
              "the record only.", ""]
    L += ["| # | Prediction | Prior | Holds? |", "|---|---|---|---|"]
    for r in gr["rows"]:
        h = "n/a" if r["holds"] is None else ("YES" if r["holds"] else "no")
        L.append(f"| {r['id']} | {r['prediction']} | {int(round(100 * r['prior']))}% | {h} |")
    b = gr["brier"]
    L += ["", f"Brier score of the registered priors: {_fmt(b['brier'])} over {b['n']} rows "
          f"(a 50% forecaster scores 0.250).", "", "## Key numbers", ""]
    prim = _block(res, "primary") or {}
    L.append(_mf_line("H1, all land tetrapods (7.0 + 27.5 n)", _block(prim, "H_rcz")))
    L.append(_mf_line("H2, nine event stages left out", _block(prim, "leave_events_out")))
    L.append(_mf_line("H3, normal scores", _block(prim, "H_rcz_ranks")))
    L.append(_mf_line("H6, marine ephemeris (9.5 + 27.16 n)", _block(prim, "H_mb2014")))
    L.append(_mf_line("H7, stage durations", _block(prim, "duration_artifact")))
    pg = _block(prim, "periodogram") or {}
    L.append(f"- H4 periodogram: global peak {_fmt(pg.get('global_peak_period_myr'), 2)} Myr "
             f"(in 26-29 band: {_fmt(pg.get('global_peak_in_band'))}); band-max p_AR1 "
             f"{_fmt(pg.get('p_band_ar1'), 4)}")
    ev = _block(prim, "event_test") or {}
    L.append(f"- H5 event test: {ev.get('hits', 'n/a')} of {ev.get('n_peaks', 'n/a')} peaks "
             f"within 3 Myr, p_binom {_fmt(ev.get('p_binom'), 4)}")
    jt = _block(prim, "H_rcz_jitter") or {}
    L.append(f"- Timescale jitter: r 5-50-95% {_fmt(jt.get('r_p05'))} / {_fmt(jt.get('r_p50'))}"
             f" / {_fmt(jt.get('r_p95'))}; fraction r > 0 {_fmt(jt.get('frac_r_pos'))}")
    for key, name in (("nonmammal", "Non-mammal tetrapods"),
                      ("window_260", "Window 2.5-260 Ma"),
                      ("foote_q", "Foote boundary-crosser rate q")):
        blk = _block(res, key) or {}
        L.append(_mf_line(f"{name} (decision {blk.get('decision', 'n/a')})",
                          _block(blk, "H_rcz")))
    san = _block(prim, "sanity") or {}
    cks = san.get("checks") if isinstance(san.get("checks"), dict) else {}
    L.append("- Sanity T1: " + (", ".join(f"{k} {_fmt(v)}" for k, v in sorted(cks.items()))
                                or "n/a"))
    ck = gr.get("checks") or {}
    L.append(f"- Bin check: window bins {ck.get('n_bins')}, H1 n {ck.get('n_h1')}, "
             f"H2 n {ck.get('n_h2')}, dropped {', '.join(ck.get('dropped_bins') or []) or 'none'}")
    miss = ck.get("event_stages_missing_from_window") or []
    if miss:
        L.append(f"- Event stages with no defined rate (not in the window): {', '.join(miss)}")
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
    with open(os.path.join(args.out_dir, "graded3.json"), "w", encoding="utf-8") as f:
        json.dump(gr, f, indent=1, default=str)
    with open(os.path.join(args.out_dir, "GRADED-3.md"), "w", encoding="utf-8") as f:
        f.write(render_md(gr, res))
    print(json.dumps({"status": gr["status"], "decision": gr["decision"],
                      "headline": gr["headline"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
