# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""PREREG-3: an occurrence-based test of the claimed 27.5-Myr cycle in the
extinction of non-marine (land) tetrapods (Rampino, Caldeira & Zhu 2021,
Historical Biology 33:3084; Geoscience Frontiers 12:101245).

Pure functions (no network); job_tetra.py does the downloads. The statistics
are the ones registered for marine extinction in PREREG-1 (cycle.py).
Made with AI agents working in MV2's name.

Amendment 1 (2026-11-17, before any tetrapod data were downloaded): the
registered query value envtype=terr is rejected by the PBDB API ("bad value
'terr' for parameter 'envtype'", zero records). The API's code for "Any
terrestrial environment" is envtype=terrestrial, which is used below. See
PREREG-3-AMENDMENT-1.md.
"""
import numpy as np

import cycle

# PRIMARY ephemeris. Rampino, Caldeira & Zhu (2021, Geoscience Frontiers
# 12:101245): period 27.5 Myr; "the most recent cluster of geological events
# was approximately 7 million years ago" (NYU press release of the paper);
# moving-window pulse ages below. t0 = 7.0 Ma.
RCZ2021 = {"period": 27.5, "t0": 7.0,
           "source": "Rampino, Caldeira & Zhu 2021, Geoscience Frontiers 12:101245 "
                     "(27.5-Myr pulse, most recent ~7 Ma) and Historical Biology "
                     "33:3084 (27.5-Myr cycle in non-marine tetrapod extinctions)"}
# Secondary: the marine ephemeris of PREREG-1 (cross-realm coherence).
MB2014 = cycle.EPHEMERIDES["MB2014"]

# 'Peaks are centered at about 7, 35, 61, 95, 124, 140, 186, 196, 213, and
# 255 Ma' (RCZ2021, 10-Myr moving window).
PULSES_RCZ2021 = (7.0, 35.0, 61.0, 95.0, 124.0, 140.0, 186.0, 196.0, 213.0, 255.0)

# The nine non-marine tetrapod extinction events of RCZ2021 Table 1 (Ma).
TETRAPOD_EVENTS = (
    ("Late Miocene", 7.25), ("Eocene-Oligocene", 33.9), ("Cretaceous-Paleogene", 66.0),
    ("Jurassic-Cretaceous", 143.1), ("Triassic-Jurassic", 201.4), ("mid-Norian", 215.0),
    ("Carnian-Norian", 227.3), ("Permian-Triassic", 251.9), ("end-Guadalupian", 259.8),
)

WINDOW_T = (2.5, 300.0)     # bin-top ages (Ma) in the primary test ("last 300 My")
WINDOW_260 = (2.5, 260.0)   # secondary window ("last 260 Myr")
BAND_T = (26.0, 29.0)       # periodogram band around 27.5 Myr
EVENT_TOL = 1.5             # Myr: an event belongs to a bin whose top is this close
MIN_BINS = 40               # sanity: bins with a defined rate in the window
MIN_RUN = 12                # below this many bins nothing is computed

PBDB = "https://paleobiodb.org/data1.2/occs/list.csv"
COMMON = "envtype=terrestrial&idreso=genus&idqual=genus_certain&pres=regular&show=class"
QUERIES = {
    "tetra_nonmammal": f"{PBDB}?base_name=Tetrapoda%5EMammalia&{COMMON}",
    "tetra_mammalia": f"{PBDB}?base_name=Mammalia&{COMMON}",
}


def event_bins(stages, ages=None, tol=EVENT_TOL):
    """Stage indices holding the events: the bin whose top is nearest the event
    age if within `tol` Myr, else the bin containing the event age."""
    if ages is None:
        ages = [a for _, a in TETRAPOD_EVENTS]
    out = []
    for a in ages:
        best, bd = None, None
        for i, s in enumerate(stages):
            d = abs(s[2] - a)
            if d <= tol and (bd is None or d < bd):
                best, bd = i, d
        if best is None:
            for i, s in enumerate(stages):
                if s[2] < a <= s[1]:
                    best = i
                    break
        if best is not None and best not in out:
            out.append(best)
    return sorted(out)


def circular_phase(ages, period):
    """Circular mean of ages modulo period (in [0, P)) and resultant length R."""
    th = 2 * np.pi * np.asarray(ages, float) / period
    c, s = float(np.mean(np.cos(th))), float(np.mean(np.sin(th)))
    return float((np.arctan2(s, c) % (2 * np.pi)) * period / (2 * np.pi)), float(np.hypot(c, s))


def analyze(stages, mu, n_sur, rng, window=WINDOW_T, full=True):
    """Run the registered tests on one extinction-rate series (one value per
    stage, at the stage TOP, NaN where undefined).

    Returns (results, (freqs, power)); the periodogram arrays are None when
    nothing was computed or full is False.
    """
    names = [s[0] for s in stages]
    tops = np.array([s[2] for s in stages], float)
    mu = np.asarray(mu, float)
    lo, hi = window
    k = [i for i in range(len(stages)) if lo <= tops[i] <= hi and np.isfinite(mu[i])]
    t = tops[k]
    y = mu[k]
    res = {"window": list(window), "n_bins": len(k), "window_bins": [names[i] for i in k]}
    im = names.index("Maastrichtian")
    maas = bool(len(k) > 0 and im in k and y[k.index(im)] >= np.quantile(y, 0.75))
    checks = {"maastrichtian_top_quartile": maas, "enough_bins": len(k) >= MIN_BINS}
    res["sanity"] = {"checks": checks, "passed": all(checks.values())}
    if len(k) < MIN_RUN:
        res["decision"] = "DATA_LIMITED"
        return res, (None, None)

    P, T0 = RCZ2021["period"], RCZ2021["t0"]
    res["H_rcz"] = cycle.matched_filter_test(t, y, P, T0, n_sur, rng)
    ev = set(event_bins(stages))
    k2 = np.array([j for j, i in enumerate(k) if i not in ev], int)
    res["leave_events_out"] = cycle.matched_filter_test(t[k2], y[k2], P, T0, n_sur, rng)
    res["leave_events_out"]["dropped_bins"] = [names[i] for i in k if i in ev]
    res["H_mb2014"] = cycle.matched_filter_test(t, y, MB2014["period"], MB2014["t0"], n_sur, rng)
    pf = pp = None
    if full:
        res["H_rcz_ranks"] = cycle.matched_filter_test(t, cycle.normal_scores(y), P, T0,
                                                       n_sur, rng)
        res["H_rcz_jitter"] = cycle.jitter_r(t, y, P, T0, 500, rng)
        pg = cycle.periodogram_test(t, y, n_sur, rng, band=BAND_T)
        pf, pp = pg.pop("freqs"), pg.pop("power")
        res["periodogram"] = pg
        res["event_test"] = cycle.event_test(t, y, P, T0)
        dur = np.array([stages[i][1] - stages[i][2] for i in k], float)
        res["duration_artifact"] = cycle.matched_filter_test(t, dur, P, T0, n_sur, rng)

    h = res["H_rcz"]
    lo_ = res["leave_events_out"]
    if not res["sanity"]["passed"]:
        res["decision"] = "DATA_LIMITED"
    elif h["r"] > 0 and h["p_ar1"] < 0.05:
        res["decision"] = ("SUPPORTED_BEYOND_EVENTS" if (lo_["r"] > 0 and lo_["p_ar1"] < 0.05)
                           else "SUPPORTED_ONLY_BY_EVENTS")
    else:
        res["decision"] = "NOT_SUPPORTED"
    return res, (pf, pp)
