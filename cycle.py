# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Core functions for a pre-registered, out-of-sample test of the claimed
~27-Myr periodicity in marine extinction (Raup & Sepkoski 1984; Melott &
Bambach 2010-2014).

Everything here is pure (no network); job.py does the downloads.
Made with AI agents working in MV2's name.
"""
import bisect
import math
from statistics import NormalDist

import numpy as np

# International stages, GTS2020 ages as served by the Paleobiology Database
# (data1.2/intervals/list.txt?scale=1&type=age, read 2026-11-17).
# (name, base_ma, top_ma), oldest first. PBDB serves no Pridoli "age", so the
# Pridoli series (422.7-419.62 Ma) is added as one bin; all Quaternary ages
# are merged into one bin.
STAGES_GTS2020 = [
    ("Fortunian", 538.8, 529.0), ("Stage 2", 529.0, 521.0),
    ("Stage 3", 521.0, 514.5), ("Stage 4", 514.5, 506.5),
    ("Wuliuan", 506.5, 504.5), ("Drumian", 504.5, 500.5),
    ("Guzhangian", 500.5, 497.0), ("Paibian", 497.0, 494.2),
    ("Jiangshanian", 494.2, 491.0), ("Stage 10", 491.0, 486.85),
    ("Tremadocian", 486.85, 477.1), ("Floian", 477.1, 471.3),
    ("Dapingian", 471.3, 469.4), ("Darriwilian", 469.4, 458.2),
    ("Sandbian", 458.2, 452.8), ("Katian", 452.8, 445.2),
    ("Hirnantian", 445.2, 443.1), ("Rhuddanian", 443.1, 440.5),
    ("Aeronian", 440.5, 438.6), ("Telychian", 438.6, 432.9),
    ("Sheinwoodian", 432.9, 430.6), ("Homerian", 430.6, 426.7),
    ("Gorstian", 426.7, 425.0), ("Ludfordian", 425.0, 422.7),
    ("Pridoli", 422.7, 419.62), ("Lochkovian", 419.62, 413.02),
    ("Pragian", 413.02, 410.62), ("Emsian", 410.62, 393.47),
    ("Eifelian", 393.47, 387.95), ("Givetian", 387.95, 382.31),
    ("Frasnian", 382.31, 372.15), ("Famennian", 372.15, 358.86),
    ("Tournaisian", 358.86, 346.7), ("Visean", 346.7, 330.3),
    ("Serpukhovian", 330.3, 323.4), ("Bashkirian", 323.4, 315.2),
    ("Moscovian", 315.2, 307.0), ("Kasimovian", 307.0, 303.7),
    ("Gzhelian", 303.7, 298.9), ("Asselian", 298.9, 293.52),
    ("Sakmarian", 293.52, 290.1), ("Artinskian", 290.1, 283.3),
    ("Kungurian", 283.3, 274.4), ("Roadian", 274.4, 266.9),
    ("Wordian", 266.9, 264.28), ("Capitanian", 264.28, 259.857),
    ("Wuchiapingian", 259.857, 254.14), ("Changhsingian", 254.14, 251.902),
    ("Induan", 251.902, 250.8), ("Olenekian", 250.8, 247.0),
    ("Anisian", 247.0, 241.464), ("Ladinian", 241.464, 237.0),
    ("Carnian", 237.0, 227.3), ("Norian", 227.3, 205.7),
    ("Rhaetian", 205.7, 201.4), ("Hettangian", 201.4, 199.5),
    ("Sinemurian", 199.5, 192.9), ("Pliensbachian", 192.9, 184.2),
    ("Toarcian", 184.2, 174.7), ("Aalenian", 174.7, 170.9),
    ("Bajocian", 170.9, 168.2), ("Bathonian", 168.2, 165.3),
    ("Callovian", 165.3, 161.5), ("Oxfordian", 161.5, 154.8),
    ("Kimmeridgian", 154.8, 149.2), ("Tithonian", 149.2, 143.1),
    ("Berriasian", 143.1, 137.05), ("Valanginian", 137.05, 132.6),
    ("Hauterivian", 132.6, 125.77), ("Barremian", 125.77, 121.4),
    ("Aptian", 121.4, 113.2), ("Albian", 113.2, 100.5),
    ("Cenomanian", 100.5, 93.9), ("Turonian", 93.9, 89.8),
    ("Coniacian", 89.8, 85.7), ("Santonian", 85.7, 83.6),
    ("Campanian", 83.6, 72.2), ("Maastrichtian", 72.2, 66.0),
    ("Danian", 66.0, 61.66), ("Selandian", 61.66, 59.24),
    ("Thanetian", 59.24, 56.0), ("Ypresian", 56.0, 48.07),
    ("Lutetian", 48.07, 41.03), ("Bartonian", 41.03, 37.71),
    ("Priabonian", 37.71, 33.9), ("Rupelian", 33.9, 27.3),
    ("Chattian", 27.3, 23.04), ("Aquitanian", 23.04, 20.45),
    ("Burdigalian", 20.45, 15.98), ("Langhian", 15.98, 13.82),
    ("Serravallian", 13.82, 11.63), ("Tortonian", 11.63, 7.246),
    ("Messinian", 7.246, 5.333), ("Zanclean", 5.333, 3.6),
    ("Piacenzian", 3.6, 2.58), ("Quaternary", 2.58, 0.0),
]

# Bins whose top is one of the "Big Five" (end-Ordovician has two pulses).
BIG_FIVE_BINS = ("Katian", "Hirnantian", "Frasnian", "Changhsingian",
                 "Rhaetian", "Maastrichtian")
BIG_FIVE_EVENTS = (("Katian", "Hirnantian"), ("Frasnian",), ("Changhsingian",),
                   ("Rhaetian",), ("Maastrichtian",))

# Published ephemerides (period P in Myr, t0 = age in Ma of the most recent
# maximum of extinction). MB2014 is the PRIMARY hypothesis.
EPHEMERIDES = {
    "MB2014": {"period": 27.16, "t0": 9.5,
               "source": "Melott & Bambach 2014, Paleobiology 40:177-196 "
                         "(arXiv:1310.4712), window test: '27.16-Myr spacing, "
                         "with the most recent peak at 9.5 Ma'"},
    "MB2017": {"period": 27.19, "t0": 11.3,
               "source": "Melott & Bambach comment, Palaeontology "
                         "doi:10.1111/pala.12322 (fitted sine, as quoted in "
                         "search snippets; secondary)"},
    "RC2015": {"period": 27.0, "t0": 11.8,
               "source": "Rampino & Caldeira 2015, MNRAS 454:3480 "
                         "(extinctions + craters; secondary)"},
}

WINDOW = (2.5, 475.0)          # bin top ages (Ma) used in the tests
HIT_HALF_WIDTH = 3.0           # Myr, as in Melott & Bambach

# Melott & Bambach 2014 Table 1 (19 post-Cambrian events), GTS2012 age, and the
# GTS2020 age of the same boundary/event used here (see PREREG.md, part A).
MB2014_EVENTS = [
    ("End-Maastrichtian", 66.00, 66.0), ("Rhaetian event", 201.30, 201.4),
    ("Changhsingian event", 252.20, 251.902),
    ("Late Frasnian event", 372.20, 372.15),
    ("Late Hirnantian event", 443.80, 443.1),
    ("Katian/Hirnantian event", 445.20, 445.2),
    ("Pliocene events", 2.59, 2.58), ("Late Eocene", 33.90, 33.9),
    ("Late Cenomanian event", 93.90, 93.9), ("Late Tithonian", 145.00, 143.1),
    ("Early Toarcian event", 182.00, 183.0),
    ("Capitanian event", 259.80, 259.857),
    ("Mid-Serpukhovian event", 327.05, 326.85),
    ("Late Famennian events", 358.90, 358.86),
    ("Late Givetian events", 382.70, 382.31),
    ("Middle Miocene", 11.63, 11.63), ("Spathian", 247.10, 247.0),
    ("Late Moscovian", 307.00, 307.0), ("Mid-Caradocian", 452.14, 451.9),
]


# ----------------------------------------------------------------- binning
class BinAssigner:
    """Assign an occurrence age range to one stage.

    The stage is the one containing the range midpoint; the range must lie
    inside that stage within `tol` Myr on each side, else None.
    """

    def __init__(self, stages=STAGES_GTS2020, tol=0.1):
        self.stages = stages
        self.n = len(stages)
        self.tol = tol
        self.bases = [s[1] for s in stages]
        self.tops = [s[2] for s in stages]
        self._tops_asc = self.tops[::-1]

    def __call__(self, max_ma, min_ma):
        if max_ma is None or min_ma is None:
            return None
        if max_ma < min_ma:
            max_ma, min_ma = min_ma, max_ma
        m = 0.5 * (max_ma + min_ma)
        if m < 0 or m >= self.bases[0]:
            return None
        k = bisect.bisect_right(self._tops_asc, m) - 1
        if k < 0:
            return None
        i = self.n - 1 - k
        if max_ma <= self.bases[i] + self.tol and min_ma >= self.tops[i] - self.tol:
            return i
        return None


def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def genus_key(row):
    """'phylum|Genus' for a PBDB occurrence row (pbdb vocabulary), or None."""
    g = (row.get("genus") or "").strip()
    if not g or g.upper().startswith("NO_"):
        rank = (row.get("accepted_rank") or "").strip()
        name = (row.get("accepted_name") or "").strip()
        g = name.split()[0] if (name and rank in (
            "genus", "subgenus", "species", "subspecies")) else ""
    if not g:
        return None
    g = g.split()[0]
    ph = (row.get("phylum") or "").strip()
    return ph + "|" + g


class Accumulator:
    def __init__(self, n):
        self.sets = [set() for _ in range(n)]
        self.n_occ = [0] * n
        self._colls = [set() for _ in range(n)]
        self.stats = {"rows": 0, "assigned": 0, "no_genus": 0, "no_bin": 0}

    def add_rows(self, rows, assigner):
        for row in rows:
            self.stats["rows"] += 1
            key = genus_key(row)
            if key is None:
                self.stats["no_genus"] += 1
                continue
            b = assigner(to_float(row.get("max_ma")), to_float(row.get("min_ma")))
            if b is None:
                self.stats["no_bin"] += 1
                continue
            self.stats["assigned"] += 1
            self.sets[b].add(key)
            self.n_occ[b] += 1
            self._colls[b].add(row.get("collection_no"))

    @property
    def n_coll(self):
        return [len(c) for c in self._colls]


# ------------------------------------------------------------------- rates
def three_timer(sets):
    """Alroy (2008) three-timer extinction rate at the TOP of each bin.

    mu_i = ln(2t_i / 3t_i) + ln(Ps_{i+1}), with 2t_i = |S_{i-1} & S_i|,
    3t_i = |S_{i-1} & S_i & S_{i+1}|, Ps_j = 3t_j / (3t_j + pt_j),
    pt_j = |(S_{j-1} & S_{j+1}) - S_j|. Bins are ordered oldest first.
    Returns (mu, ps, counts) with NaN where undefined.
    """
    n = len(sets)
    nan = float("nan")
    ps = [nan] * n
    mu = [nan] * n
    two = [0] * n
    three = [0] * n
    part = [0] * n
    for j in range(1, n - 1):
        a, b, c = sets[j - 1], sets[j], sets[j + 1]
        ab = a & b
        two[j] = len(ab)
        three[j] = len(ab & c)
        part[j] = len((a & c) - b)
        if three[j] + part[j] > 0:
            ps[j] = three[j] / (three[j] + part[j])
    for i in range(1, n - 1):
        p = ps[i + 1]
        if two[i] > 0 and three[i] > 0 and p == p and p > 0:
            mu[i] = math.log(two[i] / three[i]) + math.log(p)
    return mu, ps, {"two_timers": two, "three_timers": three, "part_timers": part}


def foote_rates(sets):
    """Range-through boundary-crosser metrics (Foote 2000), oldest first.

    q_i = -ln(Nbt / (Nbt + NbL)); prop_i = genera last seen in bin / all genera
    ranging through the bin (Melott & Bambach's 'magnitude').
    """
    n = len(sets)
    first, last = {}, {}
    for i, s in enumerate(sets):
        for k in s:
            if k not in first:
                first[k] = i
            last[k] = i
    nbt = [0] * n
    nbl = [0] * n
    nft = [0] * n
    nfl = [0] * n
    for k, f in first.items():
        l = last[k]
        for i in range(f, l + 1):
            if f < i and l > i:
                nbt[i] += 1
            elif f < i and l == i:
                nbl[i] += 1
            elif f == i and l > i:
                nft[i] += 1
            else:
                nfl[i] += 1
    nan = float("nan")
    q = [nan] * n
    prop = [nan] * n
    for i in range(n):
        if nbt[i] > 0:
            q[i] = -math.log(nbt[i] / (nbt[i] + nbl[i]))
        tot = nbt[i] + nbl[i] + nft[i] + nfl[i]
        if tot > 0:
            prop[i] = (nbl[i] + nfl[i]) / tot
    return q, prop, {"Nbt": nbt, "NbL": nbl, "NFt": nft, "NFL": nfl}


# ----------------------------------------------------------- statistics
def detrend(t, y):
    A = np.vstack([np.ones_like(t), t]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ coef


def detrend_rows(t, Y):
    A = np.vstack([np.ones_like(t), t]).T
    P = A @ np.linalg.pinv(A)
    return Y - Y @ P


def corr_rows(Y, c):
    """Pearson correlation of each row of Y with vector c."""
    Y = np.atleast_2d(np.asarray(Y, float))
    c = np.asarray(c, float) - np.mean(c)
    c = c / np.linalg.norm(c)
    Yc = Y - Y.mean(axis=1, keepdims=True)
    return (Yc @ c) / np.linalg.norm(Yc, axis=1)


def template(t, period, t0):
    return np.cos(2 * np.pi * (np.asarray(t, float) - t0) / period)


def circ_dist(t, period, t0):
    """Distance (Myr) from age t to the nearest predicted maximum."""
    return np.abs(np.mod(np.asarray(t, float) - t0 + period / 2.0, period) - period / 2.0)


def lag1(y):
    y = np.asarray(y, float) - np.mean(y)
    return float(np.sum(y[1:] * y[:-1]) / np.sum(y * y))


def ar1_surrogates(n, phi, m, rng):
    phi = max(min(phi, 0.95), -0.95)
    e = rng.standard_normal((m, n))
    x = np.empty((m, n))
    x[:, 0] = e[:, 0]
    s = math.sqrt(1.0 - phi * phi)
    for k in range(1, n):
        x[:, k] = phi * x[:, k - 1] + s * e[:, k]
    return x


def fit_phase(t, y, period):
    """Least-squares sinusoid at fixed period; returns (amplitude, t_max in [0, P))."""
    w = 2 * np.pi / period
    A = np.vstack([np.ones_like(t), np.cos(w * t), np.sin(w * t)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    a, b = coef[1], coef[2]
    delta = math.atan2(b, a)
    tmax = (delta / w) % period
    return float(math.hypot(a, b)), float(tmax)


def matched_filter_test(t, y, period, t0, n_sur, rng):
    """Correlation of detrended y with the predicted cosine, with two nulls.

    p_ar1: AR(1) surrogates with the series' own lag-1 autocorrelation.
    p_phase: same data, template phase drawn uniformly at random.
    Both one-sided (large positive r = support for the ephemeris).
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    yd = detrend(t, y)
    c = template(t, period, t0)
    r = float(corr_rows(yd, c)[0])
    phi = lag1(yd)
    X = detrend_rows(t, ar1_surrogates(len(t), phi, n_sur, rng))
    rs = corr_rows(X, c)
    p_ar1 = (1 + int(np.sum(rs >= r))) / (n_sur + 1)
    phases = rng.uniform(0, period, n_sur)
    C = np.cos(2 * np.pi * (t[None, :] - phases[:, None]) / period)
    rp = corr_rows(C, yd)
    p_phase = (1 + int(np.sum(rp >= r))) / (n_sur + 1)
    amp, tmax = fit_phase(t, yd, period)
    return {"n": int(len(t)), "r": r, "lag1": phi, "p_ar1": p_ar1,
            "p_phase": p_phase, "fit_amplitude": amp, "fit_tmax_ma": tmax,
            "phase_offset_myr": float(circ_dist(tmax, period, t0))}


def ls_matrices(t, freqs):
    t = np.asarray(t, float)
    w = 2 * np.pi * np.asarray(freqs, float)[:, None]
    tau = np.arctan2(np.sum(np.sin(2 * w * t[None, :]), axis=1),
                     np.sum(np.cos(2 * w * t[None, :]), axis=1)) / (2 * w[:, 0])
    arg = w * (t[None, :] - tau[:, None])
    C = np.cos(arg)
    S = np.sin(arg)
    return C, S, np.sum(C * C, axis=1), np.sum(S * S, axis=1)


def ls_power(Y, mats):
    """Normalised Lomb-Scargle power for each row of Y."""
    C, S, cc, ss = mats
    Y = np.atleast_2d(np.asarray(Y, float))
    Yc = Y - Y.mean(axis=1, keepdims=True)
    var = Yc.var(axis=1, ddof=1)
    return ((Yc @ C.T) ** 2 / cc + (Yc @ S.T) ** 2 / ss) / (2 * var[:, None])


def periodogram_test(t, y, n_sur, rng, pmin=10.0, pmax=100.0, band=(26.0, 28.5)):
    t = np.asarray(t, float)
    yd = detrend(t, np.asarray(y, float))
    freqs = np.arange(1.0 / pmax, 1.0 / pmin + 1e-12, 0.0001)
    pw = ls_power(yd, ls_matrices(t, freqs))[0]
    bfreqs = np.linspace(1.0 / band[1], 1.0 / band[0], 41)
    bm = ls_matrices(t, bfreqs)
    band_obs = float(np.max(ls_power(yd, bm)[0]))
    phi = lag1(yd)
    exceed = 0
    done = 0
    while done < n_sur:
        m = min(2000, n_sur - done)
        X = detrend_rows(t, ar1_surrogates(len(t), phi, m, rng))
        exceed += int(np.sum(np.max(ls_power(X, bm), axis=1) >= band_obs))
        done += m
    # local maxima of the full periodogram, strongest first
    loc = [k for k in range(1, len(pw) - 1) if pw[k] >= pw[k - 1] and pw[k] >= pw[k + 1]]
    loc.sort(key=lambda k: -pw[k])
    top = [{"period_myr": float(1.0 / freqs[k]), "power": float(pw[k])} for k in loc[:8]]
    gk = int(np.argmax(pw))
    in_band = band[0] <= 1.0 / freqs[gk] <= band[1]
    return {"band": list(band), "band_max_power": band_obs,
            "p_band_ar1": (1 + exceed) / (n_sur + 1),
            "global_peak_period_myr": float(1.0 / freqs[gk]),
            "global_peak_in_band": bool(in_band), "top_peaks": top,
            "freqs": freqs, "power": pw}


def peaks(y, q=0.75):
    """Indices of local maxima (series order) at or above the q-quantile."""
    y = np.asarray(y, float)
    thr = float(np.quantile(y, q))
    return [k for k in range(1, len(y) - 1)
            if y[k] > y[k - 1] and y[k] > y[k + 1] and y[k] >= thr]


def binom_sf(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p)."""
    if k <= 0:
        return 1.0
    return float(sum(math.comb(n, j) * p ** j * (1 - p) ** (n - j) for j in range(k, n + 1)))


def event_test(t, y, period, t0, half=HIT_HALF_WIDTH):
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    yd = detrend(t, y)
    idx = peaks(yd)
    d = circ_dist(t[idx], period, t0) if idx else np.array([])
    hits = int(np.sum(d <= half))
    p0 = 2 * half / period
    return {"n_peaks": len(idx), "hits": hits, "p0": p0,
            "p_binom": binom_sf(hits, len(idx), p0) if idx else 1.0,
            "peak_ages_ma": [float(x) for x in t[idx]],
            "peak_offsets_myr": [float(x) for x in d]}


def normal_scores(y):
    y = np.asarray(y, float)
    ranks = np.argsort(np.argsort(y, kind="mergesort"), kind="mergesort")
    nd = NormalDist()
    return np.array([nd.inv_cdf((r + 0.5) / len(y)) for r in ranks])


def jitter_r(t, y, period, t0, n_draws, rng):
    """Correlation r when every bin top age is perturbed by its timescale error."""
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    sig = np.where(t < 66.0, 0.25, np.where(t < 252.0, 0.75, 1.0))
    out = np.empty(n_draws)
    for k in range(n_draws):
        tj = t + rng.standard_normal(len(t)) * sig
        out[k] = corr_rows(detrend(tj, y), template(tj, period, t0))[0]
    return {"r_p05": float(np.quantile(out, 0.05)), "r_p50": float(np.median(out)),
            "r_p95": float(np.quantile(out, 0.95)), "frac_r_pos": float(np.mean(out > 0))}


def redated_hits(period=27.16, t0=9.5, half=HIT_HALF_WIDTH):
    rows = []
    for name, a12, a20 in MB2014_EVENTS:
        d12 = float(circ_dist(a12, period, t0))
        d20 = float(circ_dist(a20, period, t0))
        rows.append({"event": name, "age_gts2012": a12, "age_gts2020": a20,
                     "offset_2012": d12, "offset_2020": d20,
                     "hit_2012": d12 <= half, "hit_2020": d20 <= half})
    return rows
