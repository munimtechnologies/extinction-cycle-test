# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""PREREG-2: out-of-sample test of the ~62-Myr cycle in marine genus diversity
(Rohde & Muller 2005; Melott & Bambach 2011, 2014) on the 2026 PBDB at
GTS2020 stage resolution. Pure functions; job62.py does the I/O.
Made with AI agents working in MV2's name.
"""
import numpy as np

import cycle

# Melott & Bambach 2014 (Paleobiology 40:177, arXiv:1310.4712): "The latest
# peak in the 62-Myr periodicity occurred 27.8 Myr ago." Period 62 Myr
# (Rohde & Muller 2005: 62 +/- 3 Myr). Maxima of DIVERSITY at 27.8 + 62 n Ma.
P62 = 62.0
T0_62 = 27.8
WINDOW62 = (5.0, 520.0)        # stage MIDPOINTS used (Melott & Bambach 2011: 5-520 Ma)
SHORT_LIVED_MAX = 45.0         # Myr, Rohde & Muller's short-lived genera
BAND62 = (58.0, 66.0)


# ------------------------------------------------------------ diversity
def range_through(sets, mids=None, max_dur=None):
    """Range-through genus richness per bin (bins oldest first).

    A genus counts in every bin from its first to its last bin. If max_dur
    is given, only genera with mids[first] - mids[last] <= max_dur count.
    """
    n = len(sets)
    first, last = {}, {}
    for i, s in enumerate(sets):
        for k in s:
            if k not in first:
                first[k] = i
            last[k] = i
    out = [0] * n
    for k, f in first.items():
        l = last[k]
        if max_dur is not None and (mids[f] - mids[l]) > max_dur:
            continue
        for i in range(f, l + 1):
            out[i] += 1
    return out


# ------------------------------------------------------------ statistics
def poly_basis(t, deg):
    t = np.asarray(t, float)
    sd = float(t.std()) or 1.0
    z = (t - t.mean()) / sd
    return np.vstack([z ** k for k in range(deg + 1)]).T


def detrend_poly(t, y, deg=3):
    A = poly_basis(t, deg)
    coef, *_ = np.linalg.lstsq(A, np.asarray(y, float), rcond=None)
    return np.asarray(y, float) - A @ coef


def detrend_rows_poly(t, Y, deg=3):
    A = poly_basis(t, deg)
    P = A @ np.linalg.pinv(A)
    return np.asarray(Y, float) - np.asarray(Y, float) @ P


def matched_filter_poly(t, y, period, t0, n_sur, rng, deg=3):
    """Like cycle.matched_filter_test, with a polynomial (default cubic) detrend."""
    t = np.asarray(t, float)
    yd = detrend_poly(t, y, deg)
    c = cycle.template(t, period, t0)
    r = float(cycle.corr_rows(yd, c)[0])
    phi = cycle.lag1(yd)
    X = detrend_rows_poly(t, cycle.ar1_surrogates(len(t), phi, n_sur, rng), deg)
    rs = cycle.corr_rows(X, c)
    p_ar1 = (1 + int(np.sum(rs >= r))) / (n_sur + 1)
    phases = rng.uniform(0, period, n_sur)
    C = np.cos(2 * np.pi * (t[None, :] - phases[:, None]) / period)
    rp = cycle.corr_rows(C, yd)
    p_phase = (1 + int(np.sum(rp >= r))) / (n_sur + 1)
    amp, tmax = cycle.fit_phase(t, yd, period)
    return {"n": int(len(t)), "r": r, "lag1": phi, "p_ar1": p_ar1, "p_phase": p_phase,
            "fit_amplitude": amp, "fit_tmax_ma": tmax,
            "phase_offset_myr": float(cycle.circ_dist(tmax, period, t0)), "detrend_deg": deg}


def periodogram_poly(t, y, n_sur, rng, pmin=20.0, pmax=200.0, band=BAND62, deg=3):
    t = np.asarray(t, float)
    yd = detrend_poly(t, y, deg)
    freqs = np.arange(1.0 / pmax, 1.0 / pmin + 1e-12, 0.0001)
    pw = cycle.ls_power(yd, cycle.ls_matrices(t, freqs))[0]
    bm = cycle.ls_matrices(t, np.linspace(1.0 / band[1], 1.0 / band[0], 41))
    band_obs = float(np.max(cycle.ls_power(yd, bm)[0]))
    phi = cycle.lag1(yd)
    exceed, done = 0, 0
    while done < n_sur:
        m = min(2000, n_sur - done)
        X = detrend_rows_poly(t, cycle.ar1_surrogates(len(t), phi, m, rng), deg)
        exceed += int(np.sum(np.max(cycle.ls_power(X, bm), axis=1) >= band_obs))
        done += m
    loc = [k for k in range(1, len(pw) - 1) if pw[k] >= pw[k - 1] and pw[k] >= pw[k + 1]]
    loc.sort(key=lambda k: -pw[k])
    gk = int(np.argmax(pw))
    return {"band": list(band), "band_max_power": band_obs,
            "p_band_ar1": (1 + exceed) / (n_sur + 1),
            "global_peak_period_myr": float(1.0 / freqs[gk]),
            "global_peak_in_band": bool(band[0] <= 1.0 / freqs[gk] <= band[1]),
            "top_peaks": [{"period_myr": float(1.0 / freqs[k]), "power": float(pw[k])}
                          for k in loc[:8]],
            "freqs": freqs, "power": pw}


def jitter_r_poly(t, y, period, t0, n_draws, rng, deg=3):
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    sig = np.where(t < 66.0, 0.25, np.where(t < 252.0, 0.75, 1.0))
    out = np.empty(n_draws)
    for k in range(n_draws):
        tj = t + rng.standard_normal(len(t)) * sig
        out[k] = cycle.corr_rows(detrend_poly(tj, y, deg), cycle.template(tj, period, t0))[0]
    return {"r_p05": float(np.quantile(out, 0.05)), "r_p50": float(np.median(out)),
            "r_p95": float(np.quantile(out, 0.95)), "frac_r_pos": float(np.mean(out > 0))}


# ------------------------------------------------- declining-phase tests
def declining(t, period=P62, t0=T0_62):
    """True where age t lies in a declining-diversity half-cycle, i.e. within
    period/2 Myr AFTER (younger than) a predicted maximum at t0 + n*period."""
    u = np.mod(t0 - np.asarray(t, float), period)
    return u < period / 2.0


def phase_contrast(t, y, n_sur, rng, period=P62, t0=T0_62):
    """Mean normal score of linearly detrended y in declining minus increasing
    half-cycles; one-sided p against templates with uniformly random t0."""
    t = np.asarray(t, float)
    z = cycle.normal_scores(cycle.detrend(t, np.asarray(y, float)))

    def stat(D):
        D = np.atleast_2d(D)
        nd = D.sum(axis=1)
        ni = (~D).sum(axis=1)
        sd = np.where(nd > 0, (D @ z) / np.maximum(nd, 1), 0.0)
        si = np.where(ni > 0, ((~D) @ z) / np.maximum(ni, 1), 0.0)
        return np.where((nd > 0) & (ni > 0), sd - si, 0.0)

    D0 = declining(t, period, t0)
    obs = float(stat(D0)[0])
    phases = rng.uniform(0, period, n_sur)
    Dn = np.mod(phases[:, None] - t[None, :], period) < period / 2.0
    null = stat(Dn)
    return {"n": int(len(t)), "n_declining": int(D0.sum()), "diff_normal_scores": obs,
            "p_random_phase": (1 + int(np.sum(null >= obs))) / (n_sur + 1)}


def mb_events_declining():
    rows = []
    for name, a12, a20 in cycle.MB2014_EVENTS:
        rows.append({"event": name, "age_gts2012": a12, "age_gts2020": a20,
                     "declining_2012": bool(declining(a12)),
                     "declining_2020": bool(declining(a20))})
    return rows
