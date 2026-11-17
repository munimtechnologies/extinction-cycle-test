# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Run PREREG-2 (the 62-Myr diversity cycle) end to end.

    python job62.py --data-dir /data --out-dir /out

Uses the same PBDB queries and cache layout as job.py (PREREG-1), so both
tests read the same download when they share a data volume.
Made with AI agents working in MV2's name.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import shutil
import sys
import time
import urllib.request

import numpy as np

import cycle
import cycle62 as c6

PBDB = "https://paleobiodb.org/data1.2/occs/list.csv"
GROUPS = ["Porifera", "Cnidaria", "Bryozoa", "Brachiopoda", "Mollusca",
          "Echinodermata", "Trilobita", "Graptolithina"]
REQUIRED = {"Brachiopoda", "Mollusca", "Echinodermata", "Trilobita"}
COMMON = "envtype=marine&idreso=genus&idqual=genus_certain&pres=regular&show=class"
UA = "MV2-extinction-cycle/0.2 (Munim, Inc.; mv2@munimtech.com)"


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url, dest, timeout=1800, tries=3):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return "cached"
    tmp = dest + ".part"
    last = ""
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as f:
                shutil.copyfileobj(r, f, 1 << 20)
            os.replace(tmp, dest)
            return "downloaded"
        except Exception as e:  # noqa: BLE001
            last = repr(e)
            print(f"fetch failed ({k + 1}/{tries}) {url}: {last}", flush=True)
            time.sleep(30 * (k + 1))
    raise RuntimeError(last)


def get_pbdb(data_dir, prov):
    d = os.path.join(data_dir, "pbdb")
    os.makedirs(d, exist_ok=True)
    paths = []
    for g in GROUPS:
        url = f"{PBDB}?base_name={g}&{COMMON}"
        dest = os.path.join(d, f"{g}.csv")
        try:
            how = fetch(url, dest)
            with open(dest, encoding="utf-8", errors="replace") as f:
                head = f.readline()
            if "occurrence_no" not in head:
                os.remove(dest)
                raise RuntimeError("unexpected header: " + head[:300])
            prov[g] = {"url": url, "status": how, "bytes": os.path.getsize(dest),
                       "md5": md5(dest),
                       "file_mtime_utc": dt.datetime.fromtimestamp(
                           os.path.getmtime(dest), dt.timezone.utc).isoformat()}
            paths.append((g, dest))
        except Exception as e:  # noqa: BLE001
            prov[g] = {"url": url, "error": repr(e)[:500]}
        print(g, prov[g].get("status", prov[g].get("error")), flush=True)
    return paths


def sel(xs, vals, lo, hi):
    keep = [i for i in range(len(xs)) if lo <= xs[i] <= hi and vals[i] == vals[i]]
    keep.sort(key=lambda i: xs[i])
    return keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--n-sur", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20261118)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    prov = {"started_utc": now(), "pbdb": {}, "groups": GROUPS, "query_common": COMMON,
            "n_surrogates": args.n_sur, "seed": args.seed}
    res = {"prereg": "PREREG-2.md (extinction-cycle-test)", "status": "ok"}

    paths = get_pbdb(args.data_dir, prov["pbdb"])
    if not REQUIRED <= {g for g, _ in paths}:
        res["status"] = "data_incomplete"
        return _write(args.out_dir, res, prov, 2)

    assigner = cycle.BinAssigner()
    acc = cycle.Accumulator(assigner.n)
    for _, p in paths:
        with open(p, newline="", encoding="utf-8", errors="replace") as f:
            acc.add_rows(csv.DictReader(f), assigner)
    prov["parse"] = acc.stats

    st = cycle.STAGES_GTS2020
    names = [s[0] for s in st]
    mids = [(s[1] + s[2]) / 2.0 for s in st]
    tops = [s[2] for s in st]
    s_all = c6.range_through(acc.sets)
    s_short = c6.range_through(acc.sets, mids, c6.SHORT_LIVED_MAX)
    s_sib = [len(s) for s in acc.sets]
    ncoll = acc.n_coll
    logc = [math.log(c) if c > 0 else float("nan") for c in ncoll]
    mu, _, _ = cycle.three_timer(acc.sets)

    with open(os.path.join(args.out_dir, "series62.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["bin", "base_ma", "top_ma", "mid_ma", "range_through", "range_through_short",
                    "sampled_in_bin", "n_coll", "mu3t", "declining_62"])
        for i, s in enumerate(st):
            w.writerow([s[0], s[1], s[2], mids[i], s_all[i], s_short[i], s_sib[i], ncoll[i],
                        mu[i], bool(c6.declining(mids[i]))])

    # Sanity gate S2: the end-Permian, end-Cretaceous and end-Ordovician drops show.
    ix = {n: i for i, n in enumerate(names)}
    checks = {"induan_lt_changhsingian": s_all[ix["Induan"]] < s_all[ix["Changhsingian"]],
              "danian_lt_maastrichtian": s_all[ix["Danian"]] < s_all[ix["Maastrichtian"]],
              "rhuddanian_lt_katian": s_all[ix["Rhuddanian"]] < s_all[ix["Katian"]]}
    res["sanity"] = {"checks": checks, "passed": all(checks.values())}

    lo, hi = c6.WINDOW62
    k = sel(mids, s_all, lo, hi)
    t = np.array([mids[i] for i in k])
    res["window_bins"] = [names[i] for i in k]

    def series(v):
        kk = sel(mids, v, lo, hi)
        return np.array([mids[i] for i in kk]), np.array([v[i] for i in kk], float)

    # H6 (primary): range-through richness vs the MB2014 62-Myr ephemeris.
    y = np.array([s_all[i] for i in k], float)
    res["H6_range_through"] = c6.matched_filter_poly(t, y, c6.P62, c6.T0_62, args.n_sur, rng)
    pg = c6.periodogram_poly(t, y, args.n_sur, rng)
    fr, pw = pg.pop("freqs"), pg.pop("power")
    res["H6_periodogram"] = pg
    with open(os.path.join(args.out_dir, "periodogram_div.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["freq_per_myr", "period_myr", "power"])
        for a, b in zip(fr, pw):
            w.writerow([f"{a:.5f}", f"{1 / a:.4f}", f"{b:.6f}"])
    res["H6_jitter"] = c6.jitter_r_poly(t, y, c6.P62, c6.T0_62, 1000, rng)
    res["H6_quadratic_detrend"] = c6.matched_filter_poly(t, y, c6.P62, c6.T0_62,
                                                         args.n_sur, rng, deg=2)

    # Secondary series.
    for label, v in (("short_lived", s_short), ("sampled_in_bin", s_sib),
                     ("log_collections", logc)):
        ts, ys = series(v)
        res[label] = c6.matched_filter_poly(ts, ys, c6.P62, c6.T0_62, args.n_sur, rng)
    ts, ys = series(s_short)
    sp = c6.periodogram_poly(ts, ys, args.n_sur, rng)
    sp.pop("freqs"), sp.pop("power")
    res["short_lived_periodogram"] = sp

    # H7: extinction (mu3t, PREREG-1 window, bin tops) higher in declining phases.
    kt = sel(tops, mu, cycle.WINDOW[0], cycle.WINDOW[1])
    res["H7_mu3t_declining"] = c6.phase_contrast(np.array([tops[i] for i in kt]),
                                                 np.array([mu[i] for i in kt]), args.n_sur, rng)

    res["partA2_mb_events"] = c6.mb_events_declining()
    res["partA2_counts"] = {"gts2012": sum(r["declining_2012"] for r in res["partA2_mb_events"]),
                            "gts2020": sum(r["declining_2020"] for r in res["partA2_mb_events"]),
                            "n": 19}

    h6 = res["H6_range_through"]
    res["decision_H6"] = ("SANITY_FAILED" if not res["sanity"]["passed"] else
                          "SUPPORTED" if (h6["r"] > 0 and h6["p_ar1"] < 0.05) else
                          "NOT_SUPPORTED")
    prov["finished_utc"] = now()
    return _write(args.out_dir, res, prov, 0)


def _write(out_dir, res, prov, code):
    with open(os.path.join(out_dir, "results62.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    with open(os.path.join(out_dir, "provenance62.json"), "w", encoding="utf-8") as f:
        json.dump(prov, f, indent=1, default=str)
    print(json.dumps({"status": res.get("status"), "decision_H6": res.get("decision_H6"),
                      "sanity": res.get("sanity")}, default=str), flush=True)
    return code


if __name__ == "__main__":
    sys.exit(main())
