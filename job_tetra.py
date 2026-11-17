# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Run PREREG-3 (27.5-Myr cycle in land-tetrapod extinction) end to end.

    python job_tetra.py --data-dir /data --out-dir /out

Downloads PBDB terrestrial tetrapod genus occurrences (two files: mammals and
all other tetrapods) into <data-dir>/pbdb/, bins them into GTS2020 stages,
computes three-timer extinction rates and runs the registered tests.
Made with AI agents working in MV2's name.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import shutil
import sys
import time
import urllib.request

import numpy as np

import cycle
import tetra

UA = "MV2-extinction-cycle/0.3 (Munim, Inc.; mv2@munimtech.com)"


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


def get_data(data_dir, prov):
    d = os.path.join(data_dir, "pbdb")
    os.makedirs(d, exist_ok=True)
    paths = {}
    for label, url in tetra.QUERIES.items():
        dest = os.path.join(d, f"{label}.csv")
        try:
            how = fetch(url, dest)
            with open(dest, encoding="utf-8", errors="replace") as f:
                head = f.readline()
            if "occurrence_no" not in head:
                os.remove(dest)
                raise RuntimeError("unexpected header: " + head[:300])
            prov[label] = {"url": url, "status": how, "bytes": os.path.getsize(dest),
                           "md5": md5(dest),
                           "file_mtime_utc": dt.datetime.fromtimestamp(
                               os.path.getmtime(dest), dt.timezone.utc).isoformat()}
            paths[label] = dest
        except Exception as e:  # noqa: BLE001
            prov[label] = {"url": url, "error": repr(e)[:500]}
        print(label, prov[label].get("status", prov[label].get("error")), flush=True)
    return paths


def accumulate(files, assigner):
    acc = cycle.Accumulator(assigner.n)
    for p in files:
        with open(p, newline="", encoding="utf-8", errors="replace") as f:
            acc.add_rows(csv.DictReader(f), assigner)
    return acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--n-sur", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20261119)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    prov = {"started_utc": now(), "pbdb": {}, "queries": tetra.QUERIES,
            "n_surrogates": args.n_sur, "seed": args.seed}
    res = {"prereg": "PREREG-3.md (extinction-cycle-test)", "status": "ok",
           "ephemeris": tetra.RCZ2021}

    paths = get_data(args.data_dir, prov["pbdb"])
    if set(paths) != set(tetra.QUERIES):
        res["status"] = "data_incomplete"
        return _write(args.out_dir, res, prov, 2)

    st = cycle.STAGES_GTS2020
    assigner = cycle.BinAssigner()
    acc_all = accumulate([paths["tetra_nonmammal"], paths["tetra_mammalia"]], assigner)
    acc_nm = accumulate([paths["tetra_nonmammal"]], assigner)
    prov["parse_all"] = acc_all.stats
    prov["parse_nonmammal"] = acc_nm.stats

    mu, _, cnt = cycle.three_timer(acc_all.sets)
    mu_nm, _, _ = cycle.three_timer(acc_nm.sets)
    q, prop, _ = cycle.foote_rates(acc_all.sets)
    ev = set(tetra.event_bins(st))
    P, T0 = tetra.RCZ2021["period"], tetra.RCZ2021["t0"]

    with open(os.path.join(args.out_dir, "series_tetra.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["bin", "base_ma", "top_ma", "n_occ", "n_coll", "sampled_in_bin",
                    "two_timers", "three_timers", "part_timers", "mu3t_all",
                    "mu3t_nonmammal", "foote_q_all", "prop_all", "rcz_event_bin",
                    "rcz_template"])
        for i, s in enumerate(st):
            w.writerow([s[0], s[1], s[2], acc_all.n_occ[i], acc_all.n_coll[i],
                        len(acc_all.sets[i]), cnt["two_timers"][i], cnt["three_timers"][i],
                        cnt["part_timers"][i], mu[i], mu_nm[i], q[i], prop[i], i in ev,
                        float(cycle.template(s[2], P, T0))])

    # Part A (no PBDB data): where the published events and pulses sit.
    ages = [a for _, a in tetra.TETRAPOD_EVENTS]
    d = cycle.circ_dist(ages, P, T0)
    res["partA_events"] = [{"event": n, "age_ma": a, "offset_myr": float(x),
                            "within_3": bool(x <= 3.0)}
                           for (n, a), x in zip(tetra.TETRAPOD_EVENTS, d)]
    res["partA_pulse_phase"] = dict(zip(("t0_ma", "R"),
                                        tetra.circular_phase(tetra.PULSES_RCZ2021, P)))
    res["partA_event_phase"] = dict(zip(("t0_ma", "R"), tetra.circular_phase(ages, P)))

    # Primary: all land tetrapods, window 2.5-300 Ma.
    prim, (pf, pp) = tetra.analyze(st, mu, args.n_sur, rng)
    res["primary"] = prim
    if pf is not None:
        with open(os.path.join(args.out_dir, "periodogram_tetra.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["freq_per_myr", "period_myr", "power"])
            for a, b in zip(pf, pp):
                w.writerow([f"{a:.5f}", f"{1 / a:.4f}", f"{b:.6f}"])
    # Secondary series.
    res["nonmammal"], _ = tetra.analyze(st, mu_nm, args.n_sur, rng, full=False)
    res["window_260"], _ = tetra.analyze(st, mu, args.n_sur, rng,
                                         window=tetra.WINDOW_260, full=False)
    res["foote_q"], _ = tetra.analyze(st, q, args.n_sur, rng, full=False)
    res["decision"] = prim["decision"]
    prov["finished_utc"] = now()
    return _write(args.out_dir, res, prov, 0)


def _write(out_dir, res, prov, code):
    with open(os.path.join(out_dir, "results_tetra.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    with open(os.path.join(out_dir, "provenance_tetra.json"), "w", encoding="utf-8") as f:
        json.dump(prov, f, indent=1, default=str)
    print(json.dumps({"status": res.get("status"), "decision": res.get("decision"),
                      "sanity": res.get("primary", {}).get("sanity")}, default=str), flush=True)
    return code


if __name__ == "__main__":
    sys.exit(main())
