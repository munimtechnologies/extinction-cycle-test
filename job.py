# Copyright 2026 MV2, Munim, Inc. (mv2@munimtech.com)
# SPDX-License-Identifier: Apache-2.0
"""Run the pre-registered extinction-cycle test (PREREG.md) end to end.

    python job.py --data-dir /data --out-dir /out

Downloads PBDB marine invertebrate genus occurrences (cached in --data-dir),
builds GTS2020 stage-level extinction-rate series and runs the tests fixed in
PREREG.md. Also downloads Fan et al. (2020)'s CC0 dataset and writes an
INVENTORY ONLY (file list and a few header lines; no statistics), so its
parser can be written before any of its numbers are looked at.
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
import zipfile

import numpy as np

import cycle

PBDB = "https://paleobiodb.org/data1.2/occs/list.csv"
GROUPS = ["Porifera", "Cnidaria", "Bryozoa", "Brachiopoda", "Mollusca",
          "Echinodermata", "Trilobita", "Graptolithina"]
REQUIRED = {"Brachiopoda", "Mollusca", "Echinodermata", "Trilobita"}
COMMON = "envtype=marine&idreso=genus&idqual=genus_certain&pres=regular&show=class"
FAN_URL = ("https://zenodo.org/records/4977790/files/"
           "Dataset-PaleozoicMarineDiversity.zip?download=1")
FAN_MD5 = "8ee5ccd5b8bdea48ae8cd2310ea36b93"
UA = "MV2-extinction-cycle/0.1 (Munim, Inc.; mv2@munimtech.com)"


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
        except Exception as e:  # noqa: BLE001 - report and retry
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


def fan_inventory(data_dir, out_dir, prov):
    d = os.path.join(data_dir, "fan2020")
    os.makedirs(d, exist_ok=True)
    dest = os.path.join(d, "Dataset-PaleozoicMarineDiversity.zip")
    lines = []
    try:
        prov["status"] = fetch(FAN_URL, dest)
        prov["md5"] = md5(dest)
        prov["md5_ok"] = prov["md5"] == FAN_MD5
        with zipfile.ZipFile(dest) as z:
            infos = z.infolist()
            prov["n_entries"] = len(infos)
            lines.append(f"# {len(infos)} entries; md5 {prov['md5']} (expected {FAN_MD5})")
            for zi in infos:
                lines.append(f"{zi.file_size}\t{zi.filename}")
            files = [zi for zi in infos if not zi.is_dir()]
            keys = ("read", "fad", "lad", "range", "compos", "sequence", "species",
                    "divers", "result", "output", "time", "age", "sol", "best")
            chosen = [zi for zi in files if "/" not in zi.filename.strip("/")
                      or any(k in zi.filename.lower() for k in keys)]
            chosen += sorted(files, key=lambda zi: -zi.file_size)[:10]
            seen = set()
            lines.append("\n# HEADS (first 6 lines, 300 chars each) -- format only, no statistics")
            for zi in chosen:
                if zi.filename in seen or len(seen) >= 40:
                    continue
                seen.add(zi.filename)
                with z.open(zi) as f:
                    raw = f.read(20000)
                text = raw.decode("utf-8", errors="replace").splitlines()[:6]
                lines.append(f"\n## {zi.filename} ({zi.file_size} bytes)")
                lines.extend(x[:300] for x in text)
    except Exception as e:  # noqa: BLE001
        prov["error"] = repr(e)[:500]
        lines.append("ERROR " + prov["error"])
    with open(os.path.join(out_dir, "fan_inventory.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def window(names, tops, metric):
    lo, hi = cycle.WINDOW
    keep = [i for i in range(len(tops))
            if lo <= tops[i] <= hi and metric[i] == metric[i]]
    keep.sort(key=lambda i: tops[i])
    return ([names[i] for i in keep], np.array([tops[i] for i in keep]),
            np.array([metric[i] for i in keep]))


def run_tests(t, y, n_sur, rng, with_periodogram=True):
    res = {}
    for key, e in cycle.EPHEMERIDES.items():
        res[key] = {"matched_filter": cycle.matched_filter_test(t, y, e["period"], e["t0"], n_sur, rng),
                    "events": cycle.event_test(t, y, e["period"], e["t0"])}
    if with_periodogram:
        pg = cycle.periodogram_test(t, y, n_sur, rng)
        res["periodogram"] = {k: v for k, v in pg.items() if k not in ("freqs", "power")}
        res["_pg_arrays"] = (pg["freqs"], pg["power"])
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--n-sur", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20261117)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    prov = {"started_utc": now(), "pbdb": {}, "fan2020": {}, "groups": GROUPS,
            "query_common": COMMON, "n_surrogates": args.n_sur, "seed": args.seed}

    paths = get_pbdb(args.data_dir, prov["pbdb"])
    fan_inventory(args.data_dir, args.out_dir, prov["fan2020"])
    results = {"prereg": "PREREG.md (extinction-cycle-test v0.1.0)", "status": "ok"}
    got = {g for g, _ in paths}
    if not REQUIRED <= got:
        results["status"] = "data_incomplete: missing " + ",".join(sorted(REQUIRED - got))
        _write(args.out_dir, results, prov)
        return 2

    assigner = cycle.BinAssigner()
    acc = cycle.Accumulator(assigner.n)
    per_group = {}
    for g, p in paths:
        before = dict(acc.stats)
        with open(p, newline="", encoding="utf-8", errors="replace") as f:
            acc.add_rows(csv.DictReader(f), assigner)
        per_group[g] = {k: acc.stats[k] - before[k] for k in acc.stats}
    prov["parse"] = {"total": acc.stats, "per_group": per_group}

    names = [s[0] for s in cycle.STAGES_GTS2020]
    tops = [s[2] for s in cycle.STAGES_GTS2020]
    mu, ps, c3 = cycle.three_timer(acc.sets)
    q, prop, cf = cycle.foote_rates(acc.sets)
    with open(os.path.join(args.out_dir, "series.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["bin", "base_ma", "top_ma", "n_occ", "n_coll", "n_genera", "two_timers",
                    "three_timers", "part_timers", "ps", "mu3t", "Nbt", "NbL", "NFt", "NFL",
                    "foote_q", "prop_ext"])
        for i, s in enumerate(cycle.STAGES_GTS2020):
            w.writerow([s[0], s[1], s[2], acc.n_occ[i], acc.n_coll[i], len(acc.sets[i]),
                        c3["two_timers"][i], c3["three_timers"][i], c3["part_timers"][i],
                        ps[i], mu[i], cf["Nbt"][i], cf["NbL"][i], cf["NFt"][i], cf["NFL"][i],
                        q[i], prop[i]])

    # Sanity gate (PREREG S1): >= 4 of the Big Five among the 10 highest mu3t bins.
    nw, tw, yw = window(names, tops, mu)
    top10 = [nw[k] for k in np.argsort(-yw)[:10]]
    hit5 = sum(any(b in top10 for b in ev) for ev in cycle.BIG_FIVE_EVENTS)
    results["sanity"] = {"top10_mu3t_bins": top10, "big_five_in_top10": hit5,
                         "passed": hit5 >= 4, "n_bins_in_window": len(tw)}

    results["partA_redating"] = cycle.redated_hits()
    results["partA_hits"] = {
        "gts2012": sum(r["hit_2012"] for r in results["partA_redating"]),
        "gts2020": sum(r["hit_2020"] for r in results["partA_redating"]), "n": 19}

    # Primary: mu3t, all ephemerides, periodogram.
    prim = run_tests(tw, yw, args.n_sur, rng)
    fr, pw = prim.pop("_pg_arrays")
    with open(os.path.join(args.out_dir, "periodogram_mu3t.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["freq_per_myr", "period_myr", "power"])
        for a, b in zip(fr, pw):
            w.writerow([f"{a:.5f}", f"{1 / a:.4f}", f"{b:.6f}"])
    e = cycle.EPHEMERIDES["MB2014"]
    prim["jitter_MB2014"] = cycle.jitter_r(tw, yw, e["period"], e["t0"], 1000, rng)
    results["mu3t"] = prim
    results["mu3t_bins"] = [{"bin": a, "top_ma": float(b), "mu3t": float(c)}
                            for a, b, c in zip(nw, tw, yw)]

    # Secondary: ranks (normal scores) and leave-Big-Five-out.
    results["mu3t_normal_scores"] = run_tests(tw, cycle.normal_scores(yw), args.n_sur, rng,
                                              with_periodogram=False)
    keep = [k for k, nm in enumerate(nw) if nm not in cycle.BIG_FIVE_BINS]
    results["mu3t_without_big_five"] = run_tests(tw[keep], yw[keep], args.n_sur, rng,
                                                 with_periodogram=False)
    results["mu3t_without_big_five"]["n_removed"] = len(tw) - len(keep)

    # Secondary metrics.
    for label, metric in (("foote_q", q), ("prop_ext", prop)):
        _, tm, ym = window(names, tops, metric)
        r = run_tests(tm, ym, args.n_sur, rng)
        r.pop("_pg_arrays")
        results[label] = r

    # Decision (PREREG H1).
    mf = results["mu3t"]["MB2014"]["matched_filter"]
    results["decision_H1"] = ("SUPPORTED" if (results["sanity"]["passed"] and mf["r"] > 0
                                              and mf["p_ar1"] < 0.05)
                              else ("SANITY_FAILED" if not results["sanity"]["passed"]
                                    else "NOT_SUPPORTED"))
    prov["finished_utc"] = now()
    _write(args.out_dir, results, prov)
    return 0


def _write(out_dir, results, prov):
    with open(os.path.join(out_dir, "results.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1, default=float)
    with open(os.path.join(out_dir, "provenance.json"), "w", encoding="utf-8") as f:
        json.dump(prov, f, indent=1, default=str)
    print(json.dumps({"status": results.get("status"), "decision_H1": results.get("decision_H1"),
                      "sanity": results.get("sanity")}, default=str), flush=True)


if __name__ == "__main__":
    sys.exit(main())
