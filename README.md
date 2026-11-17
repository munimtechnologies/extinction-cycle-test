# extinction-cycle-test

**Is there a 27-million-year clock in the history of life?** This is an open,
pre-registered, out-of-sample test of the claim that marine mass extinctions
recur every ~27 Myr (Raup & Sepkoski 1984; Melott & Bambach 2010-2014).

Status (2026-11-17): **pre-registered; results pending.** The plan, the decision
rule and our probability for each outcome are fixed in [PREREG.md](PREREG.md),
committed before any data was downloaded. The results will be added here when the
job finishes. They will be published whichever way they come out.

## Why now

- Melott & Bambach (Paleobiology 2014) published a fixed ephemeris: extinction
  maxima at **9.5 + 27.16 n Ma**. A published ephemeris tested on data it
  never saw is a real prediction.
- Since then, the Paleobiology Database has grown and the 2020 Geologic Time
  Scale has moved many boundaries. The 2014 PBDB series used ~11-Myr bins, which
  barely resolve a 27-Myr period. We use single stages (~4-5 Myr).
- Fan et al. (Science 2020) published an independent Paleozoic record at
  ~26-kyr resolution (CC0). We found no published periodicity test on it.

## What is in this repository

| File | What it does |
|---|---|
| `PREREG.md` | Hypotheses, tests, sanity gate, our priors (P1-P9), and part C for Fan et al. |
| `cycle.py` | GTS2020 stage table, binning, Alroy three-timer and Foote rates, matched-filter test of a fixed ephemeris with AR(1) and random-phase nulls, Lomb-Scargle, event test, timescale jitter |
| `job.py` | Downloads PBDB marine invertebrate genus occurrences and the Fan et al. archive, runs everything, writes `results.json`, `series.csv`, `periodogram_mu3t.csv`, `provenance.json`, `fan_inventory.txt` |
| `test_cycle.py` | Offline unit tests (`python -m unittest -v test_cycle`) |

Run it yourself: `pip install numpy && python -m unittest -v test_cycle && python job.py --data-dir data --out-dir out`.

## Data and credit

- The Paleobiology Database (paleobiodb.org), data1.2 API. Please cite the PBDB
  and its contributors.
- Fan, J.-x., Shen, S.-z. et al. 2020, *Science* 367:272, data CC0 at
  https://zenodo.org/records/4977790.
- Melott, A. L. & Bambach, R. K. 2014, *Paleobiology* 40:177-196 (the
  ephemeris and the 19-event table re-dated in part A).

## Disclosure and licence

Made with AI agents working in MV2's name. Every number will be produced by the
code in this repository from public data. Code: Apache-2.0. Text: CC BY 4.0.

MV2, Munim, Inc., mv2@munimtech.com
