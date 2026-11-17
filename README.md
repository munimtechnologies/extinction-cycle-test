# extinction-cycle-test

**Is there a clock in the history of life?** This repository holds open,
pre-registered, out-of-sample tests of two famous claims:

1. Marine mass extinctions recur every ~27 Myr (Raup & Sepkoski 1984; Melott &
   Bambach 2010-2014). Plan: [PREREG.md](PREREG.md) (PREREG-1).
2. Marine genus diversity rises and falls every ~62 Myr (Rohde & Muller 2005,
   Nature; Melott & Bambach 2011, 2014). Plan: [PREREG-2.md](PREREG-2.md).

Status (2026-11-17): **both pre-registered; results pending.** Each plan, its
decision rule and our probability for each outcome were committed before any data
were analysed. Results will be added here when the jobs finish. They will be
published whichever way they come out.

## Why now

- Melott & Bambach published fixed ephemerides: extinction maxima at
  **9.5 + 27.16 n Ma**, and diversity maxima at **27.8 + 62 n Ma**. A published
  ephemeris tested on data it never saw is a real prediction.
- Since then, the Paleobiology Database has grown and the 2020 Geologic Time
  Scale has moved many boundaries. The 2014 PBDB series used ~11-Myr bins, which
  barely resolve a 27-Myr period. We use single stages (~4-5 Myr).
- Fan et al. (Science 2020) published an independent Paleozoic record at
  ~26-kyr resolution (CC0). We found no published periodicity test on it.

## What is in this repository

| File | What it does |
|---|---|
| `PREREG.md` | 27-Myr test: hypotheses, tests, sanity gate, our priors (P1-P9), part C for Fan et al. |
| `PREREG-2.md` | 62-Myr test: hypotheses, sanity gate, our priors (Q1-Q9) |
| `cycle.py` | GTS2020 stage table, binning, Alroy three-timer and Foote rates, matched-filter test of a fixed ephemeris with AR(1) and random-phase nulls, Lomb-Scargle, event test, timescale jitter |
| `cycle62.py` | Range-through richness, cubic-detrended matched filter and periodogram, declining-phase contrast |
| `job.py`, `job62.py` | Download PBDB marine invertebrate genus occurrences and run each test end to end |
| `test_cycle.py`, `test_cycle62.py` | Offline unit tests (`python -m unittest -v test_cycle test_cycle62`) |

Run it yourself: `pip install numpy && python -m unittest -v test_cycle test_cycle62 && python job.py --data-dir data --out-dir out && python job62.py --data-dir data --out-dir out`.

## Data and credit

- The Paleobiology Database (paleobiodb.org), data1.2 API. Please cite the PBDB
  and its contributors.
- Fan, J.-x., Shen, S.-z. et al. 2020, *Science* 367:272, data CC0 at
  https://zenodo.org/records/4977790.
- Melott, A. L. & Bambach, R. K. 2014, *Paleobiology* 40:177-196 (both
  ephemerides and the 19-event table re-dated in parts A and A-ii).
- Rohde, R. A. & Muller, R. A. 2005, *Nature* 434:208-210.

## Disclosure and licence

Made with AI agents working in MV2's name. Every number will be produced by the
code in this repository from public data. Code: Apache-2.0. Text: CC BY 4.0.

MV2, Munim, Inc., mv2@munimtech.com
