# extinction-cycle-test

**Is there a clock in the history of life?** This repository holds open,
pre-registered, out-of-sample tests of three famous claims:

1. Marine mass extinctions recur every ~27 Myr (Raup & Sepkoski 1984; Melott &
   Bambach 2010-2014). Plan: [PREREG.md](PREREG.md) (PREREG-1).
2. Marine genus diversity rises and falls every ~62 Myr (Rohde & Muller 2005,
   Nature; Melott & Bambach 2011, 2014). Plan: [PREREG-2.md](PREREG-2.md).
3. Extinctions of land animals (tetrapods) recur every ~27.5 Myr, with the
   most recent pulse ~7 Ma (Rampino, Caldeira & Zhu 2021, Historical Biology;
   Geoscience Frontiers). Plan: [PREREG-3.md](PREREG-3.md).

Status (2026-11-17): **all three tests are done.** Each plan, its decision rule and our
probability for each outcome were committed before any data were analysed, and each grader
was committed before its job produced output. Every result is published whichever way it came
out. **[SYNTHESIS.md](SYNTHESIS.md)** reads these three together with two sister tests on
temperature: *take away the Big Five mass extinctions and the famous signals in extinction go
quiet; the 62-Myr diversity rhythm survives, but it is in the rock record too.*

**Result 1 ([RESULTS-1.md](RESULTS-1.md)): RESTS_ON_MASS_EXTINCTIONS.** The published
27.16-Myr ephemeris lines up with the 2026 PBDB marine extinction record at stage resolution
(r = 0.24, p = 0.021, fitted phase 9.51 Ma), but only because of the Big Five mass extinctions
it was fitted to. Without those six stages: r = 0.12, p = 0.17. The record's strongest
periodicity is ~38 Myr. Part C (Fan et al. 2020) could not be run: the archive holds the
inputs to their optimisation, not the dated species ranges.

**Result 2 ([RESULTS-2.md](RESULTS-2.md)): SUPPORTED; headline
PREDICTS_DIVERSITY_AND_SAMPLING.** The 2014 ephemeris (diversity maxima at 27.8 + 62 n Ma)
predicts 2026 PBDB range-through genus richness on GTS2020 over 91 stages: r = 0.58,
p < 0.0001 (AR(1) null); the strongest period is 60.6 Myr and the fitted phase is 2.6 Myr
from the published one. All nine registered predictions held. But the number of fossil
collections per stage follows the same rhythm (r = 0.31, p = 0.015), so this test cannot tell
a rhythm in life from a rhythm in the rock record.

**Result 3 ([RESULTS-3.md](RESULTS-3.md)): DATA_LIMITED.** Only 38 of 57 stages (gate: 40)
have a defined three-timer extinction rate for land tetrapods, so the registered test makes no
claim. For the record, no version of the series shows the cycle: r = 0.08 (p = 0.31) on the
primary series, and r = 0.007 (p = 0.48) on Foote's rate, which is defined in all 57 stages.

## Why now

- Melott & Bambach published fixed ephemerides: extinction maxima at
  **9.5 + 27.16 n Ma**, and diversity maxima at **27.8 + 62 n Ma**. Rampino et
  al. published pulses ~27.5 Myr apart, the latest ~7 Ma. A published
  ephemeris tested on data it never saw is a real prediction.
- Since then, the Paleobiology Database has grown and the 2020 Geologic Time
  Scale has moved many boundaries. The 2014 PBDB series used ~11-Myr bins, which
  barely resolve a 27-Myr period. We use single stages (~4-5 Myr).
- The land-tetrapod claim rests on a list of ten event ages. We found no test of
  it on an occurrence-based extinction-rate series, and none that leaves out the
  events the cycle was fitted to (searches rechecked 2026-11-17: the published
  critique, Omerbashich 2021, and the reply both work on the event ages only).
- Fan et al. (Science 2020) published an independent Paleozoic record at
  ~26-kyr resolution (CC0). We found no published periodicity test on it.

## What is in this repository

| File | What it does |
|---|---|
| `SYNTHESIS.md` | The five tests (these three plus the two temperature tests) read together |
| `PREREG.md` | 27-Myr test: hypotheses, tests, sanity gate, our priors (P1-P9), part C for Fan et al. |
| `PRE-RESULT-NOTE-1.md`, `grade_cycle1.py` | The mechanical reading of test 1, published before its result |
| `RESULTS-1.md`, `results/prereg1/` | Test 1: graded result, per-stage series, periodogram, provenance |
| `PREREG-2.md` | 62-Myr test: hypotheses, sanity gate, our priors (Q1-Q9) |
| `PRE-RESULT-NOTE-2.md`, `grade_cycle2.py` | The mechanical reading of test 2, published before its result (discloses what test 1's public series already reveal) |
| `RESULTS-2.md`, `results/prereg2/` | Test 2: graded result, per-stage richness and collections, periodogram, provenance |
| `PREREG-3.md`, `PREREG-3-AMENDMENT-1.md` | Land-tetrapod 27.5-Myr test: hypotheses, leave-events-out test, our priors (T1-T11); the corrected PBDB environment code |
| `PRE-RESULT-NOTE-3.md`, `grade_tetra.py` | The mechanical reading of test 3, published before its result |
| `RESULTS-3.md`, `results/prereg3/` | Test 3: graded result, per-stage counts and rates, periodogram, provenance |
| `cycle.py` | GTS2020 stage table, binning, Alroy three-timer and Foote rates, matched-filter test of a fixed ephemeris with AR(1) and random-phase nulls, Lomb-Scargle, event test, timescale jitter |
| `cycle62.py` | Range-through richness, cubic-detrended matched filter and periodogram, declining-phase contrast |
| `tetra.py` | Rampino et al. ephemeris, pulses and events, event-stage mapping, the PREREG-3 analysis |
| `job.py`, `job62.py` | Download PBDB marine invertebrate genus occurrences and run tests 1 and 2 end to end |
| `job_tetra.py` | Downloads PBDB terrestrial tetrapod genus occurrences and runs test 3 end to end |
| `test_cycle.py`, `test_cycle62.py`, `test_tetra.py`, `test_grade_cycle1.py`, `test_grade_cycle2.py`, `test_grade_tetra.py` | Offline unit tests (`python -m unittest -v test_cycle test_cycle62 test_tetra test_grade_cycle1 test_grade_cycle2 test_grade_tetra`) |
| `test_graded_results.py` | Runs `grade_cycle2.py` and `grade_tetra.py` on the published outputs in `results/` and checks they give the labels in RESULTS-2 and RESULTS-3 |

Run it yourself: `pip install numpy && python -m unittest -v test_cycle test_cycle62 test_tetra && python job.py --data-dir data --out-dir out && python job62.py --data-dir data --out-dir out && python job_tetra.py --data-dir data --out-dir out`.

## Sister projects

- [extinction-temperature-test](https://github.com/munimtechnologies/extinction-temperature-test):
  does temperature change predict extinction (Song et al. 2021, the 5.2 C threshold)?
  Result: only through the mass extinctions.
- [diversity-temperature-test](https://github.com/munimtechnologies/diversity-temperature-test):
  is a warmer Earth richer in marine life (Mayhew et al. 2012)? Result: not supported.

## Data and credit

- The Paleobiology Database (paleobiodb.org), data1.2 API. Please cite the PBDB
  and its contributors.
- Fan, J.-x., Shen, S.-z. et al. 2020, *Science* 367:272, data CC0 at
  https://zenodo.org/records/4977790.
- Melott, A. L. & Bambach, R. K. 2014, *Paleobiology* 40:177-196 (both
  ephemerides and the 19-event table re-dated in parts A and A-ii).
- Rohde, R. A. & Muller, R. A. 2005, *Nature* 434:208-210.
- Rampino, M. R., Caldeira, K. & Zhu, Y. 2021, *Historical Biology* 33:3084-3090,
  and *Geoscience Frontiers* 12:101245.

## Disclosure and licence

Made with AI agents working in MV2's name. Every number is produced by the
code in this repository from public data. Code: Apache-2.0. Text: CC BY 4.0.

MV2, Munim, Inc., mv2@munimtech.com
