# RESULTS-2: the 62-Myr marine diversity cycle, out of sample

Graded 2026-11-17 from the output of job `extinction-cycle-62-v1` (ran 20:06-20:10 UTC).
Plan: [PREREG-2.md](PREREG-2.md) (committed 08:33 UTC). Reading:
[grade_cycle2.py](grade_cycle2.py) and [PRE-RESULT-NOTE-2.md](PRE-RESULT-NOTE-2.md)
(committed 19:17 UTC, before the job started). The grader's rules were applied to
`results62.json` by hand, exactly as coded; `test_graded_results.py` runs the grader itself
on the published output and checks that it gives the same labels.
Made with AI agents working in MV2's name.

## Verdict

**SUPPORTED. Registered headline: PREDICTS_DIVERSITY_AND_SAMPLING.**

The published ephemeris (Melott & Bambach 2014: diversity maxima at 27.8 + 62 n Ma, fitted to
older compendia and timescales) predicts the 2026 Paleobiology Database marine genus record on
the 2020 Geologic Time Scale at single-stage resolution, and does so strongly: r = 0.58 over 91
stages (5-520 Ma), p < 0.0001 against 20,000 AR(1) surrogates. The best-fitting phase is
30.4 Ma, 2.6 Myr from the published 27.8 Ma. The strongest period between 20 and 200 Myr is
60.6 Myr. All nine registered predictions came true, including the five we gave 40% or less.

**But the same rhythm is in the sampling record.** The log number of PBDB collections per stage
follows the same ephemeris (r = 0.31, p_AR1 = 0.015, fitted maximum 27.1 Ma). As registered, we
therefore cannot say whether the 62-Myr rhythm is a rhythm in life or a rhythm in how much
fossil-bearing rock was laid down and sampled (or a common cause of both, such as sea level).
This matches earlier work linking the cycle to the sedimentary record (Peters 2008; Melott et
al. 2012, J. Geol. 120:217; Meyers & Peters 2011).

## Graded predictions

| # | Prediction | Prior | Holds? |
|---|---|---|---|
| Q1 | Sanity gate S2 passes (richness drops across the end-Permian, end-Cretaceous and end-Ordovician) | 85% | YES |
| Q2 | H6 r > 0 (range-through richness vs cos(2 pi (t - 27.8)/62)) | 70% | YES |
| Q3 | H6 SUPPORTED: r > 0 and p_AR1 < 0.05 | 40% | YES |
| Q4 | The fitted 62-Myr phase is within +/-8 Myr of 27.8 Ma | 45% | YES |
| Q5 | The global periodogram peak (20-200 Myr) lies in 58-66 Myr | 30% | YES |
| Q6 | The 58-66 Myr band p_AR1 < 0.05 | 30% | YES |
| Q7 | Short-lived genera: r > 0 and p_AR1 < 0.05 | 40% | YES |
| Q8 | log(collections): r > 0 and p_AR1 < 0.05 (the cycle is also in sampling) | 30% | YES |
| Q9 | H7 (extinction higher in declining half-cycles) p_random_phase < 0.05 | 40% | YES |

Brier score of our registered priors: 0.329 over 9 rows (a 50% forecaster scores 0.250). We
were too sceptical: the cycle is far stronger in today's data than we expected.

## Key numbers

| Series (91 stages, cubic detrend unless noted) | r | p_AR1 | p_phase | fitted max (Ma) |
|---|---|---|---|---|
| Range-through richness (H6) | 0.582 | 0.00005 | 0.085 | 30.40 |
| Range-through richness, quadratic detrend | 0.415 | 0.0012 | 0.078 | 30.24 |
| Short-lived genera (<= 45 Myr) | 0.582 | 0.0001 | 0.019 | 28.38 |
| Sampled-in-bin richness | 0.468 | 0.0002 | 0.103 | 31.00 |
| log(collections) | 0.308 | 0.015 | 0.022 | 27.12 |

- p = (k + 1)/20,001, where k is the number of the 20,000 surrogates at least as large:
  0.00005 means none reached the observed r, 0.0001 means one did.
- Periodogram (richness): global peak 60.6 Myr, band-max p_AR1 0.00005; the next peaks are far
  weaker (158.7 Myr, power 4.2 vs 18.4; 24.3 Myr, 3.9). Short-lived genera: peak 61.7 Myr,
  band p_AR1 0.001.
- Timescale jitter (stage ages perturbed): r from 0.566 to 0.594 (5-95%); r > 0 in every draw.
- H7: 50 of 83 extinction-rate stages fall in declining half-cycles; their mean normal score is
  0.52 higher; p = 0.034 against random phases.
- Part A-ii (not a test; disclosed in PREREG-2): 16 of Melott & Bambach's 19 events fall in
  declining phases on GTS2020 (15 on GTS2012).

**p_phase is the honest caveat.** It asks whether the *published* phase fits better than a
random phase of the same period. For total richness it is 0.085: the 62-Myr period is very
strong, but its fitted phase is nearly as good anywhere within a few Myr, so the published
phase itself is not singled out at p < 0.05 for total richness (it is for short-lived genera,
0.019, and collections, 0.022). Q4 (phase within 8 Myr) was the registered phase test and holds.

## What this means, and its limits

- This is a **confirmation, not a discovery**: the 62-Myr cycle was reported in Sepkoski's
  compendium (Rohde & Muller 2005) and in the PBDB at ~11-Myr resolution (Melott 2008; Melott &
  Bambach 2011). What is new is that a phase fixed in 2014 predicts a 2026 database, on a
  timescale it was never fitted to, at about twice the resolution (91 stages averaging
  5.7 Myr, against ~11-Myr bins), under a plan published
  before the data were analysed.
- Related recent work (found in a prior-art recheck while grading, 2026-11-17): Boulila et
  al. (2023, PNAS 120:e2221149120) analysed PBDB diversity and Macrostrat sediment records
  together and reported robust ~28 and ~36 Myr cycles, noting the 62-Myr cycle in both
  biodiversity and rocks; Boulila, Peters et al. (2025, Communications Earth & Environment 6,
  doi:10.1038/s43247-025-02370-6) propose a tectonically driven ~60-Myr redox cycle that
  paces marine biodiversity. Our sampling result is consistent with a cycle shared by rocks
  and life. Neither paper tests a fixed, previously published phase out of sample.
- The parallel signal in collection counts means this test cannot separate biology from the
  rock record. A test that could: the same ephemeris on a sampling-standardised richness
  series with collections partialled out, or on an independent record of sediment volume
  (Macrostrat). That is a natural next pre-registration.
- Range-through richness is not sample-standardised (registered limit). The AR(1) null treats
  unevenly spaced stages as evenly spaced. The 5-520 Ma window holds only about eight cycles.
- Hand checks: the window holds 91 stages and H6 n = 91; `declining_62` in `series62.csv` was
  re-derived for Stage 3 (mid 517.75 Ma, 6.05 Myr after the 523.8 Ma maximum: declining) and
  Jiangshanian (492.6 Ma, 31.2 Myr after it: rising), both as written.

## Files

`results/prereg2/`: `results62.json` (every number above), `series62.csv` (per-stage
richness, collections, extinction rate, phase), `periodogram_div.csv`, `provenance62.json`
(PBDB query URLs, md5 sums, times).

Code Apache-2.0; text CC BY 4.0. MV2, Munim, Inc., mv2@munimtech.com
