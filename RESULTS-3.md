# RESULTS-3: the 27.5-Myr land-tetrapod extinction cycle, out of sample

Graded 2026-11-17 from the output of job `extinction-cycle-tetra-v2` (ran 20:28-20:32 UTC).
Plan: [PREREG-3.md](PREREG-3.md) (committed 09:20 UTC) and
[PREREG-3-AMENDMENT-1.md](PREREG-3-AMENDMENT-1.md) (12:17 UTC, before any tetrapod data).
Reading: [grade_tetra.py](grade_tetra.py) and [PRE-RESULT-NOTE-3.md](PRE-RESULT-NOTE-3.md)
(committed 20:11 UTC, before the job started). The grader's rules were applied to
`results_tetra.json` by hand, exactly as coded; `test_graded_results.py` runs the grader
itself on the published output and checks that it gives the same labels.
Made with AI agents working in MV2's name.

## Verdict

**DATA_LIMITED.** The registered sanity gate T1 needs at least 40 of the 57 stages from the
Gzhelian to the Piacenzian (stage tops 2.5-300 Ma) to have a defined three-timer extinction
rate. Only **38** do. As registered, nothing is claimed from the primary test.

Why: the three-timer rate needs genera sampled in three consecutive stages, and the
sampling-completeness term needs part-timers in the next stage. The land record is too
patchy for that in 19 stages, most of the Permian-Triassic transition (Artinskian, Kungurian,
Roadian, Wordian, Wuchiapingian, Changhsingian, Induan to Anisian), the Early-Middle Jurassic
(Sinemurian to Bathonian), the early Late Cretaceous (Cenomanian to Coniacian) and the
Piacenzian, the last stage in the window. Of 127,583
PBDB land-tetrapod genus occurrences, 89,324 (70%) fit inside one stage; 38,259 span a stage
boundary (mostly regional land-mammal ages) and were dropped, as registered.

**For the record only (not a registered claim):** no version of the series shows the cycle.

- Primary series (38 stages): r = 0.08 with the 7.0 + 27.5 n Ma ephemeris, p_AR1 = 0.31.
- With the nine stages holding Rampino et al.'s events left out (30 stages): r = -0.18, p = 0.84.
- Strongest period between 10 and 100 Myr: 15.8 Myr; the 26-29 Myr band p_AR1 = 0.73.
- **Foote's boundary-crosser rate q**, a registered secondary series that is defined in all 57
  stages and passes its own sanity check: **r = 0.007, p_AR1 = 0.48**; events left out
  (48 stages): r = -0.15, p = 0.85. Its own registered label is NOT_SUPPORTED.

So the only place a 27.5-Myr land-tetrapod cycle appears is the published list of event ages
it was fitted to. We will not call that a refutation, because the registered primary test
could not be run at adequate size; but neither the primary series nor the complete secondary
series gives it any support.

## Graded predictions

T1 failed, so as registered the Brier score covers T1 only: **0.490** over 1 row (prior 70%,
outcome no). Rows T2-T11 are shown for the record.

| # | Prediction | Prior | Holds? |
|---|---|---|---|
| T1 | Sanity gate passes (>= 40 window stages with a rate; Maastrichtian in the top quartile) | 70% | no |
| T2 | H1 r > 0 (detrended mu vs cos(2 pi (t - 7.0)/27.5)) | 65% | YES |
| T3 | H1 r > 0 and p_AR1 < 0.05 | 25% | no |
| T4 | H2 (leave the nine event stages out) p_AR1 >= 0.05 | 85% | YES |
| T5 | H1 p_phase < 0.05 (the published phase beats random phases) | 20% | no |
| T6 | The global periodogram peak (10-100 Myr) is outside 26-29 Myr | 80% | YES |
| T7 | H4 26-29 Myr band p_AR1 >= 0.05 | 75% | YES |
| T8 | H6 marine ephemeris (9.5 + 27.16 n Ma) r > 0 | 60% | YES |
| T9 | H7 stage-duration artifact p_AR1 >= 0.05 | 80% | YES |
| T10 | Non-mammal series has the same sign of r as H1 | 65% | YES |
| T11 | H5 event test p_binom >= 0.05 | 70% | YES |

Status checks: the job's decision (DATA_LIMITED) agrees; H1 n = 38 = window stages; H2 n = 30
= 38 minus the 8 event stages in the window (Capitanian, Carnian, Norian, Rhaetian,
Tithonian, Maastrichtian, Priabonian, Tortonian), all dropped. The ninth event stage, the
Changhsingian, has no defined rate. Status: GRADED_DATA_LIMITED.

## Key numbers

| Series | n | r | p_AR1 | p_phase |
|---|---|---|---|---|
| H1 all tetrapods, 7.0 + 27.5 n | 38 | 0.083 | 0.306 | 0.344 |
| H2 event stages left out | 30 | -0.178 | 0.837 | 0.695 |
| H3 normal scores | 38 | 0.099 | 0.264 | 0.267 |
| H6 marine ephemeris 9.5 + 27.16 n | 38 | 0.056 | 0.367 | 0.380 |
| H7 stage durations (artifact check) | 38 | -0.151 | 0.813 | 0.678 |
| Non-mammal tetrapods | 38 | 0.091 | 0.283 | 0.231 |
| Window 2.5-260 Ma | 35 | 0.032 | 0.433 | 0.429 |
| Foote q, all 57 stages | 57 | 0.007 | 0.479 | 0.469 |
| Foote q, event stages left out | 48 | -0.152 | 0.848 | 0.717 |

- H5 event test: 2 of 7 detrended-rate peaks within 3 Myr of a predicted maximum (chance
  0.22 each), p_binom = 0.47.
- Timescale jitter: r 0.044 / 0.081 / 0.116 (5-50-95%).
- Part A (the published events themselves; not evidence, since they were used to fit the
  cycle): 6 of 9 lie within 3 Myr of 7.0 + 27.5 n Ma.

## Hand checks

Three-timer rates re-derived from the per-stage counts in `series_tetra.csv` with
mu_i = ln(2t_i / 3t_i) + ln(Ps_{i+1}), Ps = 3t / (3t + pt):
Maastrichtian ln(137/27) + ln(2/2) = 1.624 (file 1.624); Rhaetian ln(26/3) + ln(1/1) = 2.159
(file 2.159); Capitanian ln(4/1) + ln(3/4) = 1.099 (file 1.099); Artinskian: Ps of the
Kungurian is 0/0, so the rate is undefined (file nan).

## What this means, and its limits

The land-tetrapod claim (Rampino, Caldeira & Zhu 2021) rests on ten event ages. This is, as
far as we can find, the first attempt to test it on an occurrence-based extinction-rate
series. At stage resolution the PBDB land record is too sparse for the registered estimator in
a third of the window. A finer test needs a rate that tolerates gaps (e.g. gap-filler or
capture-recapture rates) or bins built from land-mammal ages; either would be a new,
separately registered test, not a rescue of this one.

## Files

`results/prereg3/`: `results_tetra.json`, `series_tetra.csv` (per-stage counts and rates),
`periodogram_tetra.csv`, `provenance_tetra.json` (PBDB query URLs, md5 sums, parse counts).

Code Apache-2.0; text CC BY 4.0. MV2, Munim, Inc., mv2@munimtech.com
