# Pre-result note 2 for PREREG-2 (2026-11-17, about 19:40 UTC)

Written by MV2 of Munim, Inc. (made with AI agents working in his name) **before the job
`extinction-cycle-62-v1` had produced any output**. Nothing in `PREREG-2.md` is changed.

## 1. The reading is now code

`grade_cycle2.py` (tests: `test_grade_cycle2.py`, standard library only) turns the job's
`results62.json` into the Q1-Q9 table, the H6 decision and the registered headline:

| Row | Key in results62.json | TRUE when |
|---|---|---|
| Q1 | `sanity.passed` | true |
| Q2 | `H6_range_through.r` | > 0 |
| Q3 | `H6_range_through.r`, `.p_ar1` | r > 0 and p_ar1 < 0.05 |
| Q4 | `H6_range_through.phase_offset_myr` | <= 8.0 |
| Q5 | `H6_periodogram.global_peak_in_band` | true (58-66 Myr) |
| Q6 | `H6_periodogram.p_band_ar1` | < 0.05 |
| Q7 | `short_lived.r`, `.p_ar1` | r > 0 and p_ar1 < 0.05 |
| Q8 | `log_collections.r`, `.p_ar1` | r > 0 and p_ar1 < 0.05 |
| Q9 | `H7_mu3t_declining.p_random_phase` | < 0.05 |

Every threshold needs a finite number. Headline: `PREDICTS_DIVERSITY` if Q3 holds and Q8
fails; `PREDICTS_DIVERSITY_AND_SAMPLING` if Q3 and Q8 both hold (PREREG-2: "if Q8 also
holds, we will say plainly that the same rhythm is in the sampling record");
`NOT_SUPPORTED` if Q3 fails; `SANITY_FAILED` if Q1 fails. For `NOT_SUPPORTED` the grader
also says what failed (PREREG-2: "whether the period or only the phase failed"):
`PERIOD_FAILED` (neither Q5 nor Q6), `PHASE_FAILED` (Q5 or Q6, but not Q4) or `WEAK`
(Q5 or Q6, and Q4, but H6 not significant).

## 2. What we have seen: the PREREG-2 data are already public

PREREG-1's job finished first and its per-stage table is now published
(`results/prereg1/series.csv`, see RESULTS-1.md). We read that table in full before
writing this note. It reads the same PBDB download (same queries, same cached files) and
the same stage assignment as PREREG-2's job. It holds Foote's four boundary counts per
stage (Nbt, NbL, NFt, NFL), and **their sum is exactly PREREG-2's range-through richness**
(every genus with first stage <= i <= last stage). So the main PREREG-2 series is in
effect already public. The H7 series (three-timer extinction rate, PREREG-1's window) is
published too.

We did **not** compute any PREREG-2 statistic: no detrended richness, no correlation with
the 62-Myr template, no periodogram, no H7 contrast. The decision rule, the thresholds and
our priors were all fixed in PREREG-2.md long before, and this note changes none of them.
But a reader should know that the grader was written after the data could be seen.

One row is settled by those counts, and we report it now:

| Check (S2) | Range-through richness (Nbt+NbL+NFt+NFL) | Holds? |
|---|---|---|
| Induan < Changhsingian | 438 < 861 | yes |
| Danian < Maastrichtian | 1406 < 1601 | yes |
| Rhuddanian < Katian | 720 < 1559 | yes |

So Q1 (sanity gate S2) will pass unless the job reads different data from PREREG-1's
(for example, if the cached files were replaced). We will compare its series62.csv with
these counts when it finishes.

MV2, Munim, Inc., mv2@munimtech.com
