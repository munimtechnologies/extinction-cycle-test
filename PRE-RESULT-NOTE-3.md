# PRE-RESULT-NOTE-3: how PREREG-3 will be read, fixed before its result

Dated 2026-11-17 (about 20:20 UTC), by MV2 of Munim, Inc. Made with AI agents
working in his name.

PREREG-3 tests whether land-tetrapod extinction follows the 27.5-Myr clock of
Rampino, Caldeira & Zhu (2021), maxima at 7.0 + 27.5 n Ma, **outside** the
nine events the clock was fitted to. Its job (`job_tetra.py`, compute job
`extinction-cycle-tetra-v2`) is queued and has produced no output. This note and
`grade_tetra.py` fix, in code, how its `results_tetra.json` will be turned into
a verdict, so that the reading cannot be bent after the numbers are known. The
same was done for tests 1 and 2 (`grade_cycle1.py`, `grade_cycle2.py`).

## What had been seen when this was written

- **No tetrapod occurrence data had been analysed.** The only tetrapod rows
  ever fetched were the format checks of PREREG-3-AMENDMENT-1 (`limit=3` and
  `limit=5`, a handful of genus names); nothing was computed from them.
- Our own published outputs reveal nothing about this series: test 1
  (`results/prereg1/`) and test 2 use **marine** invertebrates only
  (`envtype=marine`); the tetrapod queries use `envtype=terrestrial` and
  different taxa, so no occurrence is shared.
- Test 1's verdict is public (`RESULTS-1.md`, RESTS_ON_MASS_EXTINCTIONS: the
  marine 27-Myr signal disappears without the Big Five). That is relevant
  background for this test but not data from it. **The T1-T11 priors in
  PREREG-3 are not changed.** They were registered before test 1 was graded.
- Part A of PREREG-3 (where the nine published event ages fall on the
  ephemeris: 6 of 9 within 3 Myr) uses only the published ages and was
  already stated in PREREG-3. It is not evidence, because the clock was
  fitted to those events.

## Reading choices fixed now

| Row | Read as (keys written by `tetra.analyze`) |
|---|---|
| T1 | `primary.sanity.passed` is true |
| T2 | `primary.H_rcz.r` > 0 |
| T3 | `primary.H_rcz`: r > 0 and p_ar1 < 0.05 |
| T4 | `primary.leave_events_out.p_ar1` >= 0.05, read literally as registered |
| T5 | `primary.H_rcz.p_phase` < 0.05 |
| T6 | `primary.periodogram.global_peak_in_band` is false (band 26-29 Myr) |
| T7 | `primary.periodogram.p_band_ar1` >= 0.05 |
| T8 | `primary.H_mb2014.r` > 0 |
| T9 | `primary.duration_artifact.p_ar1` >= 0.05 |
| T10 | non-mammal r and all-tetrapod r both finite, non-zero, same sign; FALSE if the non-mammal series was too short to test |
| T11 | `primary.event_test.p_binom` >= 0.05 |

- A NaN never satisfies an inequality (neither "< 0.05" nor ">= 0.05").
- A row whose block is missing is n/a and left out of the Brier score.
- If the sanity gate T1 fails, the result is DATA_LIMITED: rows are shown for
  the record, nothing is claimed, and the Brier score covers T1 only.

## Decision and headline

The decision is the registered rule: DATA_LIMITED / SUPPORTED_BEYOND_EVENTS /
SUPPORTED_ONLY_BY_EVENTS / NOT_SUPPORTED. The headline adds PREREG-3's
condition for saying the land record "independently supports" the clock:

| Decision | Headline |
|---|---|
| SUPPORTED_BEYOND_EVENTS, no duration artifact (T9 holds) | **INDEPENDENT_SUPPORT** |
| SUPPORTED_BEYOND_EVENTS, duration artifact (T9 fails) | SUPPORT_BUT_DURATION_ARTIFACT |
| SUPPORTED_ONLY_BY_EVENTS | RESTS_ON_FITTED_EVENTS |
| NOT_SUPPORTED | NOT_SUPPORTED |
| DATA_LIMITED | DATA_LIMITED |

Only INDEPENDENT_SUPPORT would be reported as new evidence for a pacemaker of
extinction on land. Our registered priors put it near 4% (T3 25% times the
complement of T4, 15%, if the two were independent; they are not quite, since
H2 is a subset of H1's data).

## Consistency checks (they change the status, never the label)

- CHECK_CODE: the grader's decision differs from the one written by the job,
  or a registered event stage (Tortonian, Priabonian, Maastrichtian,
  Tithonian, Rhaetian, Norian, Carnian, Changhsingian, Capitanian) has a rate
  in the window but was not dropped by the leave-events-out test.
- CHECK_DATA: the H1 sample size differs from the number of window stages,
  or the H2 sample size differs from that number minus the dropped stages.

Any CHECK_* status is published as it is, with the reason found.

## Files

- `grade_tetra.py`: the reading above, standard library only.
  `python grade_tetra.py --results results_tetra.json --out-dir .` writes
  `graded3.json` and `GRADED-3.md`.
- `test_grade_tetra.py`: 18 offline tests (`python -m unittest -v test_grade_tetra`).

Code Apache-2.0; text CC BY 4.0. MV2, Munim, Inc., mv2@munimtech.com
