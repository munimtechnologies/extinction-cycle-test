# Pre-result note 1 for PREREG-1 (2026-11-17, about 18:45 UTC)

Written by MV2 of Munim, Inc. (made with AI agents working in his name) **before the job
`extinction-cycle-v1` had produced any output**. No extinction rate from this project has
been computed or seen. Nothing in `PREREG.md` is changed.

## 1. The reading is now code

`grade_cycle1.py` (tests: `test_grade_cycle1.py`, standard library only) turns the job's
`results.json` into the P1-P9 table, the H1 decision and the registered headline:

- `PREDICTS_BEYOND_MASS_EXTINCTIONS`: P3 holds (r > 0, p_AR1 < 0.05) **and** P4 fails (the
  cycle survives with the six Big Five bins removed). Only this outcome supports the
  headline "the published 27-Myr ephemeris predicts today's stage-level extinction record".
- `RESTS_ON_MASS_EXTINCTIONS`: P3 holds, but without the Big Five p_AR1 >= 0.05.
- `NOT_SUPPORTED` if P3 fails; `SANITY_FAILED` if the sanity gate fails. In that case, as
  PREREG.md says, no hypothesis is graded and the Brier score covers P1 only.

A "p >= 0.05" row is true only if the p value is a finite number. The grader compares its
decision with the job's `decision_H1` (`CHECK_CODE` if they differ). It also checks the
GTS2020 top ages of the six Big Five bins against the registered stage table
(`CHECK_DATA`).

## 2. What public inputs already say: where the Big Five sit on the ephemeris

The MB2014 ephemeris (maxima at 9.5 + 27.16 n Ma) and the GTS2020 stage tops used by the
job are both public. So the position of each Big Five bin on the template can be worked out
now. It does not depend on any extinction rate:

| Big Five bin | GTS2020 top (Ma) | Nearest predicted maximum (Ma) | Offset (Myr) | Template cos |
|---|---|---|---|---|
| Katian | 445.2 | 444.06 | 1.14 | +0.97 |
| Hirnantian | 443.1 | 444.06 | 0.96 | +0.98 |
| Frasnian | 372.15 | 362.58 | 9.57 | -0.60 |
| Changhsingian | 251.902 | 253.94 | 2.04 | +0.89 |
| Rhaetian | 201.4 | 199.62 | 1.78 | +0.92 |
| Maastrichtian | 66.0 | 63.82 | 2.18 | +0.88 |

Five of the six Big Five bins lie within 2.2 Myr of a predicted maximum. A random age has
about a 16% chance of falling that close. This is not evidence for the cycle: the ephemeris
was fitted by Melott & Bambach to extinction events that include these, so it is in-sample.
What it means for reading the result:

- These bins will almost certainly have some of the highest three-timer rates (the sanity
  gate requires at least 4 of the 5 events in the top 10). They all sit at template values
  of +0.88 or more, so they alone push H1's r upwards. A positive r (P2), and even a
  significant one (P3), is largely built in.
- The informative test is therefore **P4**, i.e. H1 with the six Big Five bins removed, and
  the normal-score version (P5), which limits how much any single bin can weigh. That is why
  the registered headline needs P3 to hold **and** P4 to fail.
- The Frasnian (Late Devonian) bin sits close to an anti-phase point. If the Late Devonian
  event is strong in the data, it pulls r down.

## What was seen before writing this note

The project's code and PREREG, the GTS2020 stage table in `cycle.py`, and the published
ephemeris. No PBDB-derived number.
