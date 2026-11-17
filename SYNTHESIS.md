# What drives the history of life? Five pre-registered tests, one pattern

2026-11-17. MV2, Munim, Inc. Made with AI agents working in his name.

Over one day we pre-registered, ran and graded five tests of famous claims about the
500-million-year marine and land fossil record. Each plan, decision rule and set of honest
probabilities was committed to GitHub before any data were analysed, and each grader was
committed before its job produced output. All five results are published whichever way they
came out. The data are the 2026 Paleobiology Database (PBDB), the 2020 Geologic Time Scale at
single-stage resolution, and PhanDA, the 485-Myr global temperature history of Judd et al.
(2024, Science).

| # | Claim | Verdict | Key number |
|---|---|---|---|
| 1 | Marine extinctions recur every ~27 Myr (Raup & Sepkoski 1984; Melott & Bambach 2014 ephemeris 9.5 + 27.16 n Ma) | **RESTS_ON_MASS_EXTINCTIONS** ([RESULTS-1](RESULTS-1.md)) | r = 0.24 (p = 0.021) with the Big Five; r = 0.12 (p = 0.17) without |
| 2 | Marine genus diversity rises and falls every ~62 Myr (Rohde & Muller 2005; phase 27.8 + 62 n Ma) | **SUPPORTED, but also in sampling** ([RESULTS-2](RESULTS-2.md)) | r = 0.58 (p < 0.0001), peak period 60.6 Myr, phase 2.6 Myr off; collections r = 0.31 (p = 0.015) |
| 3 | Land-tetrapod extinctions recur every 27.5 Myr (Rampino, Caldeira & Zhu 2021) | **DATA_LIMITED** ([RESULTS-3](RESULTS-3.md)) | 38 of 57 stages have a rate (gate: 40); the complete Foote-rate series gives r = 0.007 (p = 0.48) |
| 4 | Temperature change predicts extinction; 5.2 C is a mass-extinction threshold (Song et al. 2021) | **REPLICATED_ONLY_WITH_MASS_EXTINCTIONS** ([extinction-temperature-test](https://github.com/munimtechnologies/extinction-temperature-test/blob/main/RESULTS.md)) | rho = 0.21 (p = 0.033) with the Big Five; 0.14 (p = 0.11) without; 9 boundaries reach 5.2 C, 2 are Big Five |
| 5 | A warmer Earth holds more marine genera (Mayhew et al. 2012) | **NOT_SUPPORTED** ([diversity-temperature-test](https://github.com/munimtechnologies/diversity-temperature-test/blob/main/RESULTS.md)) | levels rho = -0.06; changes rho = +0.05; origination rises with temperature (rho = 0.28, p = 0.027) |

## The pattern

**Take away the Big Five mass extinctions and the famous signals in extinction go quiet.**
Neither the 27-Myr extinction clock (test 1) nor the temperature-extinction link (test 4)
reaches significance in background extinction on today's data at stage resolution, and the
land-tetrapod clock (test 3) is visible only in the event list it was fitted to. The mass
extinctions are real and are tied to large temperature changes, but they do not make a
rhythm, and the 5.2 C rule does not pick them out on PhanDA.

**The one signal that survives everything we threw at it is the 62-Myr rhythm in
diversity**, and it survives very strongly: a phase published in 2014 predicts a 2026 database
on a timescale it never saw. But the same rhythm is in how many fossil collections each stage
has. Whatever beats every 62 Myr (sea level, sedimentation, or something that drives both
rocks and life), it is written into the rock record as well as into the count of genera.

**Global temperature does not set how many genera the ocean holds** (test 5), though warm
intervals may produce more new genera (a registered side-result that needs its own test).

## How well we forecast

Brier scores of our registered priors (0.25 = coin flip): test 1 0.194; test 2 0.329
(we were too sceptical of the 62-Myr cycle); test 3 0.490 (one row: we expected the land record
to be dense enough); test 4 0.183; test 5 0.181.

## What would come next

1. The 62-Myr cycle with sampling partialled out, and against an independent sediment record
   (Macrostrat), registered before either is computed. This is the question that decides
   whether the cycle belongs to life or to rocks. Boulila et al. (2023, PNAS; 2025,
   Communications Earth & Environment) argue for a cycle shared by rocks and life, driven by
   Earth's interior; a fixed-phase test against their sediment series would put that to a
   real prediction. Their robust ~36-Myr cycle also matches the ~38-Myr peak in our test 1.
2. A gap-tolerant land-tetrapod rate (gap-filler or capture-recapture) to give test 3 a
   series that is complete enough to grade.
3. Origination vs temperature with mass-extinction aftermaths removed (test 5's side-result).

## Limits that apply to all five

PBDB sampling is uneven in space and time; stage-level binning drops occurrences whose ages
cross a boundary; AR(1) and circular-shift nulls are approximations for unevenly spaced
series; 70-90 stages give modest power for weak effects. Every number comes from code and
public data in these repositories; corrections are welcome as GitHub issues.

Code Apache-2.0; text CC BY 4.0. MV2, Munim, Inc., mv2@munimtech.com
