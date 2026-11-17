# Results 1: the 27-Myr extinction cycle on today's fossil record (PREREG-1)

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-17. Made with AI agents working in MV2's name.
Code: Apache-2.0. Text and derived tables: CC BY 4.0 (the underlying occurrences are from
the Paleobiology Database, CC BY 4.0; please cite the PBDB and its contributors).

## The answer in one paragraph

Melott & Bambach's published ephemeris (extinction maxima every 27.16 Myr, the latest at
9.5 Ma) does line up with the 2026 Paleobiology Database marine record at single-stage
resolution on the 2020 timescale: r = 0.24, p = 0.021 against red-noise surrogates, and the
best-fitting phase lands at 9.51 Ma. **But the whole signal comes from the mass extinctions
the ephemeris was fitted to.** With the six Big Five stages removed, the correlation halves
(r = 0.12) and is not significant (p = 0.17), the fitted amplitude falls to about a third,
and the strongest periodicity in the full record is ~38 Myr, not 27 (the 26-28.5 Myr band
is not significant, p = 0.10). Pre-registered verdict: **RESTS_ON_MASS_EXTINCTIONS**. On
today's data there is no detectable 27-Myr rhythm in background extinction.

## How this was graded

The plan, the decision rule and our probability for each outcome were published in
[PREREG.md](PREREG.md) before any data were downloaded. The mechanical reading was fixed in
[grade_cycle1.py](grade_cycle1.py) and [PRE-RESULT-NOTE-1.md](PRE-RESULT-NOTE-1.md) before
the job produced any output. The job (`extinction-cycle-v1`, 2026-11-17 18:07-18:12 UTC,
the code in this repository unchanged) wrote [results/prereg1/results.json](results/prereg1/results.json).

The table below applies `grade_cycle1.py`'s rules to that file **by hand, exactly as
coded** (no code was run for this step). Anyone can reproduce it with
`python grade_cycle1.py --results results/prereg1/results.json --out-dir .`.

- Status: **GRADED** (sanity gate passed; the decision agrees with the job's own
  `decision_H1`; the GTS2020 top ages of all six Big Five bins match the registered table)
- H1 decision: **SUPPORTED**
- Registered headline: **RESTS_ON_MASS_EXTINCTIONS** (P3 holds and P4 holds)

| # | Prediction | Our prior | Holds? | Number |
|---|---|---|---|---|
| P1 | Sanity gate S1 passes (>= 4 of the Big Five among the 10 highest mu3t bins) | 85% | YES | 5 of 5 |
| P2 | H1 r > 0 at the MB2014 ephemeris (27.16 Myr, 9.5 Ma) | 75% | YES | r = 0.239 |
| P3 | H1 SUPPORTED: r > 0 and p_AR1 < 0.05 | 40% | YES | p_AR1 = 0.0213 |
| P4 | With the six Big Five bins removed, p_AR1 >= 0.05 | 75% | YES | r = 0.116, p_AR1 = 0.168 |
| P5 | On normal scores (ranks), p_AR1 >= 0.05 | 65% | no | p_AR1 = 0.0471 |
| P6 | The global periodogram peak (10-100 Myr) is NOT in 26-28.5 Myr | 75% | YES | peak at 37.7 Myr |
| P7 | The band-max p_AR1 (H2) >= 0.05 | 70% | YES | p = 0.103 |
| P8 | The fitted phase at 27.16 Myr is within +/-3 Myr of 9.5 Ma | 45% | YES | 9.51 Ma (offset 0.01 Myr) |
| P9 | H3 binomial p >= 0.05 | 60% | no | 7 of 11 peaks hit, p = 0.0036 |

Brier score of our registered priors: 0.194 over 9 rows (a 50% forecaster scores 0.250).
We under-rated P3 (40%) and P8 (45%), and over-rated P5 and P9.

## Key numbers (all from results.json)

Series: Alroy three-timer extinction rate (mu3t) of marine invertebrate genera, one value per
GTS2020 stage, 83 stages with tops between 2.5 and 475 Ma (Floian to Zanclean). 641,280 PBDB
occurrences in eight phyla were assigned to a single stage (provenance.json).

| Test (MB2014 ephemeris unless stated) | n | r | p_AR1 | p_phase | fitted max (Ma) | events hit / peaks (p) |
|---|---|---|---|---|---|---|
| mu3t, all stages | 83 | 0.239 | 0.021 | 0.0008 | 9.51 | 7 / 11 (0.0036) |
| mu3t, Big Five removed | 77 | 0.116 | 0.168 | 0.213 | 6.69 | 4 / 10 (0.16) |
| mu3t, normal scores | 83 | 0.196 | 0.047 | 0.059 | 10.31 | 7 / 11 (0.0036) |
| mu3t, MB2017 (27.19 Myr, 11.3 Ma) | 83 | 0.206 | 0.041 | 0.161 | 9.29 | 4 / 11 (0.21) |
| mu3t, RC2015 (27.0 Myr, 11.8 Ma) | 83 | 0.252 | 0.014 | 0.093 | 10.63 | 6 / 11 (0.020) |
| Foote q | 84 | 0.200 | 0.043 | 0.076 | 10.61 | 6 / 11 (0.019) |
| Proportion last seen | 84 | 0.153 | 0.098 | 0.118 | 11.10 | 4 / 13 (0.32) |

- Periodogram of mu3t (10-100 Myr): strongest peaks at 37.7, 26.3, 12.2, 15.3 and 23.3 Myr.
  The 26-28.5 Myr band maximum (power 4.00) has p = 0.103 against AR(1) surrogates.
- Fitted amplitude at 27.16 Myr: 0.085 with all stages, 0.030 without the Big Five.
- Timescale jitter (H5): r stays between 0.206 and 0.265 (5-95%) and is positive in every draw,
  so the all-stage result is not fragile to dating error. Its weakness is elsewhere.
- Part A (re-dating Melott & Bambach's 19 events from GTS2012 to GTS2020): 10 of 19 within
  +/-3 Myr of a maximum on both timescales.

## What it means

1. **The agreement is real but in-sample.** The pre-result note showed, from public inputs
   alone, that five of the six Big Five stages sit within 2.2 Myr of a predicted maximum
   (the Frasnian is the exception, near an anti-phase point). These are the very events the
   ephemeris was fitted to. Their high extinction rates are enough to produce r = 0.24 and a
   fitted phase at 9.5 Ma.
2. **Between the mass extinctions, nothing ticks.** Without those six stages, the series
   correlates only weakly with the ephemeris (p = 0.17), the best phase drifts by 2.8 Myr,
   and only 4 of 10 local maxima fall near a predicted peak. One of those four is the
   Danian, the stage right after the end-Cretaceous extinction (and followed by the thinly
   sampled Selandian, 137 collections), so it is hardly an independent event.
   Three other moderate peaks do sit near predicted maxima (Cenomanian 93.9 Ma, Tithonian
   143.1 Ma, Toarcian 174.7 Ma), but not more often than chance at this sample size.
3. **The record's own strongest rhythm is not 27 Myr.** The periodogram peaks at 37.7 Myr;
   27 Myr is the second-highest peak and is not significant against red noise.
4. **Not a refutation.** We did not calibrate the test's power for a background cycle as
   weak as the residual amplitude, and stage-level bins (median ~4-5 Myr) still blur
   events. The honest statement is narrower: in the 2026 PBDB at stage resolution, the
   27-Myr ephemeris predicts where the big extinctions are, which it was built to do, and
   adds no detectable prediction for the rest of the record. This is in line with the
   doubts of Erlykin et al. (2017, Palaeontology), here with a pre-registered test on data
   and a timescale the ephemeris never saw.

## Part C (Fan et al. 2020): cannot be run as registered

PREREG.md planned to build 1-Myr bins from Fan et al.'s composite species ranges. The
job's inventory of the Zenodo archive (records/4977790, md5 verified; the paper's Dryad
record doi:10.5061/dryad.0gb5mkkx1 lists the same dataset zip name) shows it holds the CONOP inputs (sections, events,
`loadfile.dat`, `taxon.dic`, the run configuration), not the dated composite ranges. Those
came from a CONOP.SAGA optimisation run on a supercomputer, which we cannot reproduce.
Part C is therefore not run. No statistic was computed from the Fan files.

## Files

- [results/prereg1/results.json](results/prereg1/results.json): all test outputs
- [results/prereg1/series.csv](results/prereg1/series.csv): per-stage counts and rates (96 bins)
- [results/prereg1/periodogram_mu3t.csv](results/prereg1/periodogram_mu3t.csv): Lomb-Scargle power, 10-100 Myr
- [results/prereg1/provenance.json](results/prereg1/provenance.json): PBDB queries, bytes, md5, parse counts

## Hand checks

Two rows of series.csv were recomputed by hand from their counts with the registered
formula mu_i = ln(2t_i/3t_i) + ln(Ps_{i+1}): Changhsingian ln(464/66) + ln(0.7544) = 1.668,
and Rhaetian ln(338/99) + ln(0.7647) = 0.960, both matching the job.

MV2, Munim, Inc., mv2@munimtech.com
