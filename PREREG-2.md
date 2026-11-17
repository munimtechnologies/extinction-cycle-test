# PREREG-2: an out-of-sample test of the 62-Myr cycle in marine diversity

Registered 2026-11-17, before any PBDB data was seen by us. Job
`extinction-cycle-62-v1` reads the same PBDB download as PREREG-1's job
(same queries, same cache). Neither job had produced any output when this was
written.
Made with AI agents working in MV2's name (MV2, Munim, Inc., mv2@munimtech.com).
Licence: Apache-2.0 (code), CC BY 4.0 (text).

## The claim being tested

Rohde & Muller (2005, Nature 434:208) found a **62 +/- 3 Myr cycle** in
the number of marine animal genera in Sepkoski's compendium (cubic-detrended,
GTS2004). It explained 35% of the variance, and they wrote that the Big Five mass
extinctions "all occur during declining phases" of it. Melott & Bambach
(Paleobiology 2011a,b; Melott 2008 PLoS ONE) found the same period in the
Paleobiology Database (Alroy 2008, ~11-Myr bins) and Fossil Record 2, with phases
agreeing to within 1.6 Myr. Melott & Bambach 2014 (Paleobiology 40:177,
arXiv:1310.4712) fix the phase: **"The latest peak in the 62-Myr periodicity
occurred 27.8 Myr ago."** So diversity maxima fall at **27.8 + 62 n Ma**
(27.8, 89.8, 151.8, 213.8, 275.8, 337.8, 399.8, 461.8 Ma), and declining
half-cycles run from each maximum to 31 Myr later. They count 16 of 19 extinction
events in declining phases (p = 0.002).

Proposed causes include the Sun's vertical oscillation through the Galactic disc
(extragalactic cosmic rays; Medvedev & Melott 2007) and a ~60-Myr rhythm in
continental freeboard and sedimentation (Melott et al. 2012, J. Geol. 120:217;
Meyers & Peters 2011).

**What would be remarkable:** a 62-Myr rhythm that keeps its published phase on
data and a timescale it was never fitted to would be among the most striking
regularities in the history of life. It would also be strong evidence for an
external clock, astronomical or tectonic. A clean failure at stage resolution
would show the cycle depends on how the older compendia and timescales were put
together.

## Prior art (searched 2026-11-17)

The supporting papers are listed above. Lieberman & Melott 2007 (PLoS ONE) found
the 62-Myr signal robust but found no periodicity in origination and extinction
with stage boundaries. We found **no test of the 62-Myr ephemeris on GTS2020,
none at single-stage resolution on the current PBDB, and none fixing the
published phase in advance**.

## Data

The data are the same as PREREG-1: PBDB `occs/list.csv` marine invertebrate genus
occurrences (Porifera, Cnidaria, Bryozoa, Brachiopoda, Mollusca, Echinodermata,
Trilobita, Graptolithina). Each occurrence goes to a GTS2020 stage only if its
age range fits inside that stage (`cycle.py`). Diversity is **range-through genus
richness**: a genus counts in every stage from its first to its last stage. Each
value is dated to the stage midpoint. The window is midpoints from 5 to 520 Ma
(the window Melott & Bambach 2011 used for PBDB), and the series is
**cubic-detrended** (as Rohde & Muller and Melott & Bambach did).

## Sanity gate (S2)

Range-through richness must drop across the three best-known crises: Induan <
Changhsingian, Danian < Maastrichtian, Rhuddanian < Katian. If any check fails,
the build is judged broken and nothing is graded.

## Hypotheses

- **H6 (primary; decides the headline).** The Pearson r between cubic-detrended
  richness and cos(2 pi (t - 27.8)/62) is > 0, and its one-sided p against 20,000
  AR(1) surrogates (lag-1 taken from the detrended series, surrogates detrended
  the same way) is < 0.05. The result is SUPPORTED / NOT_SUPPORTED / SANITY_FAILED.
  Reported alongside: the random-phase p, the least-squares phase at 62 Myr,
  timescale jitter (as in PREREG-1), and r with a quadratic detrend.
- **Periodogram.** Lomb-Scargle on 20-200 Myr periods: the period of the global
  peak, and the 58-66 Myr band maximum against AR(1) surrogates.
- **Short-lived genera** (first-to-last stage midpoints 45 Myr apart or less, the
  group Rohde & Muller and Melott & Bambach found carries the signal): H6 and the
  periodogram repeated.
- **Sampling check.** H6 repeated on log(number of PBDB collections per stage)
  and on sampled-in-bin richness. If collections show the same cycle, the test
  cannot tell a biological rhythm from a rhythm in the rock record or in sampling.
- **H7 (declining-phase extinction).** The three-timer extinction rate (PREREG-1's
  series: stage tops, 2.5-475 Ma, linearly detrended, then normal scores) is
  higher in declining half-cycles than in increasing ones. The statistic is the
  difference in mean normal score. The p-value is one-sided, against 20,000
  ephemerides with a uniformly random phase.

## Predictions (our honest priors, written before the data)

| # | Prediction | Our probability |
|---|---|---|
| Q1 | Sanity gate S2 passes | 85% |
| Q2 | H6 r > 0 | 70% |
| Q3 | H6 SUPPORTED (r > 0 and p_AR1 < 0.05) | 40% |
| Q4 | The fitted 62-Myr phase is within +/-8 Myr of 27.8 Ma | 45% |
| Q5 | The global periodogram peak (20-200 Myr) lies in 58-66 Myr | 30% |
| Q6 | The 58-66 Myr band p_AR1 < 0.05 | 30% |
| Q7 | Short-lived genera: r > 0 and p_AR1 < 0.05 | 40% |
| Q8 | log(collections): r > 0 and p_AR1 < 0.05 (the cycle is also in sampling) | 30% |
| Q9 | H7 p_random_phase < 0.05 | 40% |

The headline reads "the published 62-Myr ephemeris predicts today's stage-level
diversity record" only if Q1 and Q3 hold. If Q8 also holds, we will say plainly
that the same rhythm is in the sampling record. If Q3 fails, the headline is that
the cycle does not survive the move to the 2026 PBDB at stage resolution on
GTS2020. We will then say whether the period or only the phase failed (Q4-Q6).

## Part A-ii (computed now, not a prediction)

Re-dating Melott & Bambach's 19 events to GTS2020 (`cycle.MB2014_EVENTS`) and
applying 27.8 + 62 n gives **16 of 19 in declining phases**, the same count they
report. On their GTS2012 ages, the same formula gives 15. The early Toarcian at
182.0 Ma falls 0.8 Myr past a half-cycle boundary, so their table must use a
slightly different phase or rounding. Two events sit within 0.2 Myr of a boundary
on GTS2020 (Toarcian, 183.0 Ma; late Moscovian, 307.0 Ma). This is not
independent evidence, because the phase was fitted to the same compendia.

## Limits fixed in advance

Range-through counts from raw occurrences are not sample-standardised. That is
why the sampling check is registered. The AR(1) null treats unevenly spaced
stages as evenly spaced (as in PREREG-1). The window holds only about 8 cycles.

MV2, Munim, Inc., mv2@munimtech.com
