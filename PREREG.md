# PREREG-1: an out-of-sample test of the 27-Myr extinction cycle

Registered 2026-11-17, before any of the data below was downloaded or analysed.
Made with AI agents working in MV2's name (MV2, Munim, Inc., mv2@munimtech.com).
Licence: Apache-2.0 (code), CC BY 4.0 (text).

## The claim being tested

Raup & Sepkoski (1984) reported that marine extinctions over the last 250 Myr
peak every ~26 Myr. Melott & Bambach (2010, 2013, 2014) extended the claim to
the whole Phanerozoic at ~27 Myr, from the Sepkoski compendium and from John
Alroy's 2008 Paleobiology Database (PBDB) rates on the 2012 geologic timescale
(GTS2012). Their Paleobiology 2014 paper (40:177-196, arXiv:1310.4712) fixes an
ephemeris: **maxima every 27.16 Myr, the most recent at 9.5 Ma**, and reports
that 10 of 19 marked extinction intervals lie within +/-3 Myr of those maxima
(p ~ 0.004 by their binomial count).

Two things are new since that ephemeris was published, and neither was available
when it was fitted:

1. **The data.** The 2026 PBDB holds far more marine occurrences than Alroy's
   2008 download, and the 2014 PBDB series used 48 bins averaging ~11 Myr, so a
   27-Myr period sat close to its Nyquist limit (~21 Myr). Single-stage bins
   (median ~4-5 Myr) resolve it properly.
2. **The timescale.** The Geologic Time Scale 2020 (as served by PBDB) moved
   many boundaries, e.g. Jurassic/Cretaceous 145.0 -> 143.1 Ma, Permian/Triassic
   252.2 -> 251.902 Ma, Toarcian base 182.7 -> 184.2 Ma.

A fixed, published ephemeris tested on new data is a genuine out-of-sample
prediction. If the cycle is real, today's record should show extinction
maxima at 9.5 + 27.16 n Ma. If the cycle is an artefact of coarse bins, of
the old timescale or of a few large events, it should not.

**What would be remarkable:** a confirmed, phase-stable 27-Myr pulse in Earth's
extinction record would point to an astronomical clock (the Sun's vertical
oscillation through the Galactic disc has been proposed) and would be one of
the most important regularities in the history of life. A clean failure would
retire a 40-year-old claim that is still cited.

## Prior art (searched 2026-11-17)

- Raup & Sepkoski 1984 (PNAS 81:801); Stigler & Wagner 1987/88 (timescale
  artefact critique).
- Melott & Bambach 2010 (MNRAS Letters, "Nemesis reconsidered"), 2013 (ApJ
  773:6), 2014 (Paleobiology 40:177): periodicity survives the GTS2012 revision.
- Erlykin et al. 2017/2018 (Palaeontology): events not independent; evidence
  not compelling. Melott & Bambach reply (Palaeontology, doi:10.1111/pala.12322).
- Rampino & Caldeira 2015 (MNRAS 454:3480): 27.0 +/- 0.7 Myr in extinctions
  and craters, most recent maximum ~11.8 Ma.
- Lieberman & Melott 2007 (PLoS ONE): 62 Myr robust, 27 Myr ambiguous.
- We found **no published test on GTS2020, none at single-stage resolution on
  the current PBDB, and none on Fan et al. (2020)'s independent ~26-kyr
  resolution Paleozoic record** (Science 367:272; data CC0 at
  zenodo.org/records/4977790).

## Data (part B, primary)

PBDB `occs/list.csv`, one request per group, `base_name` in {Porifera,
Cnidaria, Bryozoa, Brachiopoda, Mollusca, Echinodermata, Trilobita,
Graptolithina}, with `envtype=marine&idreso=genus&idqual=genus_certain&
pres=regular&show=class`. Brachiopoda, Mollusca, Echinodermata and Trilobita
are required; if one fails, the run stops with status `data_incomplete`.

Taxon key = phylum + accepted genus. An occurrence goes to the GTS2020 stage that
contains the midpoint of its age range, and only if the whole range lies inside
that stage (0.1 Myr tolerance). Stages are PBDB's international ages (GTS2020),
plus Pridoli as one bin and the Quaternary as one bin (96 bins; `cycle.py`).

## Metrics

- **Primary:** Alroy (2008) three-timer extinction rate at the top of each stage,
  mu_i = ln(2t_i/3t_i) + ln(Ps_{i+1}), dated to the stage's top (GTS2020).
- Secondary: Foote (2000) boundary-crosser rate q; the proportion of genera
  last seen in the stage (Melott & Bambach's "magnitude").

Analysis window: stage tops from 2.5 to 475 Ma (Melott & Bambach truncate at
475 Ma). Each series is detrended with a least-squares straight line.

## Sanity gate (S1)

At least 4 of the Big Five (end-Ordovician = Katian or Hirnantian; Frasnian;
Changhsingian; Rhaetian; Maastrichtian) must be among the 10 highest mu values
in the window. If not, the data build is judged broken and no hypothesis is graded.

## Hypotheses and tests

**H1 (primary, decides the headline).** With mu, the Pearson correlation r
between the detrended series and cos(2 pi (t - 9.5)/27.16) is positive, and
its one-sided p-value against 20,000 AR(1) surrogates (lag-1 coefficient taken
from the detrended series, surrogates detrended the same way) is below 0.05.
Result: SUPPORTED / NOT_SUPPORTED / SANITY_FAILED.

Reported alongside (not deciding):
- p_phase: the same data against 20,000 templates with uniformly random phase.
- The least-squares phase at 27.16 Myr and its offset from 9.5 Ma.
- The same test for the MB2017 (27.19 Myr, 11.3 Ma) and RC2015 (27.0 Myr,
  11.8 Ma) ephemerides.
- **H2:** Lomb-Scargle periodogram (periods 10-100 Myr): the period of the
  global peak, and the maximum power in the 26-28.5 Myr band against AR(1)
  surrogates (band-maximum statistic).
- **H3:** event test: local maxima of detrended mu at or above its 75th
  percentile; the count within +/-3 Myr of predicted maxima vs binomial
  p0 = 6/27.16.
- **H4 (robustness):** H1 on normal scores of mu (ranks), and H1 with the six
  Big Five bins removed.
- **H5:** H1 r under timescale jitter (each top age +N(0, sigma), sigma =
  0.25 / 0.75 / 1.0 Myr for Cenozoic / Mesozoic / Paleozoic; 1,000 draws).
- Secondary metrics q and proportion: H1-H3 repeated.

## Predictions (our honest priors, written before the data)

| # | Prediction | Our probability |
|---|---|---|
| P1 | Sanity gate S1 passes | 85% |
| P2 | H1 r > 0 (four of the Big Five lie near the MB2014 maxima) | 75% |
| P3 | H1 SUPPORTED (r > 0 and p_AR1 < 0.05) | 40% |
| P4 | With the six Big Five bins removed, p_AR1 >= 0.05 | 75% |
| P5 | On normal scores (ranks), p_AR1 >= 0.05 | 65% |
| P6 | The global periodogram peak (10-100 Myr) is NOT in 26-28.5 Myr | 75% |
| P7 | The band-max p_AR1 (H2) >= 0.05 | 70% |
| P8 | The fitted phase at 27.16 Myr is within +/-3 Myr of 9.5 Ma | 45% |
| P9 | H3 binomial p >= 0.05 | 60% |

The headline reads "the published 27-Myr ephemeris predicts today's stage-level
extinction record" only if P3 holds AND P4 fails. If P3 holds only because of the
Big Five (P4 holds), we will say so: the cycle then rests on where a few mass
extinctions happen to fall, not on a rhythm in background extinction.

## Part A (already computed, not a prediction)

Re-dating Melott & Bambach's 19 Table 1 events from GTS2012 to the GTS2020 ages
of the same boundaries (`cycle.MB2014_EVENTS`; the Toarcian event and the
mid-Caradoc event are approximate, and neither is near a maximum either way)
leaves the count at 10 of 19 within +/-3 Myr. The re-dating alone does not break
the window result. This is not an independent test, because the ephemeris was fitted
to these same events.

## Part C (registered now, run in a later job): Fan et al. 2020

Fan et al. (2020) date 11,268 Cambrian-Early Triassic species from 3,112
sections (mostly South China) by constrained optimisation (CONOP.SAGA) at about
26 kyr resolution. It shares no data with the Sepkoski compendium or the PBDB
rates the ephemeris was fitted to. Plan, fixed now: from the composite species
ranges (first and last appearance ages in Ma), build 1-Myr bins over 250-470 Ma.
Per bin, compute the proportion of species last appearing in the bin among all
species ranging through it, and Foote's q. Then run H1 (MB2014 ephemeris; AR(1)
null from the 1-Myr series) and H2 on the proportion, and on q as secondary.
Predicted maxima in this window: 253.94, 281.10, 308.26, 335.42, 362.58,
389.74, 416.90 and 444.06 Ma. Our priors: r > 0 (55%); p_AR1 < 0.05 (25%).
Job v1 only lists the files in the archive and the first lines of each
(`fan_inventory.txt`) and computes no statistics, so the parser can be written
before any Fan number is seen.

## What would change our mind

If H1 holds, P4 fails (the cycle survives without the Big Five) and Part C
also shows r > 0 with p_AR1 < 0.05, that is strong out-of-sample support for a
real 27-Myr rhythm, and we will say so plainly. If H1 fails at stage resolution,
the 2014 result was probably shaped by the ~11-Myr bins and the old timescale.

MV2, Munim, Inc., mv2@munimtech.com
