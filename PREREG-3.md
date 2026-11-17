# PREREG-3: does land-animal extinction follow the claimed 27.5-Myr clock?

Registered 2026-11-17, before any PBDB tetrapod data were downloaded or looked
at. Code: `tetra.py`, `job_tetra.py`, tests `test_tetra.py` (statistics from
`cycle.py`, as registered in PREREG-1). Made with AI agents working in MV2's
name (MV2, Munim, Inc.).

## The claim

Rampino, Caldeira & Zhu (2021, *Historical Biology* 33:3084, "A 27.5-My
underlying periodicity detected in extinction episodes of non-marine tetrapods")
listed ten extinction episodes of land tetrapods (amphibians, reptiles, birds,
mammals) over the last 300 Myr. Circular-spectral and Fourier analysis of
their **ages** gave a period of ~27.5 Myr at 99% confidence. Their follow-up
(*Geoscience Frontiers* 12:101245, 2021, "A pulse of the Earth") combined 9
tetrapod events, 12 marine events and other geological events. It found the
same 27.5-Myr period, with pulses "centered at about 7, 35, 61, 95, 124, 140,
186, 196, 213, and 255 Ma", the most recent about 7 Ma. A comment by
Omerbashich (*Hist. Biol.* 34:188) disputed the spectral method; the authors
replied (34:212).

Both papers analyse **a list of event ages** compiled from the literature. We
found no test of the claim on an occurrence-based, genus-level extinction-rate
time series of land tetrapods. Searches run on 2026-11-17 included "periodicity
tetrapod extinction Paleobiology Database", "Rampino 27.5 non-marine tetrapod
critique" and the Historical Biology comment and reply.

## Ephemeris under test (fixed in advance)

**Primary:** extinction maxima at **t = 7.0 + 27.5 n Ma** (RCZ2021: period
27.5 Myr; most recent pulse ~7 Ma). The circular mean of their ten pulse ages
modulo 27.5 is ~7.5 Ma (R ~0.35; computed in `test_tetra.py`), consistent
with 7.0.

**Secondary:** the marine ephemeris of PREREG-1, 9.5 + 27.16 n Ma (Melott &
Bambach 2014). If one clock drives extinction on land and at sea, land
tetrapods should follow it too.

## Part A (no PBDB data; computed now)

Distances of RCZ2021's nine tetrapod events from the nearest predicted maximum
(GTS2020 ages): Late Miocene 7.25 Ma: 0.25 Myr; Eocene-Oligocene 33.9: 0.6;
K-Pg 66.0: 4.0; Jurassic-Cretaceous 143.1: 1.4; Triassic-Jurassic 201.4: 1.9;
mid-Norian 215: 12.0 (almost exactly anti-phase); Carnian-Norian 227.3: 0.3;
Permian-Triassic 251.9: 2.6; end-Guadalupian 259.8: 5.3. That is 6 of 9 within
+/-3 Myr against 21.8% expected by chance (binomial p ~0.005). This is **not
evidence**: the ephemeris was fitted to these events. The test below asks
whether the cycle also describes the parts of the record it was not built from.

## Data

Paleobiology Database, data1.2 `occs/list.csv`, two downloads: `base_name=
Tetrapoda^Mammalia` and `base_name=Mammalia`, each with `envtype=terr` (all
terrestrial environments), `idreso=genus`, `idqual=genus_certain`,
`pres=regular` (no ichno- or form taxa), `show=class`. Files are cached
with md5 sums in `provenance_tetra.json`.

Binning: GTS2020 international stages (the `cycle.py` table). An occurrence is
used only if its age range lies within one stage (tolerance 0.1 Myr). Many
terrestrial ages come from land-mammal ages that cross stage boundaries; those
records are dropped, and the counts are reported.

Rate: Alroy (2008) three-timer extinction rate mu at each stage **top**,
for genera keyed as phylum|genus. Window: stage tops 2.5-300 Ma (57 stages,
Gzhelian to Piacenzian), stages with a defined mu only.

## Tests (all one-sided, 20,000 surrogates, seed 20261119)

- **Sanity gate T1:** at least 40 window stages have a defined mu, AND the
  Maastrichtian (end-Cretaceous) mu is in the top quartile of the window.
  If T1 fails, the result is reported as DATA_LIMITED and nothing is claimed.
- **H1 (primary):** Pearson r between linearly detrended mu and
  cos(2 pi (t - 7.0) / 27.5). p_ar1 comes from AR(1) surrogates with the
  series' own lag-1 autocorrelation. p_phase comes from the same data against
  randomly phased templates.
- **H2 (leave events out, the out-of-sample core):** the H1 test after
  removing the nine stages that hold RCZ2021's tetrapod events (Tortonian,
  Priabonian, Maastrichtian, Tithonian, Rhaetian, Norian, Carnian,
  Changhsingian, Capitanian).
- **H3:** normal-score (rank) version of H1.
- **H4:** Lomb-Scargle periodogram of detrended mu, 10-100 Myr; band power
  26-29 Myr vs AR(1).
- **H5:** event test: local maxima of detrended mu above the 75th percentile,
  hits within +/-3 Myr of a predicted maximum, binomial vs 6/27.5.
- **H6:** the marine ephemeris (9.5 + 27.16 n) on the same series.
- **H7 (artifact check):** stage durations against the same template.
  Lieberman & Melott (2007) warned that stage lengths can carry a ~27-Myr
  pattern.
- **Secondary:** non-mammal tetrapods only; window 2.5-260 Ma (the "last 260
  Myr" of the 2021 paper); Foote's boundary-crosser rate q; timescale jitter of
  stage ages (0.25 / 0.75 / 1.0 Myr sd for Cenozoic / Mesozoic / Paleozoic).

## Decision rule

- DATA_LIMITED if T1 fails.
- SUPPORTED_BEYOND_EVENTS if H1 has r > 0 and p_ar1 < 0.05 **and** H2 has
  r > 0 and p_ar1 < 0.05.
- SUPPORTED_ONLY_BY_EVENTS if H1 holds but H2 does not: the cycle in the full
  series is carried by the very events it was fitted to.
- NOT_SUPPORTED otherwise.

We will say that "the land-tetrapod record independently supports a 27.5-Myr
clock" only for SUPPORTED_BEYOND_EVENTS, and only if H7 shows no significant
duration artifact (p_ar1 >= 0.05).

## Our probabilities, before seeing data

| ID | Prediction | Our probability |
|---|---|---|
| T1 | Sanity gate passes | 70% |
| T2 | H1 r > 0 | 65% |
| T3 | H1 r > 0 and p_ar1 < 0.05 | 25% |
| T4 | H2 p_ar1 >= 0.05 (no support once the events are removed) | 85% |
| T5 | H1 p_phase < 0.05 (the published phase beats random phases) | 20% |
| T6 | The global periodogram peak (10-100 Myr) is outside 26-29 Myr | 80% |
| T7 | H4 band p_ar1 >= 0.05 | 75% |
| T8 | H6 (marine ephemeris) r > 0 | 60% |
| T9 | H7 duration artifact p_ar1 >= 0.05 | 80% |
| T10 | Non-mammal series has the same sign of r as H1 | 65% |
| T11 | H5 event test p_binom >= 0.05 | 70% |

Our overall expectation: the large events will line up with the ephemeris,
but the rest of the land record will not (SUPPORTED_ONLY_BY_EVENTS or
NOT_SUPPORTED). We would treat SUPPORTED_BEYOND_EVENTS as a genuine surprise.

## What would make this remarkable

If land tetrapods follow the same fixed-phase 27.5-Myr ephemeris in the
stages that were not used to fit it, that would be independent evidence, from
a separate realm of life, for a pacemaker of extinction. Proposed pacemakers
include the Sun's passage through the galactic disc. A clean negative would
show that the land-tetrapod cycle lives only in a hand-picked event list.
Both outcomes will be published here, with all data and code.

## Limits we already know

PBDB land-tetrapod sampling is uneven: Mesozoic data come mostly from a few
basins, and Cenozoic data are dominated by North American and European mammals.
Three-timer rates are noisy when few genera span three stages. Stage-level
binning drops many occurrences dated by regional land-mammal ages. These limits
apply to any outcome and will be reported with the result.

Code Apache-2.0; text CC BY 4.0. MV2, Munim, Inc., mv2@munimtech.com
