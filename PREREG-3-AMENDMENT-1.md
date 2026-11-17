# PREREG-3, Amendment 1: corrected PBDB environment code (before any data)

Dated 2026-11-17 (about 12:30 UTC), by MV2 of Munim, Inc. Made with AI agents
working in his name.

## What changed

PREREG-3 registered the land-tetrapod downloads with `envtype=terr`, meaning
"all terrestrial environments". A code review before the job ran found that
the Paleobiology Database API (data1.2) **rejects** that value. A test request
(`base_name=Mammalia&envtype=terr&...&limit=3`) returned the warning
`bad value 'terr' for parameter 'envtype'` and **zero records**. The job would
therefore have stopped with `data_incomplete` and produced no result.

The API's documented code for "Any terrestrial environment" is
`envtype=terrestrial`. A test request with that value returned the normal CSV
header with no warning. The same was true for `base_name=Tetrapoda^Mammalia`.

The only change is:

| | Registered | Now used |
|---|---|---|
| PBDB `envtype` | `terr` | `terrestrial` |

This is what PREREG-3 already said it meant ("all terrestrial environments").
Nothing else changes: the taxa, the filters (`idreso=genus`,
`idqual=genus_certain`, `pres=regular`, `show=class`), the binning, the
ephemeris (7.0 + 27.5 n Ma), the tests, the decision rule and our
probabilities T1-T11 are all as registered.

## What we had seen

No tetrapod occurrence data had been downloaded or analysed. The two test
requests used `limit=3` and `limit=5`. They returned a handful of rows
(for example *Brabovus*, *Saturninia*, *Chaoyangsaurus*, *Caudipteryx*). We
looked at these only to confirm the CSV format and computed nothing from them.

## Files

- `tetra.py`: `COMMON` now uses `envtype=terrestrial`. A note in the module
  docstring points here.
- `test_tetra.py`: `test_queries` now checks for `envtype=terrestrial&` and
  rejects `envtype=terr&`.
- `job_tetra.py`: records this amendment in `provenance_tetra.json` and in
  the `prereg` field of `results_tetra.json`. It also uses per-process
  temporary download names, as the other jobs do. The analysis is unchanged.

The marine tests (PREREG-1 and PREREG-2) use `envtype=marine`, which the API
accepts. They are not affected.

Code Apache-2.0; text CC BY 4.0. MV2, Munim, Inc., mv2@munimtech.com
