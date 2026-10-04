# Data inventory

Created 2026-10-04 by inspecting the files (headers and a few columns only; no records printed).

## Location
On 2026-10-04, at the user's request, both `.dta` files were moved from the repo root to `data/raw/` (made read-only), and the Lynch & Zang PDF to `refs/`. Contents unchanged.

## Files

| File | Source | Size | Modified | Rows | Cols |
|---|---|---|---|---|---|
| `randhrs1992_2022v1.dta` | RAND HRS Longitudinal File 2022 (V1), Stata | 1.74 GB | 2025-05-08 | 45,234 (one per respondent; `hhidpn` unique) | 19,880 |
| `cogfinalimp_9522wide.dta` | Langa-Weir classification of cognitive function, imputed, wide, 1995–2022 | 12.8 MB | 2026-09-09 | 42,890 | 286 |

### RAND HRS (`randhrs1992_2022v1.dta`)
- ID: `hhidpn` (numeric). Wide format, `r{w}`-prefixed variables for waves 1–16 (waves 4–16 = 1998–2022).
- Variables confirmed present: `r1iwstat`–`r16iwstat`, `radyear`, `radmonth`, `rabyear`, `ragender`, `raracem`, `rahispan`, `raeduc`, `rabplace`, `hacohort`.
- `r#iwstat` counts, waves 4–16 (codes: 0 inapplicable, 1 resp alive, 4 nonresp alive, 5 nonresp died this wave, 6 nonresp died prior wave, 7 nonresp dropped from sample):

| wave | 0 | 1 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| 4 (1998) | 18774 | 21384 | 2173 | 1345 | 1312 | 246 |
| 5 (2000) | 18628 | 19578 | 2486 | 1441 | 2656 | 445 |
| 6 (2002) | 18468 | 18165 | 2250 | 1571 | 4097 | 683 |
| 7 (2004) | 14923 | 20129 | 2330 | 1297 | 5641 | 914 |
| 8 (2006) | 14765 | 18469 | 2196 | 1384 | 6925 | 1495 |
| 9 (2008) | 14640 | 17217 | 2144 | 1298 | 8305 | 1630 |
| 10 (2010) | 8025 | 22034 | 2240 | 1609 | 9603 | 1723 |
| 11 (2012) | 7850 | 20554 | 2265 | 1202 | 11212 | 2151 |
| 12 (2014) | 7703 | 18747 | 2551 | 1343 | 12414 | 2476 |
| 13 (2016) | 3012 | 20912 | 3340 | 1480 | 13757 | 2733 |
| 14 (2018) | 2886 | 17146 | 5391 | 1221 | 15237 | 3353 |
| 15 (2020) | 2790 | 15723 | 5049 | 1447 | 16458 | 3767 |
| 16 (2022) | — | 15856 | 5748 | 1382 | 17905 | 4343 |

  Code 4 (nonresponse, alive) is sizable in 2018–2022 (5–6k per wave); relevant to the missing-wave ❓ decision.
- `hacohort` distribution (codes 0–8): 108, 8337, 2431, 13659, 2831, 5081, 5387, 4672, 2728.

### Langa-Weir (`cogfinalimp_9522wide.dta`)
- IDs: `hhid` (character, `%6s`) and `pn` (character, `%3s`). **No `hhidpn`**; linking to RAND requires constructing it (expected `hhid` + `pn` as numeric; verify in `01_import.R`).
- Wide, suffixed by year. Years present: 1995, 1996, 1998, 2000, …, 2022.
- Per-year variables (stem + year): `cogfunction` (1 = Normal, 2 = CIND, 3 = Demented, per label on `cogfunction1998`), `cogtot27_imp`, `proxy`, `interview`, plus component scores (`imrc`, `dlrc`, `ser7`, `bwc20`, `memoryp`, `numiadl`, `prxyscore`, and `f*` flags). Time-invariant: `firstiw`, `study_cohort`.
- Note: the file name and variables say "imp" (imputed). Confirm the version/source with the user.

## Open items
1. Exact versions/download dates of both files are not recorded in the files; add if known.
