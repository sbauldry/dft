# data/

Not tracked by git except this file.

- `raw/` — downloaded HRS files (RAND HRS Longitudinal File, Langa-Weir classification). Read-only. See `docs/data-inventory.md` for file names and versions.
- `derived/` — analysis datasets (`.rds`), each produced by a numbered script in `R/`.

## Derived files (`data/derived/`)
| File | Script | Contents |
|---|---|---|
| `hrs_wide.rds` | `R/01_import.R` | Selected RAND HRS variables (waves 4–16) linked to Langa-Weir `cogfunction` and `proxy` by `hhidpn` (built as `hhid` + `pn`) |
| `person.rds` | `R/02_construct_states.R` | One row per person: cohort, birth/death dates, sex, race, Hispanic, education, birthplace |
| `person_wave.rds` | `R/02_construct_states.R` | One row per person-wave: interview status/date, age, ADL count, cognition class, 4-state classification (`state`; `state_cind` for the CIND + dementia sensitivity) |
| `intervals.rds` | `R/03_transitions.R` | Person-interval analysis file: start/end state (5 = dead), age, covariates, interval length |
