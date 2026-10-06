# CLAUDE.md — Dual Function Recovery (HRS Multistate Analysis)

## Project overview
This project analyzes transitions into and out of physical and cognitive impairment, plus mortality, among older adults in the Health and Retirement Study (HRS), 1998–2022. The analysis uses the Bayesian discrete-time multistate life table (MSLT) approach of Lynch and Zang. The paper is co-authored with Samuel Nemeth. Conceptual framing and writing are handled outside this repo. This repo covers data construction, estimation, and output.

### State space (5 states)
1. Dual functional (no physical or cognitive impairment)
2. Physical impairment only
3. Cognitive impairment only
4. Both impaired
5. Dead (absorbing)

States 2–4 are transient. Recovery transitions (e.g., 4 → 2, 2 → 1) are substantively central and must not be ruled out by model structure.

## Data
All data live in `data/raw/` and are **read-only**. Never modify, overwrite, or move raw files.

- **RAND HRS Longitudinal File** (through 2022)
- **Langa-Weir Classification of Cognitive Function** (through 2022)
- The exact filenames and versions are listed in `docs/data-inventory.md`. Create that file on first run by inspecting `data/raw/`.

Data rules:
- Never commit data files. `data/` is gitignored, except for `data/README.md`.
- Do not print individual-level records beyond what is needed for debugging. Prefer summaries, counts, and `glimpse()`-style structure checks.
- Derived datasets go in `data/derived/` as `.rds`, and each one is produced by a numbered script.

## Key measurement decisions
Confirmed decisions are marked ✅ (details and rationale in `docs/decisions.md`). Open decisions are marked ❓. **Do not resolve an ❓ item on your own.** Propose options with tradeoffs and ask. Log every resolved decision with its date and rationale in `docs/decisions.md`.

- ✅ Observation window: 1998–2022 waves (RAND waves 4–16). 1998 is the baseline because it is the first wave with all original cohorts combined.
- ✅ Mortality: `r#iwstat == 5` (died since last interview) is the primary indicator. RAND `radyear`/`radmonth` give death timing when needed.
- ✅ Cognitive status: Langa-Weir classification (normal / CIND / dementia). This classification covers proxy respondents.
- ✅ Cognitive impairment cut: dementia only (Langa-Weir `cogfunction == 3`); CIND + dementia is a sensitivity analysis.
- ✅ Physical impairment: ≥1 difficulty on the RAND 5-item ADL summary (`r#adl5a > 0`); no IADL or mobility.
- ✅ Analytic sample: age 50+ per interval; all cohorts contributing ≥1 interval (entry at first interview); all household members; foreign-born kept.
- ✅ Missing waves/attrition: consecutive-wave intervals only; censor at a missed wave, re-enter on return; deaths after a gap ignored; intervals with missing state dropped.
- ✅ Covariates: age (linear + quadratic), sex, race/ethnicity, education, birth cohort. Analysis will be stratified by sex; natural-spline age possible later.
- ✅ Transition interval: treat all intervals as 2 years; interval length as covariate is a sensitivity analysis.

## Method: Lynch & Zang Bayesian MSLT
The source article is in `refs/`. **Before writing any estimation code, read it and write `docs/methods-notes.md`**. That file should summarize:
- the transition model specification,
- the estimation approach (sampler, priors, convergence checks),
- how posterior draws are turned into life tables and state expectancies,
- what quantities the method produces and how uncertainty is summarized.

The notes should also flag where our 5-state space or our data differ from the article's illustration. Treat the article as authoritative. Where it is ambiguous, say so instead of guessing.

## Code conventions
- R only. Use tidyverse style, `here::here()` for paths, and no `setwd()`.
- Use `renv` for package management. Run `renv::snapshot()` after adding packages.
- Scripts in `R/` are numbered by pipeline order, e.g. `01_import.R`, `02_construct_states.R`, `03_transitions.R`, `04_estimate.R`, `05_life_tables.R`, `06_figures.R`.
- Reusable functions go in `R/functions/` and are sourced at the top of scripts.
- Set and record seeds for all stochastic steps.
- Long model runs save posterior draws to `output/models/` so downstream scripts don't re-estimate.
- Reports and exploratory write-ups are written in Quarto (`.qmd`) in `reports/`.
- Figures use ggplot2 and are saved to `output/figures/` as PDF and PNG.

## Folder structure
```
.
├── CLAUDE.md
├── README.md
├── .gitignore
├── renv/ , renv.lock
├── data/
│   ├── README.md        # describes files; tracked
│   ├── raw/             # downloaded HRS files; read-only, untracked
│   └── derived/         # analysis datasets (.rds); untracked
├── refs/                # Lynch & Zang article, codebooks; PDFs untracked
├── R/
│   ├── functions/
│   └── 01_... 06_...R
├── output/
│   ├── models/          # posterior draws; untracked
│   ├── tables/
│   └── figures/
├── reports/             # Quarto exploratory and diagnostic reports
└── docs/
    ├── data-inventory.md
    ├── decisions.md
    └── methods-notes.md
```

## Project status and handoff (updated 2026-10-05)
Read `docs/decisions.md` (all resolved decisions with rationale), `docs/methods-notes.md` (article summary, `bayesmlogit` comparison §7, sampler pilot §8) and `docs/data-inventory.md` before continuing.

### Setting up on a new machine
1. `git clone https://github.com/sbauldry/dft`; open R in the project and run `renv::restore()` (R 4.6.x was used; `BayesLogit` and `bayesmlogit` compile/install from CRAN).
2. Copy the raw files (not in git) into `data/raw/`: `randhrs1992_2022v1.dta` (RAND HRS 1992–2022 V1, 1.7 GB) and `cogfinalimp_9522wide.dta` (Langa-Weir imputed cognition, wide). Put the Lynch & Zang PDF in `refs/` (optional).
3. Rebuild derived data (all of `data/derived/` is untracked): run `R/01_import.R`, `R/02_construct_states.R`, `R/03_transitions.R` in order (a few minutes). Expected: 205,974 intervals, 35,549 persons (men 86,095 intervals; women 119,879). Optional check: `R/03a_prevalence_check.R` (dual functional prevalence by age and sex; results looked right to the user).

### Pipeline built so far
- `R/01_import.R` reads selected RAND + Langa-Weir variables; links Langa-Weir to RAND via `hhidpn = as.numeric(paste0(hhid, pn))` (42,882 of 42,890 LW ids match).
- `R/02_construct_states.R`: 4 living states (1 functional, 2 ADL only, 3 dementia only, 4 both) from `adl5a > 0` and `cogfunction == 3`; `state_cind` is the CIND+dementia sensitivity version.
- `R/03_transitions.R`: person-intervals under the logged rules. Transition counts (men + women), start rows × end columns 1–5 (5 = dead): 1: 139418/11945/3225/1635/5914; 2: 8402/14517/464/1878/3924; 3: 1767/378/1676/1182/967; 4: 355/848/392/3493/3594. All 20 cells populated (smallest 355).
- Data flags: recovery out of dementia-only (state 3→1 about 30% per 2 years) is probably Langa-Weir classification noise near the cutoff; cognitive-only state is small (~3% of person-waves).
- `R/functions/constants.R`, `design.R` (outcome code `(start-1)*5+end`, ref = 1→1; predictors intercept, age centred 70 per decade, age², race/ethnicity 3 dummies, education 3 dummies, birth cohort per decade centred at sample mean), `pg_mlogit.R` (Pólya-Gamma Gibbs), `mh_mlogit.R` (independence Metropolis-Hastings, t proposal at posterior mode), `life_table.R` (life table from one draw; matches `bayesmlogit::mlifeTable` exactly on shared draws).
- `R/04a_validate_simulation.R` (sampler recovers simulated truth); `R/04b_compare_bayesmlogit.R` (comparison with the CRAN package `bayesmlogit`; takes ~15 min).
- `R/04_estimate.R`: independence-sampler estimation, `Rscript R/04_estimate.R <men|women> [prior_sd] [n_iter] [burn] [thin] [df]`; refuses to overwrite existing output.
- `R/04c_prior_sensitivity.R`: compares posterior means/SDs across priors sd 2.5/5/10 (tables in `output/tables/prior_sensitivity_<sex>.csv`; results in `docs/methods-notes.md` §9).
- `R/05_life_tables.R`: `Rscript R/05_life_tables.R <men|women> [prior_sd] [n_draws] [cores]`; life tables per draw (ages 50–110+, population and status-based radix), summaries to `output/tables/lifetable_summary_<sex>_psd<sd>.csv`, draw-level values at ages 50–90 to `output/models/lifetable_<sex>_psd<sd>.rds`; ~6 s per sex; refuses to overwrite.
- `R/06_figures.R`: figures from the life table summaries (`Rscript R/06_figures.R [prior_sd]`).
- `reports/diagnostics.qmd`: data, sampler convergence, prior sensitivity and life-table checks.

### Confirmed estimation plan (details in `docs/decisions.md`)
Separate models for men and women (sex not a predictor); 20 transition outcomes; priors N(0, 5²) (sensitivity N(0, 2.5²) and N(0, 10²)); life tables for ages 50–110+ in 2-year steps at birth cohort = sample mean and race/education dummies at sample proportions, population-based and status-based radix; report 95% and 84% credible intervals; unweighted. 

### ✅ Sampler choice (resolved 2026-10-05; see `docs/decisions.md`, `docs/methods-notes.md` §8)
Independence Metropolis-Hastings (t(30) proposal at the posterior mode), 4 chains × 100,000 proposals, burn-in 10,000, thin 10 (36,000 draws per sex); Pólya-Gamma Gibbs mixed poorly on full data and is kept only as a cross-check. Main runs (prior N(0, 5²)) are done and saved untracked as `output/models/draws_men_psd5.rds` and `draws_women_psd5.rds` (logs alongside; ~1–1.5 h each). Men: acceptance ~0.17, max R-hat 1.002, bulk ESS min 3,858. Women: acceptance ~0.24, max R-hat 1.001, bulk ESS min 13,936. Prior sensitivity runs (sd 2.5 and 10) are also done (all R-hat ≤ 1.002); posterior means shift ≤0.51 posterior SD (largest: college terms on rare transitions, esp. 4→3), SDs unchanged, so sd 5 stays the main prior. Never overwrite existing runs without asking.

### Status and next steps (2026-10-06)
Done: data rebuilt and verified (205,974 intervals, 35,549 persons); all six model runs (men/women × prior sd 2.5/5/10, all R-hat ≤ 1.002); life tables for all six (sd 5 total life expectancy at 50, population radix: men 27.8, women 31.5 years; prior choice changes life-table quantities by ≤0.16 years, see `docs/methods-notes.md` §9); figures (`R/06_figures.R` → `output/figures/06_*.pdf/png`: state expectancies, total expectancy, % of remaining life, status-based radix at 50; ages shown to 100); Quarto diagnostics report (`reports/diagnostics.qmd`; render with RStudio's bundled Quarto, `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render reports/diagnostics.qmd`; the rendered HTML is not tracked).
Next steps (in order; user decided 2026-10-06 that sensitivity analyses come after figures and report, which are done):
1. **CIND + dementia sensitivity.** Use `state_cind` / `start_state_cind` (already built in `R/02`/`R/03`; transition counts in `output/tables/03_transition_counts_cind.csv`). Needs `R/03_transitions.R` to write an interval file for the CIND version (check what it saves), `R/04_estimate.R` and `R/05_life_tables.R` to take an analysis tag (e.g. `main`/`cind`) so output names don't collide with existing `draws_<sex>_psd5.rds`, then run both sexes at prior sd 5 (~1–1.5 h each; run in background, do not edit the script while running), life tables, and compare with the main results (state expectancies, % of remaining life). Expect much larger cognitive states (cognitive-only ~12–13k intervals vs ~3k).
2. **Interval-length sensitivity.** Add interval length (`interval_months`, centred at 24) as a covariate; fix at 24 months for life tables (decision logged 2026-10-04). Same tagging approach as above; add a design variant in `R/functions/design.R` rather than changing the main design.
3. For each sensitivity analysis: propose the exact plan, log the decision in `docs/decisions.md` and results in `docs/methods-notes.md` (§9 onward), add a section to `reports/diagnostics.qmd` and figures to `R/06_figures.R` if useful, then commit.
4. Optional/later: Gibbs cross-check on the final model, natural-spline age, report rendering to share (HTML not tracked).
Housekeeping: raw files were found in `data/` and copied to `data/raw/` (originals left in `data/`, untracked). GitHub push needs credentials set up on the machine in use; push status should be checked with `git status -sb`.

### Tooling notes
- `bayesmlogit` (CRAN 1.0.1, Zang, Zhang & Lynch) uses a flat prior, one chain, per-observation R loops (very slow at this scale) and only varies a single `age` column, so it cannot take our age + age² model. We use our own code and the package as a validation reference.
- Background jobs: long runs should be started with `run_in_background` and a long timeout, and the script file must not be edited while it is running.

## Working style
- Be concise. Report what was done and what was found. Skip narration.
- Before coding a pipeline stage, state the plan in a few lines and wait for approval if it involves an ❓ decision.
- Check results as you go. Verify sample counts, state distributions by wave, and transition tables before modeling. Flag implausible patterns, such as transitions out of death or large unexplained sample drops.
- Never delete or overwrite outputs from a previous model run without asking.
