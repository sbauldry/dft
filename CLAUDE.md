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
- ✅ Covariates: age (linear + quadratic), sex, race/ethnicity, education, birth cohort. Possible sex-stratified models and natural-spline age later.
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

## Working style
- Be concise. Report what was done and what was found. Skip narration.
- Before coding a pipeline stage, state the plan in a few lines and wait for approval if it involves an ❓ decision.
- Check results as you go. Verify sample counts, state distributions by wave, and transition tables before modeling. Flag implausible patterns, such as transitions out of death or large unexplained sample drops.
- Never delete or overwrite outputs from a previous model run without asking.
