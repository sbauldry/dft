# Methods notes: Lynch & Zang (2022) Bayesian MSLT

Source: Lynch, S. M., & Zang, E. (2022). Bayesian multistate life table methods for large and complex state spaces. *Sociological Methodology*, 52(2), 254–286. Page references are to the journal pagination. Status: **draft, 2026-10-04**. The online supplement (trace plots, posterior draws, simulation, R code) was not available; anything it contains is marked as unknown.

## 1. Transition model specification

**Process.** Discrete-time, first-order Markov, finite state space, time-inhomogeneous (transition matrix P(a) varies with age) (pp. 260–262). Death is absorbing and reachable from every state. Waves are assumed evenly spaced (interval k; k = 2 in the illustration).

**Data layout** (p. 262–263). Long person-interval file. Each row has starting state at t, ending state at t+k, age at t, and covariates. A respondent contributes (waves − 1) rows if surviving; deaths end the record; attrition rows are censored. The ending state of one row is the starting state of the next.

**The key change from Lynch & Brown (2005)** (p. 263–264). The outcome of the multinomial logit is the *transition* (the starting-state × ending-state pair), not the ending state conditional on a starting-state covariate.
- Impossible transitions (structural zeros) are simply not outcomes.
- All covariates, including age, get a separate coefficient vector for each transition, so the age effect varies by transition without interactions.
- Model: P(y_i = j) = exp(X_i β_j) / Σ_s exp(X_i β_s), with one transition omitted as the reference (β is m × (J−1)).
- Likelihood is multinomial; posterior ∝ prior × likelihood.

**Covariates.** Any number. The illustration uses age (decimal), sex, race, Hispanic, birth cohort (year − 1900), marital status, education, birth region, current region. Time-varying covariates (current region) are observed per interval.

**Illustration's rare-cell handling** (p. 271). The A→DC transition had only 10 observations and was recoded to A→DCA because of convergence problems. Result: 42-dimensional outcome after dropping that transition and the reference.

## 2. Estimation

**Sampler choice** (pp. 264–266). The article gives rough guidance by number of possible transitions J:

| J | Sampler that works |
|---|---|
| < ~10 | Probit Gibbs (data augmentation) or logit independence sampler, both fine |
| ~10 to ~25 | Logit independence sampler better than probit Gibbs |
| > ~25 | Both become hard or infeasible; use Pólya-Gamma Gibbs for logit (Polson, Scott & Windle 2013) |

- Independence sampler: proposal is multivariate normal at the ML estimate with the ML covariance matrix. It works well when (1) n is large with few rare transitions, (2) priors are relatively noninformative, and (3) outcome dimension is small to medium (fewer than ~15–20 transitions).
- Pólya-Gamma Gibbs: used in the illustration. The authors say it converges and mixes orders of magnitude faster than the independence sampler. They adapted functions from the (unsupported) `BayesLogit` R package.

**Run settings in the illustration** (p. 272). Two chains from random starting values, 2,500 draws each, first 500 dropped, thinned to every 4th, giving 1,000 posterior draws total. Convergence and mixing were "monitored" with trace plots in the supplement.

**Priors.** Priors on β are allowed, including strong priors to force structural zeros. The article says only that they should be "relatively noninformative" for the independence sampler. **Ambiguous / not stated:** the actual prior family and variances used in the illustration.

**Convergence diagnostics.** **Ambiguous:** beyond "monitoring convergence and mixing across the two chains" and trace plots, no formal statistic (R-hat, ESS) is named.

**Survey design.** The illustration ignores weights and clustering (assumes simple random sampling). The article says a weighted bootstrap draw at each Gibbs iteration (Gunawan et al. 2017) is the preferred way to incorporate weights and, in unreported results, changed estimates only slightly.

## 3. From posterior draws to life tables

For each posterior draw g (pp. 266–268):

1. **Covariate profile.** Fix covariate values. Build the matrix Z (N age groups × m), incrementing age by k per row, from the youngest age (50) to an open-ended last group (110+). Illustration: 31 age groups, 2-year steps.
2. **Transition probabilities.** Compute Zβ(g). Convert to probabilities by row with the multinomial logit (reference transition = 1 − sum of others). These are *joint* probabilities of (start state i, end state j).
3. **Assemble matrices.** Place each row's probabilities into a d × d matrix by start/end state. Structural zeros stay zero. Append a final row of zeros with a trailing 1 (death row).
4. **Normalize rows** so each living-state row sums to 1: p(j|i) = p(i,j)/Σ_j p(i,j). This yields a right-stochastic P(a) for each age group.
5. **Radix.** Population-based: from the row sums of the *unnormalized* probabilities at the youngest age. Status-based: set all of the radix in one chosen state.
6. **Life table.**
   - l(a+1) = l(a) P(a)
   - L(a) = 0.5 k [l(a) + l(a+1)] (linear assumption)
   - Open interval: L(N) = k l(N) [I − P(N)]⁻¹, with the death row and column removed. This assumes constant transition *probabilities* beyond the last age (geometric waiting time), preferred over converting to rates.
   - T(a) = Σ_{i=a}^{N} L(i)
   - State expectancy e(a) = T(a) / Σ l(a)
7. Repeat for all G draws to get a posterior distribution of every life table quantity.

**Cohort vs period.** If birth cohort is a covariate (the authors recommend it), it must be fixed to produce a table. The result is a cohort life table informed by other cohorts, so it assumes some stationarity.

**Footnote 1 (p. 282).** Model-implied column marginals at the end of one interval need not match row marginals at the start of the next. The authors argue this is unnecessary to enforce (data from multiple cohorts and attrition don't satisfy it either) and support this with a simulation in the supplement (not seen).

## 4. Outputs and uncertainty summaries

- **Quantities:** state expectancies (years in each state) at any age, total life expectancy, sums of states (e.g., any-diabetes = D+DA+DC+DCA), percentage of remaining life in a state (computed draw by draw as years in state / total years), population-based and status-based variants.
- **Uncertainty:** posterior mean with **84% credible intervals** from sorted draws (quantiles). The 84% choice is justified because non-overlap of two 84% intervals roughly corresponds to a two-sided p < .05 test for equal SEs (Payton et al. 2003).
- **Comparisons:** posterior probability that one group's outcome is worse than another's = share of draws where the difference has a given sign (Table 4). No interval-overlap test needed.
- **Uncertainty scope:** parameter (sampling) uncertainty from the model only; no uncertainty about interval timing, measurement, or model form.

## 5. Where our setting differs from the article's illustration

| Issue | Article | This project | Consequence / flag |
|---|---|---|---|
| State space | 8 living states + death, quasi-absorbing (no recovery) | 4 living + death, **all living-to-living transitions allowed** (recovery is central) | Fewer states but no structural zeros; see below |
| Possible transitions J | 43 (42 + reference), after structural zeros | From states 1–4 to any of 5 states = **20**, minus reference = 19 | Falls in the 10–25 range, where the article says the independence sampler beats probit Gibbs but warns about rare transitions. Choice of sampler is open (see §6) |
| Rare transitions | 10 obs. A→DC recoded | Likely rare: direct 1↔4, 4→1, 2↔3 jumps | Need transition counts before deciding on recoding vs. strong priors (a ❓ decision) |
| Interpretation of states | Cumulative (conditions once acquired) | Reversible functional states | Structural-zero priors not needed; recoveries are estimated directly |
| Observation window | HRS 1998–2014 | 1998–2022 (waves 4–16) | Includes the 2020 (COVID) wave and 2022; longer span strains stationarity and the cohort/period interpretation |
| Sample | 50+, in 1998 or 2004/2010 cohorts, one per household, US-born, never lived abroad; 17,686 persons | Not decided (❓) | Article's restrictions are one possible template |
| Missing data | Dropped missing transitions (4.5%); attrition censored | Not decided (❓) | RAND `iwstat` 4 (nonresponse, alive) is common in 2018–2022 |
| Interval timing | Assumes equal; ignores actual dates (simulations suggest little effect if unrelated to covariates) | Not decided (❓) | Article does not offer a variable-interval method |
| Weights | Ignored | Not decided | Article offers weighted bootstrap within MCMC as an option |
| Cognition measure | Not used | Langa-Weir, includes proxies | Proxy-based classification may differ in error structure |
| Age effect | "age variable" included | Functional form unspecified in article | **Ambiguous**; linear age over 50–100+ is questionable for 4 states plus death |
| Age grid | 50 to 110+, 2-year steps | To be decided with sample age floor | Open-ended last interval matters for state expectancies at old ages |

**Note on state-space difference.** With no quasi-absorbing structure, the article's main reason for modelling transitions (avoiding separation from structural zeros) does not apply. The transition-as-outcome formulation is still the method's specification, and the authors report it gives nearly identical probabilities (third decimal) to a starting-state-as-covariate model (p. 267), differing only in that covariate effects vary by transition. Whether to use it as-is or to use an alternative parameterization is a modeling decision for the user.

## 6. Ambiguities and open questions

1. **Priors.** Family and scale are not stated. We will need to specify and justify; propose weakly informative normal priors and run a sensitivity check.
2. **Sampler for J = 19.** Article guidance is a range, not a rule. Options: independence sampler (simpler; risk of poor acceptance with rare transitions) vs. Pólya-Gamma Gibbs (robust; needs custom or ported code since the original package is unsupported). Test on our data and record acceptance rates/diagnostics.
3. **Convergence diagnostics.** Article shows only trace plots; we should add R-hat and effective sample size.
4. **Functional form of age and interactions** (e.g., splines; age × sex) is not specified.
5. **Radix for population-based tables.** The article only says "row sums of the unnormalized probabilities at age 0." *Resolved by the `bayesmlogit` source (2026-10-04):* `mlifeTable()` takes row sums of the unnormalized matrix at the first age **including the transition-to-death column**, sets the dead entry to 0, and uses that as the radix (it sums to 1 over living starts because the 20 joint probabilities sum to 1). Our `R/functions/life_table.R` follows this and reproduces `mlifeTable()` exactly on shared draws.
6. **Draws and thinning.** 1,000 final draws from 2 chains; we may choose differently. No guidance on how many are needed for stable interval endpoints.
7. **Treatment of death within an interval.** Linear approximation L = 0.5k[l(a)+l(a+1)] applies to all states including those who die mid-interval; no separate death-timing adjustment (e.g., our `radyear`/`radmonth`) is described.
8. **Supplement and R functions.** Article states R functions and instructions are available from the authors; obtain them to check against these notes.

## 7. The `bayesmlogit` R package (checked 2026-10-04)

CRAN package `bayesmlogit` 1.0.1 (Zang, Zhang & Lynch): `bayesmlogit()` (sampler), `mlifeTable()` (life tables), `CreateTrans()`, plus plotting/comparison helpers. Findings from reading the source and running it:
- **Algorithm.** Same Pólya-Gamma Gibbs sampler for the multinomial logit as in the article; reference outcome = last transition code.
- **Prior.** Flat/improper (prior precision fixed at 0); not user-settable. Single chain from zero starting values; no built-in convergence diagnostics beyond optional trace plots.
- **Speed.** Pólya-Gamma draws are made one observation at a time in R loops. On N = 2,000 and 20 outcomes it took 417 s for 1,500 iterations (~0.28 s/iteration); our vectorised sampler (`BayesLogit::rpg`) takes ~0.017 s/iteration on the same data and ~0.9 s/iteration on the full men's file (N = 86,095), where `bayesmlogit` would be roughly 50x slower per iteration.
- **Life tables.** `mlifeTable()` only varies a single column named `age`, so it cannot handle our age + age^2 specification. It uses the same radix, linear L(a) and geometric open interval as the article. Our implementation agrees with it to 5 decimals on shared draws.
- **Mixing.** On a stratified N = 2,000 test (100 intervals per transition) the package chain had bulk ESS ~6 of 1,000 draws (median); ours had ~55 of 4,000 across 4 chains (R-hat max 1.08). Differences between the two posterior means (~1 posterior SD) are consistent with Monte Carlo error from this poor mixing; both are within ~0.6 posterior SD of the ML estimates. Poor mixing of the Pólya-Gamma Gibbs sampler with many correlated transition outcomes is a concern for the full models (see diagnostics).

## 8. Sampler pilot on the full men's data (2026-10-05)

Model: 190 parameters (10 predictors x 19 outcomes), N = 86,095 intervals, prior N(0, 5^2), 4 chains each.

| Sampler | Run | Max R-hat | Bulk ESS median / min | Notes |
|---|---|---|---|---|
| Pólya-Gamma Gibbs (`R/functions/pg_mlogit.R`) | 1,500 iter, 500 burn-in, 1,000 kept/chain, ~0.95 s/iter, 25 min | 1.58 | 147 / 7 (of 4,000) | Slow mixing concentrated in rare transitions (4→3, 3→2, 2→3, 4→1), esp. education/race terms; lag-1 autocorrelation up to 0.97 |
| Independence MH, t(8) proposal at the posterior mode (`R/functions/mh_mlogit.R`) | 20,000 proposals, 2,000 burn-in, thin 10, 0.036 s/proposal, 15 min | 1.038 | 741 / 140 (of 7,200); tail ESS min 43 | Acceptance ~9%; mode by damped Newton in 26 s |

- The two samplers agree: posterior means differ by a median 0.06 posterior SD (max 0.56 SD, on the parameters where the Gibbs chain mixes worst); both are within 0.5 SD of the posterior mode.
- Proposal tails: acceptance by degrees of freedom was 8% (df 4), 9% (8), 15% (30), 24% (100), 31% (normal).
- The article's guidance (independence sampler fine when n is large, few rare transitions, <15–20 outcomes) is borderline for J = 19 but holds here because every transition cell has >=355 events.

## 9. Full estimation runs and prior sensitivity (2026-10-05)

Independence sampler (`R/04_estimate.R`): t(30) proposal at the posterior mode, 4 chains × 100,000 proposals, burn-in 10,000, thin 10 (36,000 draws per sex). Men N = 86,095 intervals; women N = 119,879. Draws are in `output/models/draws_<sex>_psd<sd>.rds` (untracked).

| Sex | Prior sd | Acceptance | Max R-hat | Bulk ESS median / min | Tail ESS min |
|---|---|---|---|---|---|
| Men | 2.5 | | 1.002 | 8,348 / 4,955 | 2,123 |
| Men | 5 | 0.165–0.169 | 1.002 | 7,226 / 3,858 | 1,395 |
| Men | 10 | | 1.002 | 6,135 / 3,272 | 1,205 |
| Women | 2.5 | | 1.001 | 18,731 / 14,921 | 10,643 |
| Women | 5 | 0.235–0.242 | 1.001 | 18,260 / 13,936 | 9,486 |
| Women | 10 | | 1.001 | 18,136 / 13,425 | 8,809 |

Prior sensitivity (`R/04c_prior_sensitivity.R`; per-parameter tables in `output/tables/prior_sensitivity_<sex>.csv`). Differences in posterior means relative to sd 5, in units of the sd 5 posterior SD:

| | Median abs | Max abs |
|---|---|---|
| Men, sd 2.5 | 0.022 | 0.51 |
| Men, sd 10 | 0.007 | 0.21 |
| Women, sd 2.5 | 0.025 | 0.20 |
| Women, sd 10 | 0.006 | 0.06 |

- Posterior SDs are unchanged (median ratio 0.998–1.000).
- The largest shifts are on the college-education terms for the rarest transitions, mainly 4→3 (outcome 18, 392 events): men -5.24 (sd 5), -4.73 (sd 2.5), -5.45 (sd 10), posterior SD ~1.0. Shifts for 3→2, 2→3 and 4→1 college terms are <0.17 SD.
- Conclusion: coefficient posteriors are data-dominated apart from weakly identified college terms on rare transitions; sd 5 is retained as the main prior. Life-table check below.

**Life-table sensitivity to the prior (2026-10-06).** `R/05_life_tables.R` was run for all three priors and both sexes (population and status-based radix; ages 50, 70, 90; expectancies by state, total, and % of remaining life). Relative to sd 5, the largest absolute difference in any posterior mean is 0.16 years (sd 2.5) and 0.07 years (sd 10), i.e. at most 0.14 and 0.05 of the sd 5 95% CrI half-width. Total life expectancy differs by <=0.005 years. The largest relative shift is cognitive-only expectancy at age 50 for men with a status-based radix in state 4 (0.73 vs 0.77 under sd 2.5, inside the sd 5 95% CrI of 0.51-1.09). The prior-driven shifts in the rare-transition college terms do not materially affect life-table results; sd 5 retained.
