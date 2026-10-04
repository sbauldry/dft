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
5. **Radix for population-based tables.** "Row sums of the unnormalized probabilities at age 0 (a=1)" is the only description. Whether this includes the death column (probability of dying within the first interval) is not explicit; verify against the supplement code.
6. **Draws and thinning.** 1,000 final draws from 2 chains; we may choose differently. No guidance on how many are needed for stable interval endpoints.
7. **Treatment of death within an interval.** Linear approximation L = 0.5k[l(a)+l(a+1)] applies to all states including those who die mid-interval; no separate death-timing adjustment (e.g., our `radyear`/`radmonth`) is described.
8. **Supplement and R functions.** Article states R functions and instructions are available from the authors; obtain them to check against these notes.
