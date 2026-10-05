# Build outcome codes and design matrix for the transition model (decision log: 2026-10-04).
# Outcome code t = (start - 1) * 5 + end  (1..20; 1 = functional -> functional, the reference).
# Predictors: intercept, age (centred 70, per decade), age^2, race/ethnicity (3), education (3),
# birth cohort (per decade, centred at the sample mean). Sex is the stratifier, not a predictor.
age_centre <- 70; age_scale <- 10

make_design <- function(iv, cohort_centre) {
  y <- (iv$start_state - 1L) * 5L + iv$end_state
  a <- (iv$age - age_centre) / age_scale
  X <- cbind(int = 1, age = a, age2 = a^2,
             black = as.numeric(iv$race_eth == "NH Black"),
             hisp  = as.numeric(iv$race_eth == "Hispanic"),
             othr  = as.numeric(iv$race_eth == "NH Other"),
             hs    = as.numeric(iv$educ == "HS/GED"),
             somec = as.numeric(iv$educ == "Some college"),
             coll  = as.numeric(iv$educ == "College+"),
             coh   = (iv$birth_year - cohort_centre) / 10)
  list(y = y, X = X)
}
