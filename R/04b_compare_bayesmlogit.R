# 04b_compare_bayesmlogit.R -- compare our sampler and life-table code with the bayesmlogit package
# (Lynch, Zang & Zhang; CRAN 1.0.1) on the same data.
# Data: men, stratified subsample of 100 intervals per transition cell (N = 2,000), covariates:
#   age (years since 50) and an NH Black indicator; outcome = transition code (start-1)*5 + end; ref = 20 (4->5),
#   the last category, which is what bayesmlogit uses.
# Output: output/models/validation_bayesmlogit.rds
suppressPackageStartupMessages({library(tidyverse); library(here); library(bayesmlogit); library(parallel); library(nnet); library(posterior)})
source(here("R/functions/life_table.R"))

set.seed(20261004)
iv <- readRDS(here("data/derived/intervals.rds")) |> filter(female == 0)
sub <- iv |> mutate(y = (start_state - 1L) * 5L + end_state,
                    age = age - 50, black = as.numeric(race_eth == "NH Black")) |>
  group_by(y) |> slice_sample(n = 100) |> ungroup()
y <- sub$y; X <- as.matrix(sub[, c("age", "black")])
cat("N =", nrow(sub), " categories:", n_distinct(y), "\n")

# ---- bayesmlogit (single chain, flat prior, serial R loop) ----
t_b <- system.time(fit_b <- bayesmlogit(y, X, samp = 1000, burn = 500, verbose = 100, thin = 5))
cat("bayesmlogit time (s):", round(t_b["elapsed"]), "\n")
draws_b <- as.matrix(fit_b$out)                               # 1000 x (P*(J-1)), P blocks per outcome

# ---- ours: 4 chains, same data, same reference ----
Xi <- cbind(1, X)
chains <- mclapply(1:4, function(ch) {
  t <- system.time(f <- pg_mlogit(y, Xi, J = 20, ref = 20, n_iter = 3000, burn = 1000, thin = 2,
                                  prior_sd = 5, seed = 100 + ch))
  list(fit = f, time = t["elapsed"])
}, mc.cores = 4)
cat("our time per chain (s):", round(chains[[1]]$time), "\n")
cats <- chains[[1]]$fit$cats
# our array is [iter, P, J-1]; bayesmlogit's columns are P-blocks per outcome (P fastest)
draws_o <- do.call(rbind, lapply(chains, function(c) {
  a <- c$fit$draws; do.call(cbind, lapply(seq_len(dim(a)[3]), function(j) a[, , j])) }))

cmp <- tibble(par = seq_len(ncol(draws_b)),
              outcome = cats[ceiling(par / 3)], term = c("int", "age", "black")[(par - 1) %% 3 + 1],
              mean_bm = colMeans(draws_b), sd_bm = apply(draws_b, 2, sd),
              mean_ours = colMeans(draws_o), sd_ours = apply(draws_o, 2, sd)) |>
  mutate(diff_in_sd = (mean_ours - mean_bm) / pmax(sd_bm, sd_ours))
# ---- ML estimates (reference = category 20) and MCMC diagnostics ----
ml <- multinom(relevel(factor(y), ref = "20") ~ age + black, data = data.frame(y = y, X), trace = FALSE, maxit = 500)
mlv <- as.vector(t(coef(ml)))                                  # outcome blocks of (int, age, black), cats order
cmp$mle <- mlv
cmp$bm_vs_mle_sd  <- (cmp$mean_bm   - cmp$mle) / cmp$sd_bm
cmp$our_vs_mle_sd <- (cmp$mean_ours - cmp$mle) / cmp$sd_ours
cat("\nPosterior mean minus MLE, in posterior SDs: bayesmlogit max", round(max(abs(cmp$bm_vs_mle_sd)), 2),
    " median", round(median(abs(cmp$bm_vs_mle_sd)), 2), "| ours max", round(max(abs(cmp$our_vs_mle_sd)), 2),
    " median", round(median(abs(cmp$our_vs_mle_sd)), 2), "\n")
cmp$ess_bm <- apply(draws_b, 2, function(x) ess_bulk(x))
ours_arr <- array(NA_real_, c(nrow(chains[[1]]$fit$draws), 4, ncol(draws_o)))
for (ch in 1:4) { a <- chains[[ch]]$fit$draws; ours_arr[, ch, ] <- do.call(cbind, lapply(seq_len(dim(a)[3]), function(j) a[, , j])) }
cmp$rhat_ours <- apply(ours_arr, 3, function(x) rhat(x)); cmp$ess_ours <- apply(ours_arr, 3, function(x) ess_bulk(x))
cat("bayesmlogit bulk ESS (of 1000 draws): median", round(median(cmp$ess_bm)), " min", round(min(cmp$ess_bm)), "\n")
cat("ours (4 chains, 4 x", dim(ours_arr)[1], "draws): R-hat max", round(max(cmp$rhat_ours), 3),
    " bulk ESS median", round(median(cmp$ess_ours)), " min", round(min(cmp$ess_ours)), "\n")
saveRDS(list(cmp = cmp, draws_b = draws_b, draws_o = draws_o, time_bm = t_b), here("output/models/validation_bayesmlogit.rds"))
cat("\nCoefficients compared:", nrow(cmp), "\n")
cat("max |mean diff| / sd:", round(max(abs(cmp$diff_in_sd)), 2), " median:", round(median(abs(cmp$diff_in_sd)), 2), "\n")
cat("sd ratio (ours/bm): median", round(median(cmp$sd_ours / cmp$sd_bm), 2),
    " range", paste(round(range(cmp$sd_ours / cmp$sd_bm), 2), collapse = "-"), "\n")
cat("max |mean diff| (raw):", round(max(abs(cmp$mean_ours - cmp$mean_bm)), 3), "\n")

# ---- Life tables from the SAME draws: package mlifeTable vs ours (age 50-110, k = 2) ----
nd <- 50
tr <- draws_b[round(seq(1, nrow(draws_b), length.out = nd)), ]
td <- file.path(tempdir(), "mlt"); dir.create(td)
mlifeTable(y, as.data.frame(X), trans = tr, states = 5, file_path = td,
           startages = 0, endages = 60, age.gap = 2, nums = nd)
e_pkg <- as.matrix(read.table(list.files(td, "lifetable", full.names = TRUE)[1], header = TRUE))[, 1:4]
ages <- seq(0, 60, by = 2); xbar <- c(1, NA, mean(X[, "black"]))
e_our <- t(sapply(seq_len(nd), function(i) {
  B <- matrix(tr[i, ], 3); Z <- cbind(1, ages, xbar[3])
  mslt_draw(Z, B, S = 4, ref = 20, cats = cats, k = 2)$e[1, ]
}))
cat("\nLife table, e(50) by state, 50 shared draws: max abs difference (years):",
    round(max(abs(e_pkg - e_our)), 5), "\n")
print(round(rbind(pkg = colMeans(e_pkg), ours = colMeans(e_our)), 3))

# posterior of e(50) under each sampler
e_post <- function(dr, n = 200) { idx <- round(seq(1, nrow(dr), length.out = n))
  t(sapply(idx, function(i) mslt_draw(cbind(1, ages, xbar[3]), matrix(dr[i, ], 3), 4, 20, cats)$e[1, ])) }
eb <- e_post(draws_b); eo <- e_post(draws_o)
cat("\nPosterior mean (sd) of e(50), bayesmlogit vs ours:\n")
print(round(rbind(bm_mean = colMeans(eb), ours_mean = colMeans(eo), bm_sd = apply(eb, 2, sd), ours_sd = apply(eo, 2, sd)), 3))

dir.create(here("output/models"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(cmp = cmp, e_pkg = e_pkg, e_our = e_our, time_bm = t_b, draws_b = draws_b, draws_o = draws_o),
        here("output/models/validation_bayesmlogit.rds"))
