# Simulation test: recover known coefficients; compare with nnet::multinom ML estimates.
suppressPackageStartupMessages({library(here); library(nnet)})
source(here("R/functions/pg_mlogit.R"))
set.seed(20261004)
N <- 20000; P <- 3; J <- 6
X <- cbind(1, rnorm(N), rbinom(N, 1, .4))
Btrue <- matrix(0, P, J); Btrue[, 2:J] <- matrix(rnorm(P * (J - 1), 0, 0.8), P)
pr <- exp(X %*% Btrue); pr <- pr / rowSums(pr)
y <- apply(pr, 1, function(p) sample(J, 1, prob = p))
t0 <- Sys.time()
fit <- pg_mlogit(y, X, J = J, ref = 1, n_iter = 2000, burn = 500, thin = 1, prior_sd = 5, seed = 1)
cat("time:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")
pm <- apply(fit$draws, c(2, 3), mean); psd <- apply(fit$draws, c(2, 3), sd)
ml <- multinom(factor(y) ~ X[, -1], trace = FALSE); mlb <- t(coef(ml))
cat("max |post mean - truth|:", round(max(abs(pm - Btrue[, -1])), 3), "\n")
cat("max |post mean - MLE|  :", round(max(abs(pm - mlb)), 3), "\n")
cat("truth within 95% interval:", mean(sapply(seq_len(ncol(pm)), function(j) sapply(seq_len(P), function(p) {
  q <- quantile(fit$draws[, p, j], c(.025, .975)); Btrue[p, j + 1] >= q[1] & Btrue[p, j + 1] <= q[2] }))), "\n")
