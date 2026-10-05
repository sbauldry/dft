# Polya-Gamma Gibbs sampler for the Bayesian multinomial logit (Polson, Scott & Windle 2013),
# as used by Lynch & Zang (2022). Category `ref` has coefficients fixed at 0.
# Each non-reference category j is updated given the others via the conditional binary logit
#   y_ij | beta_j, beta_-j ~ Bernoulli(logistic(x_i beta_j - c_ij)),  c_ij = log sum_{k != j} exp(x_i beta_k)
# with omega_ij ~ PG(1, x_i beta_j - c_ij) and a conjugate normal update for beta_j.
suppressPackageStartupMessages(library(BayesLogit))

# y: integer vector in 1..J; X: N x P matrix (include intercept column).
# prior_sd: scalar or length-P vector, independent N(0, prior_sd^2) on every coefficient.
# Returns list(draws = array[n_keep, P, J-1], cats = non-reference categories, ...)
pg_mlogit <- function(y, X, J = max(y), ref = 1L, n_iter = 3000, burn = 1000, thin = 2,
                      prior_sd = 5, init_sd = 0.5, seed = NULL, verbose = 0) {
  if (!is.null(seed)) set.seed(seed)
  X <- as.matrix(X); N <- nrow(X); P <- ncol(X)
  stopifnot(length(y) == N, all(y %in% seq_len(J)), !anyNA(X))
  cats <- setdiff(seq_len(J), ref)
  P0 <- diag(1 / rep_len(prior_sd, P)^2, P)
  B <- matrix(0, P, J)
  B[, cats] <- rnorm(P * length(cats), 0, init_sd)            # dispersed starting values
  eta <- X %*% B                                              # N x J linear predictors (ref column = 0)
  keep_iter <- seq(burn + thin, n_iter, by = thin)
  draws <- array(NA_real_, c(length(keep_iter), P, length(cats)))
  k <- 0L
  for (it in seq_len(n_iter)) {
    for (jj in seq_along(cats)) {
      j <- cats[jj]
      eo <- eta[, -j, drop = FALSE]
      m  <- do.call(pmax, as.data.frame(eo))
      cj <- m + log(rowSums(exp(eo - m)))                     # log-sum-exp over other categories
      psi <- eta[, j] - cj
      w  <- rpg(N, 1, psi)
      kappa <- (y == j) - 0.5
      Xw <- X * sqrt(w)
      P1 <- crossprod(Xw) + P0
      R  <- chol(P1)
      mu <- backsolve(R, forwardsolve(t(R), crossprod(X, kappa + w * cj)))
      bj <- as.vector(mu + backsolve(R, rnorm(P)))
      B[, j] <- bj
      eta[, j] <- X %*% bj
    }
    if (it %in% keep_iter) { k <- k + 1L; draws[k, , ] <- B[, cats] }
    if (verbose > 0 && it %% verbose == 0) cat("iter", it, "of", n_iter, "\n")
  }
  list(draws = draws, cats = cats, ref = ref, J = J, P = P, colnames = colnames(X))
}

# Transition probability matrices from one coefficient draw.
# Z: n_age x P design rows; B: P x (J-1) draw; returns n_age x J matrix of joint probabilities
# (rows sum to 1; column `ref` is the reference outcome).
pg_probs <- function(Z, B, J, ref, cats) {
  eta <- matrix(0, nrow(Z), J); eta[, cats] <- Z %*% B
  m <- apply(eta, 1, max); e <- exp(eta - m)
  e / rowSums(e)
}
