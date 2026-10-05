# Independence Metropolis-Hastings sampler for the Bayesian multinomial logit (Lynch & Zang 2022,
# "independence sampler"): multivariate-t proposal centred at the posterior mode with covariance =
# inverse negative Hessian of the log posterior. Same model/prior as pg_mlogit(): ref category has
# coefficients 0; independent N(0, prior_sd^2) priors on all others.

mlogit_pieces <- function(B, X, y, cats, ref, prior_sd, need_hess = FALSE) {
  N <- nrow(X); P <- ncol(X); K <- length(cats)
  eta <- X %*% B                                      # N x K
  m <- pmax(0, do.call(pmax, as.data.frame(eta)))
  lse <- m + log(exp(-m) + rowSums(exp(eta - m)))
  pr <- exp(eta - lse)
  Y <- outer(y, cats, "==") * 1
  ll <- sum(eta * Y) - sum(lse)
  lp <- ll - 0.5 * sum(B^2) / prior_sd^2
  out <- list(lp = lp)
  if (need_hess) {
    G <- crossprod(X, Y - pr) - B / prior_sd^2        # P x K gradient
    H <- matrix(0, P * K, P * K)
    for (j in seq_len(K)) for (k in j:K) {
      w <- if (j == k) pr[, j] * (1 - pr[, j]) else -pr[, j] * pr[, k]
      blk <- -crossprod(X * w, X)                      # d2 ll / dB_j dB_k'
      H[(j - 1) * P + 1:P, (k - 1) * P + 1:P] <- blk
      if (k != j) H[(k - 1) * P + 1:P, (j - 1) * P + 1:P] <- t(blk)
    }
    diag(H) <- diag(H) - 1 / prior_sd^2
    out$G <- G; out$H <- H
  }
  out
}

mlogit_mode <- function(X, y, cats, ref, prior_sd = 5, tol = 1e-6, maxit = 100, verbose = FALSE) {
  P <- ncol(X); K <- length(cats); B <- matrix(0, P, K)
  pc <- mlogit_pieces(B, X, y, cats, ref, prior_sd, need_hess = TRUE)
  for (it in seq_len(maxit)) {
    step <- matrix(solve(-pc$H, as.vector(pc$G)), P, K)
    s <- 1                                                   # backtracking line search (step halving)
    repeat {
      cand <- mlogit_pieces(B + s * step, X, y, cats, ref, prior_sd)
      if (is.finite(cand$lp) && cand$lp >= pc$lp - 1e-8) break
      s <- s / 2; if (s < 1e-10) break
    }
    B <- B + s * step
    pc <- mlogit_pieces(B, X, y, cats, ref, prior_sd, need_hess = TRUE)
    if (verbose) cat("Newton it", it, " lp =", round(pc$lp, 3), " step scale", signif(s, 2), " max|step| =", signif(max(abs(s * step)), 3), "\n")
    if (max(abs(s * step)) < tol) break
  }
  list(B = B, H = pc$H, lp = pc$lp)
}

# n_iter proposals; returns draws array [n_keep, P, K] after burn/thin, acceptance rate
mh_mlogit <- function(y, X, J = max(y), ref = 1L, n_iter = 10000, burn = 500, thin = 1, prior_sd = 5,
                      df = 8, seed = NULL, mode = NULL, verbose = 0) {
  if (!is.null(seed)) set.seed(seed)
  cats <- setdiff(seq_len(J), ref); P <- ncol(X); K <- length(cats); d <- P * K
  if (is.null(mode)) mode <- mlogit_mode(X, y, cats, ref, prior_sd)
  Rch <- chol(-mode$H)                                 # -H = R'R ; covariance = (R'R)^-1
  mu <- as.vector(mode$B)
  rprop <- function() mu + backsolve(Rch, rnorm(d)) / sqrt(rchisq(1, df) / df)
  logq <- function(th) { q <- sum((Rch %*% (th - mu))^2); -(df + d) / 2 * log1p(q / df) }
  cur <- mu; lp_cur <- mlogit_pieces(matrix(cur, P, K), X, y, cats, ref, prior_sd)$lp; lq_cur <- logq(cur)
  keep_iter <- seq(burn + thin, n_iter, by = thin)
  draws <- array(NA_real_, c(length(keep_iter), P, K)); k <- 0L; acc <- 0L
  for (it in seq_len(n_iter)) {
    cand <- rprop()
    lp_c <- mlogit_pieces(matrix(cand, P, K), X, y, cats, ref, prior_sd)$lp; lq_c <- logq(cand)
    if (log(runif(1)) < (lp_c - lp_cur) - (lq_c - lq_cur)) { cur <- cand; lp_cur <- lp_c; lq_cur <- lq_c; acc <- acc + 1L }
    if (it %in% keep_iter) { k <- k + 1L; draws[k, , ] <- matrix(cur, P, K) }
    if (verbose > 0 && it %% verbose == 0) cat("iter", it, " accept", round(acc / it, 3), "\n")
  }
  list(draws = draws, cats = cats, ref = ref, J = J, P = P, accept = acc / n_iter, mode = mode)
}
