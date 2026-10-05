# Multistate life table from one posterior draw (Lynch & Zang 2022, Step 2-3).
# States: 1..S living, S+1 dead. Transition outcome code t = (start - 1) * (S + 1) + end,
# i.e. J = S * (S + 1) outcomes (all start x end combinations for living starts).
# Z: n_age x P design matrix (rows = ages on a grid with spacing k, covariates fixed at the profile)
# B: P x (J-1) coefficient draw for the non-reference outcomes `cats`.
# radix: "population" (row sums of unnormalised probabilities at the first age; includes the
#        death column, as in bayesmlogit::mlifeTable) or an integer living state (status-based).
# Returns list(e = n_age x S matrix of state expectancies at each age, P = array of normalised matrices).
source(here::here("R/functions/pg_mlogit.R"))

mslt_draw <- function(Z, B, S, ref, cats, k = 2, radix = "population") {
  J <- S * (S + 1); D <- S + 1; n <- nrow(Z)
  pr <- pg_probs(Z, B, J, ref, cats)                       # n x J joint probabilities
  Pm <- array(0, c(D, D, n))
  for (s in seq_len(S)) Pm[s, , ] <- t(pr[, (s - 1) * D + seq_len(D)])
  Pm[D, D, ] <- 1
  l1 <- if (identical(radix, "population")) { r <- rowSums(Pm[, , 1]); r[D] <- 0; r }
        else { r <- numeric(D); r[radix] <- 1; r }
  for (a in seq_len(n)) Pm[seq_len(S), , a] <- Pm[seq_len(S), , a] / rowSums(Pm[seq_len(S), , a])
  l <- matrix(0, n, D); l[1, ] <- l1
  for (a in seq_len(n - 1)) l[a + 1, ] <- l[a, ] %*% Pm[, , a]
  L <- matrix(0, n, D)
  for (a in seq_len(n - 1)) L[a, ] <- 0.5 * k * (l[a, ] + l[a + 1, ])
  Pn <- Pm[seq_len(S), seq_len(S), n]
  L[n, seq_len(S)] <- k * (l[n, seq_len(S)] %*% solve(diag(S) - Pn))
  Tt <- apply(L[, seq_len(S), drop = FALSE], 2, function(x) rev(cumsum(rev(x))))
  e <- Tt / rowSums(l[, seq_len(S), drop = FALSE])
  list(e = e, l = l)
}
