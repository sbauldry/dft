# 04_estimate.R -- Bayesian multinomial logit for the 20 start x end transitions, independence
# Metropolis-Hastings sampler (decision log 2026-10-05). One model per sex.
# Usage: Rscript R/04_estimate.R <men|women> [prior_sd=5] [n_iter=100000] [burn=10000] [thin=10] [df=30]
# Output: output/models/draws_<sex>_psd<prior_sd>.rds  (never overwritten; script stops if it exists)
suppressPackageStartupMessages({library(tidyverse); library(here); library(parallel); library(posterior)})
source(here("R/functions/constants.R")); source(here("R/functions/design.R"))
source(here("R/functions/mh_mlogit.R"))

args     <- commandArgs(trailingOnly = TRUE)
sex      <- match.arg(if (length(args) >= 1) args[1] else "men", c("men", "women"))
prior_sd <- if (length(args) >= 2) as.numeric(args[2]) else 5
n_iter   <- if (length(args) >= 3) as.integer(args[3]) else 100000L
burn     <- if (length(args) >= 4) as.integer(args[4]) else 10000L
thin     <- if (length(args) >= 5) as.integer(args[5]) else 10L
df       <- if (length(args) >= 6) as.numeric(args[6]) else 30
n_chain  <- 4L
base_seed <- 20261005L + (sex == "women") * 1000L

out_file <- here("output/models", sprintf("draws_%s_psd%s.rds", sex, prior_sd))
if (file.exists(out_file)) stop("Refusing to overwrite existing run: ", out_file)

iv <- readRDS(here("data/derived/intervals.rds")) |> filter(female == as.integer(sex == "women"))
cohort_centre <- mean(iv$birth_year)
des <- make_design(iv, cohort_centre)
y <- des$y; X <- des$X
cat(sprintf("%s: N = %d intervals, %d persons; prior sd %s; %d chains x %d proposals (burn %d, thin %d, t df %s)\n",
            sex, nrow(iv), n_distinct(iv$hhidpn), prior_sd, n_chain, n_iter, burn, thin, df))
stopifnot(all(table(y) > 0), length(unique(y)) == 20)

# posterior mode / proposal covariance computed once, shared by all chains
t_mode <- system.time(mode <- mlogit_mode(X, y, setdiff(1:20, 1L), 1L, prior_sd))
cat("mode found in", round(t_mode["elapsed"]), "s\n")

t_run <- system.time(chains <- mclapply(seq_len(n_chain), function(ch)
  mh_mlogit(y, X, J = 20, ref = 1L, n_iter = n_iter, burn = burn, thin = thin, prior_sd = prior_sd,
            df = df, seed = base_seed + ch, mode = mode, verbose = 10000),
  mc.cores = n_chain, mc.set.seed = FALSE))
if (any(vapply(chains, inherits, TRUE, "try-error"))) stop("a chain failed")
cat("sampling time (s):", round(t_run["elapsed"]), "\n")

draws <- abind::abind(lapply(chains, `[[`, "draws"), along = 0)        # chain x iter x P x K
draws <- aperm(draws, c(2, 1, 3, 4))                                      # iter x chain x P x K
cats <- chains[[1]]$cats
pnames <- colnames(X)
dimnames(draws) <- list(NULL, NULL, pnames, paste0("t", cats))
accept <- vapply(chains, `[[`, 0, "accept")
cat("acceptance by chain:", round(accept, 3), "\n")

# convergence summary
dv <- draws; dim(dv) <- c(dim(draws)[1:2], prod(dim(draws)[3:4]))
dimnames(dv)[[3]] <- as.vector(outer(pnames, paste0("t", cats), paste, sep = ":"))
sm <- summarise_draws(as_draws_array(dv), "rhat", "ess_bulk", "ess_tail")
cat(sprintf("max R-hat %.3f; bulk ESS median %.0f, min %.0f; tail ESS min %.0f\n",
            max(sm$rhat), median(sm$ess_bulk), min(sm$ess_bulk), min(sm$ess_tail)))

saveRDS(list(draws = draws, cats = cats, ref = 1L, pnames = pnames, sex = sex, prior_sd = prior_sd,
             n_iter = n_iter, burn = burn, thin = thin, df = df, seeds = base_seed + seq_len(n_chain),
             accept = accept, cohort_centre = cohort_centre, age_centre = age_centre, age_scale = age_scale,
             N = nrow(iv), n_persons = n_distinct(iv$hhidpn), mode = mode$B, diagnostics = sm,
             run_time = t_run["elapsed"], date = Sys.time()), out_file)
cat("saved", out_file, "\n")
