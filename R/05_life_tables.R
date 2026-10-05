# 05_life_tables.R -- multistate life tables from posterior draws (Lynch & Zang 2022), one per sex.
# Usage: Rscript R/05_life_tables.R <men|women> [prior_sd=5] [n_draws=all] [cores=8]
# Profile: ages 50-110+ in 2-year steps; birth cohort at the sample mean (centred coefficient = 0);
#   race/ethnicity and education dummies at sample proportions (unweighted, this sex's intervals).
# Radix: population-based and status-based (start in each of states 1-4 at age 50).
# Quantities per draw: state expectancies e_s(a), total life expectancy, % of remaining life in each state.
# Output: output/models/lifetable_<sex>_psd<sd>.rds  (all-age summaries + draw-level values at ages 50,60,...,90)
#         output/tables/lifetable_summary_<sex>_psd<sd>.csv  (posterior mean, median, 95% and 84% CrIs)
# Never overwrites existing output.
suppressPackageStartupMessages({library(tidyverse); library(here); library(parallel)})
source(here("R/functions/constants.R")); source(here("R/functions/design.R"))
source(here("R/functions/life_table.R"))

args     <- commandArgs(trailingOnly = TRUE)
sex      <- match.arg(if (length(args) >= 1) args[1] else "men", c("men", "women"))
prior_sd <- if (length(args) >= 2) as.numeric(args[2]) else 5
n_use    <- if (length(args) >= 3 && args[3] != "all") as.integer(args[3]) else NA_integer_
cores    <- if (length(args) >= 4) as.integer(args[4]) else 8L

tag <- sprintf("%s_psd%s", sex, prior_sd)
out_rds <- here("output/models", sprintf("lifetable_%s.rds", tag))
out_csv <- here("output/tables", sprintf("lifetable_summary_%s.csv", tag))
if (is.na(n_use) && (file.exists(out_rds) || file.exists(out_csv))) stop("Refusing to overwrite existing output for ", tag)

fit <- readRDS(here("output/models", sprintf("draws_%s_psd%s.rds", sex, prior_sd)))
dm <- dim(fit$draws)                                             # iter x chain x P x K
B_all <- aperm(fit$draws, c(3, 4, 1, 2)); dim(B_all) <- c(dm[3], dm[4], dm[1] * dm[2])   # P x K x G
G <- dim(B_all)[3]
if (!is.na(n_use)) B_all <- B_all[, , round(seq(1, G, length.out = n_use)), drop = FALSE]
G <- dim(B_all)[3]

# ---- covariate profile ----
iv <- readRDS(here("data/derived/intervals.rds")) |> filter(female == as.integer(sex == "women"))
stopifnot(isTRUE(all.equal(mean(iv$birth_year), fit$cohort_centre)), nrow(iv) == fit$N)
Xs <- make_design(iv, fit$cohort_centre)$X
prof <- colMeans(Xs[, c("black", "hisp", "othr", "hs", "somec", "coll")])
ages <- seq(50, 110, by = 2); k <- 2; S <- 4
a <- (ages - fit$age_centre) / fit$age_scale
Z <- cbind(int = 1, age = a, age2 = a^2, matrix(prof, length(ages), 6, byrow = TRUE, dimnames = list(NULL, names(prof))),
           coh = 0)
stopifnot(identical(colnames(Z), fit$pnames))
cat(sprintf("%s: %d draws, %d ages; profile (black, hisp, othr, hs, somec, coll) = %s\n", sex, G, length(ages),
            paste(round(prof, 3), collapse = ", ")))

radices <- list(population = "population", s1 = 1L, s2 = 2L, s3 = 3L, s4 = 4L)
one_draw <- function(g) {
  B <- B_all[, , g]
  sapply(radices, function(r) as.vector(mslt_draw(Z, B, S, fit$ref, fit$cats, k = k, radix = r)$e),
         simplify = "array")                                       # (n_age*S) x radix
}
t_run <- system.time(res <- mclapply(seq_len(G), one_draw, mc.cores = cores))
if (any(vapply(res, inherits, TRUE, "try-error"))) stop("a draw failed")
cat("life tables computed in", round(t_run["elapsed"]), "s\n")

E <- array(unlist(res), c(length(ages), S, length(radices), G),
           dimnames = list(age = ages, state = state_labels[1:S], radix = names(radices), draw = NULL))
stopifnot(all(is.finite(E)), all(E >= 0))
TLE <- apply(E, c(1, 3, 4), sum)                                    # total life expectancy
PCT <- sweep(E, c(1, 3, 4), TLE, "/") * 100                          # % of remaining life by state

qs <- function(x) c(mean = mean(x), median = median(x), lo95 = quantile(x, .025, names = FALSE),
                    hi95 = quantile(x, .975, names = FALSE), lo84 = quantile(x, .08, names = FALSE),
                    hi84 = quantile(x, .92, names = FALSE))
summ_arr <- function(arr, quantity) {                                # arr: age x [state] x radix x draw
  d <- dim(arr); if (length(d) == 3) { arr <- array(arr, c(d[1], 1, d[2], d[3]),
                                       dimnames = list(dimnames(arr)[[1]], "Total", dimnames(arr)[[2]], NULL)) }
  apply(arr, 1:3, qs) |> as.data.frame.table(responseName = "value", stringsAsFactors = FALSE) |>
    setNames(c("stat", "age", "state", "radix", "value")) |>
    pivot_wider(names_from = stat, values_from = value) |> mutate(quantity = quantity, age = as.integer(age), .before = 1)
}
tab <- bind_rows(summ_arr(E, "expectancy"), summ_arr(TLE, "total_expectancy"), summ_arr(PCT, "pct_remaining")) |>
  mutate(sex = sex, prior_sd = prior_sd, .before = 1)
write_csv(tab, out_csv)

sel <- as.character(seq(50, 90, by = 10))
saveRDS(list(summary = tab, E_sel = E[sel, , , , drop = FALSE], ages = ages, profile = prof, radices = names(radices),
             sex = sex, prior_sd = prior_sd, n_draws = G,
             date = Sys.time()), out_rds)
cat("saved", out_csv, "and", out_rds, "\n")
print(tab |> filter(quantity == "total_expectancy", radix == "population", age %in% c(50, 70, 90)) |>
        select(age, mean, lo95, hi95, lo84, hi84))
