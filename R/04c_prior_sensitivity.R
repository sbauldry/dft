# 04c_prior_sensitivity.R -- compare posterior means/SDs across priors N(0, sd^2), sd in {2.5, 5, 10}.
# Output: output/tables/prior_sensitivity_<sex>.csv, printed summary
suppressPackageStartupMessages({library(tidyverse); library(here)})
sds <- c("2.5", "5", "10")
summ <- function(sex, psd) {
  f <- readRDS(here("output/models", sprintf("draws_%s_psd%s.rds", sex, psd)))
  d <- f$draws; dm <- dim(d)
  m <- matrix(d, dm[1] * dm[2], dm[3] * dm[4])
  tibble(par = as.vector(outer(f$pnames, colnames(d[1, 1, , , drop = TRUE]) %||% paste0("t", f$cats), paste, sep = ":")),
         mean = colMeans(m), sd = apply(m, 2, sd))
}
for (sex in c("men", "women")) {
  tab <- map(sds, \(p) summ(sex, p) |> rename_with(~ paste0(.x, "_", p), -par)) |> reduce(left_join, by = "par") |>
    mutate(d25 = (mean_2.5 - mean_5) / sd_5, d10 = (mean_10 - mean_5) / sd_5)
  write_csv(tab, here("output/tables", sprintf("prior_sensitivity_%s.csv", sex)))
  cat("\n==", sex, "(differences vs sd 5, in sd-5 posterior SDs)\n")
  print(tab |> summarise(across(c(d25, d10), list(med_abs = ~ median(abs(.x)), max_abs = ~ max(abs(.x))))))
  cat("Largest 5 by |d25| or |d10|:\n")
  print(tab |> mutate(mx = pmax(abs(d25), abs(d10))) |> slice_max(mx, n = 5) |>
          select(par, mean_2.5, mean_5, mean_10, sd_5, d25, d10), width = 200)
  cat("Posterior SD ratio vs sd 5 (median):", round(median(tab$sd_2.5 / tab$sd_5), 3), "(sd 2.5);",
      round(median(tab$sd_10 / tab$sd_5), 3), "(sd 10)\n")
}
