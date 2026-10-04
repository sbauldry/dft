# 01_import.R -- read selected RAND HRS and Langa-Weir variables, link, save wide file
# Output: data/derived/hrs_wide.rds
suppressPackageStartupMessages({library(tidyverse); library(haven); library(here)})
source(here("R/functions/constants.R"))

rand_path <- here("data/raw/randhrs1992_2022v1.dta")
lw_path   <- here("data/raw/cogfinalimp_9522wide.dta")

# ---- RAND ----
rand_vars <- c("hhidpn", "hhid", "hacohort", "rabyear", "rabmonth", "rabdate", "raddate",
               "radyear", "radmonth", "ragender", "raracem", "rahispan", "raeduc", "rabplace",
               paste0("r", waves, "iwstat"), paste0("r", waves, "iwend"),
               paste0("r", waves, "agem_e"), paste0("r", waves, "adl5a"))
rand <- read_dta(rand_path, col_select = any_of(rand_vars)) |>
  mutate(across(everything(), ~ as.numeric(zap_labels(.x))),
         hhid = NULL)  # RAND hhid is not needed; linkage is via hhidpn
stopifnot(all(rand_vars %in% c(names(rand), "hhid")), !anyDuplicated(rand$hhidpn))

# ---- Langa-Weir ----
yrs <- wave_year(waves)
lw <- read_dta(lw_path, col_select = c(hhid, pn, paste0("cogfunction", yrs), paste0("proxy", yrs))) |>
  mutate(hhidpn = as.numeric(paste0(hhid, pn))) |>
  select(-hhid, -pn) |>
  mutate(across(-hhidpn, ~ as.numeric(zap_labels(.x))))
stopifnot(!anyDuplicated(lw$hhidpn))

# ---- Link checks ----
cat("RAND rows:", nrow(rand), " LW rows:", nrow(lw), "\n")
cat("LW ids found in RAND:", sum(lw$hhidpn %in% rand$hhidpn), "of", nrow(lw), "\n")
cat("RAND ids found in LW:", sum(rand$hhidpn %in% lw$hhidpn), "of", nrow(rand), "\n")

wide <- left_join(rand, lw, by = "hhidpn")

# Cognition coverage vs. interview status (RAND iwstat == 1 expected to have a LW class)
cov <- map_dfr(waves, function(w) {
  r <- wide[[paste0("r", w, "iwstat")]] == 1 & !is.na(wide[[paste0("r", w, "iwstat")]])
  c <- wide[[paste0("cogfunction", wave_year(w))]]
  tibble(wave = w, year = wave_year(w), resp = sum(r), resp_with_cog = sum(r & !is.na(c)),
         nonresp_with_cog = sum(!r & !is.na(c)))
})
print(cov, n = Inf)

dir.create(here("data/derived"), showWarnings = FALSE)
saveRDS(wide, here("data/derived/hrs_wide.rds"))
