# 02_construct_states.R -- person-wave file with ADL, cognition and 4-state classification
# Input:  data/derived/hrs_wide.rds
# Output: data/derived/person_wave.rds ; output/tables/02_state_by_wave.csv
# States (living): 1 functional, 2 physical only (ADL>=1), 3 cognitive only (dementia), 4 both.
# Main cognitive cut: Langa-Weir dementia only (cogfunction == 3).
# Sensitivity cut:    CIND + dementia   (cogfunction >= 2) -> state_cind.
suppressPackageStartupMessages({library(tidyverse); library(here)})
source(here("R/functions/constants.R"))

wide <- readRDS(here("data/derived/hrs_wide.rds"))

person <- wide |>
  select(hhidpn, hacohort, rabyear, rabdate, raddate, radyear, radmonth,
         ragender, raracem, rahispan, raeduc, rabplace)

pw <- map_dfr(waves, function(w) {
  tibble(hhidpn = wide$hhidpn, wave = w, year = wave_year(w),
         iwstat = wide[[paste0("r", w, "iwstat")]],
         iwend  = wide[[paste0("r", w, "iwend")]],
         age    = wide[[paste0("r", w, "agem_e")]] / 12,   # decimal age at interview end
         adl5   = wide[[paste0("r", w, "adl5a")]],
         cog    = wide[[paste0("cogfunction", wave_year(w))]],
         proxy  = wide[[paste0("proxy", wave_year(w))]])
}) |>
  filter(!is.na(iwstat), iwstat != 0)

pw <- pw |>
  mutate(adl  = if_else(adl5 > 0, 1L, 0L),
         dem  = if_else(cog == 3, 1L, 0L),
         cind = if_else(cog >= 2, 1L, 0L),
         state      = if_else(iwstat == 1, 1L + adl + 2L * dem,  NA_integer_),
         state_cind = if_else(iwstat == 1, 1L + adl + 2L * cind, NA_integer_))

stopifnot(all(pw$state %in% c(1:4, NA)), all(pw$adl5 %in% c(0:5, NA)), all(pw$cog %in% c(1:3, NA)))

# ---- Checks ----
cat("Person-waves (iwstat != 0):", nrow(pw), "\n")
resp <- filter(pw, iwstat == 1)
cat("Respondent waves:", nrow(resp), " missing state:", sum(is.na(resp$state)),
    " (adl missing:", sum(is.na(resp$adl5)), " cog missing:", sum(is.na(resp$cog)), ")\n")

tab <- resp |>
  filter(!is.na(state)) |>
  count(year, state) |>
  group_by(year) |> mutate(pct = round(100 * n / sum(n), 1)) |> ungroup() |>
  mutate(state = state_labels[state])
dir.create(here("output/tables"), recursive = TRUE, showWarnings = FALSE)
write_csv(tab, here("output/tables/02_state_by_wave.csv"))
tab |> select(-n) |> pivot_wider(names_from = state, values_from = pct) |> print(n = Inf)

cat("\nAge at interview among respondents: "); print(summary(resp$age))
cat("Persons ever coded dead (iwstat 5) who respond afterwards:\n")
bad <- pw |> arrange(hhidpn, wave) |> group_by(hhidpn) |>
  summarise(bad = any(cumsum(iwstat %in% c(5, 6)) > 0 & iwstat == 1)) |> summarise(sum(bad))
print(bad)

saveRDS(person, here("data/derived/person.rds"))
saveRDS(pw,     here("data/derived/person_wave.rds"))
