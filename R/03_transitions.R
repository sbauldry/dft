# 03_transitions.R -- person-interval file and transition checks
# Input:  data/derived/person.rds, person_wave.rds
# Output: data/derived/intervals.rds ; output/tables/03_*.csv
# Rules (docs/decisions.md): interval starts at an interview with observed state and age >= 50;
# ends in death (iwstat == 5 at next wave) or a next-wave interview with observed state.
# Anything else (nonresponse, dropped, missing state) is censored/excluded.
suppressPackageStartupMessages({library(tidyverse); library(here)})
source(here("R/functions/constants.R"))

person <- readRDS(here("data/derived/person.rds"))
pw     <- readRDS(here("data/derived/person_wave.rds"))

# next-wave info attached to each starting interview
nxt <- pw |>
  transmute(hhidpn, wave = wave - 1L, next_iwstat = iwstat, next_state = state,
            next_state_cind = state_cind, next_iwend = iwend)

start <- pw |>
  filter(iwstat == 1, wave <= max(waves) - 1L) |>
  left_join(nxt, by = c("hhidpn", "wave")) |>
  mutate(next_iwstat = replace_na(next_iwstat, 0L))   # 0 = no record at next wave

# ---- Sample flow ----
flow <- start |>
  mutate(reason = case_when(
    age < 50                          ~ "a. start age < 50",
    is.na(state)                      ~ "b. missing start state",
    next_iwstat == 5                  ~ "included: ends in death",
    next_iwstat == 1 & !is.na(next_state) ~ "included: ends in observed state",
    next_iwstat == 1                  ~ "c. missing end state",
    next_iwstat == 4                  ~ "d. next wave nonresponse (alive)",
    next_iwstat == 7                  ~ "e. next wave dropped from sample",
    TRUE                              ~ "f. other (no record / died after gap)")) |>
  count(reason)
print(flow, n = Inf)
write_csv(flow, here("output/tables/03_sample_flow.csv"))

iv <- start |>
  filter(age >= 50, !is.na(state), next_iwstat == 5 | (next_iwstat == 1 & !is.na(next_state))) |>
  mutate(end_state = if_else(next_iwstat == 5, 5L, next_state),
         end_state_cind = if_else(next_iwstat == 5, 5L, next_state_cind),
         interval_months = as.numeric(as.Date(next_iwend, origin = "1960-01-01") -
                                      as.Date(iwend, origin = "1960-01-01")) / 30.4375) |>
  rename(start_state = state, start_state_cind = state_cind) |>
  left_join(person, by = "hhidpn") |>
  mutate(female = as.integer(ragender == 2),
         race_eth = case_when(rahispan == 1 ~ "Hispanic",
                              raracem == 1 ~ "NH White",
                              raracem == 2 ~ "NH Black",
                              raracem == 3 ~ "NH Other"),
         educ = case_when(raeduc == 1 ~ "<HS", raeduc %in% 2:3 ~ "HS/GED",
                          raeduc == 4 ~ "Some college", raeduc == 5 ~ "College+"),
         birth_year = rabyear)

# ---- Covariate missingness, then drop ----
cat("\nIntervals before covariate drop:", nrow(iv), " persons:", n_distinct(iv$hhidpn), "\n")
cat("Missing female/race_eth/educ/birth_year:",
    sum(is.na(iv$female)), sum(is.na(iv$race_eth)), sum(is.na(iv$educ)), sum(is.na(iv$birth_year)), "\n")
iv <- iv |> filter(!is.na(female), !is.na(race_eth), !is.na(educ), !is.na(birth_year)) |>
  mutate(race_eth = factor(race_eth, c("NH White", "NH Black", "Hispanic", "NH Other")),
         educ = factor(educ, c("<HS", "HS/GED", "Some college", "College+"))) |>
  select(hhidpn, wave, year, start_state, end_state, start_state_cind, end_state_cind,
         age, female, race_eth, educ, birth_year, hacohort, interval_months, proxy)
cat("Final intervals:", nrow(iv), " persons:", n_distinct(iv$hhidpn), "\n")

# ---- Checks ----
stopifnot(all(iv$start_state %in% 1:4), all(iv$end_state %in% 1:5), all(iv$age >= 50))
tt <- iv |> count(start_state, end_state) |>
  complete(start_state = 1:4, end_state = 1:5, fill = list(n = 0)) |>
  pivot_wider(names_from = end_state, values_from = n, names_prefix = "to_")
cat("\nTransition counts (rows = start state):\n"); print(tt)
write_csv(tt, here("output/tables/03_transition_counts.csv"))

tt_c <- iv |> count(start_state_cind, end_state_cind) |>
  complete(start_state_cind = 1:4, end_state_cind = 1:5, fill = list(n = 0)) |>
  pivot_wider(names_from = end_state_cind, values_from = n, names_prefix = "to_")
cat("\nSensitivity (CIND + dementia) transition counts:\n"); print(tt_c)
write_csv(tt_c, here("output/tables/03_transition_counts_cind.csv"))

cat("\nRow percentages (main):\n")
print(iv |> count(start_state, end_state) |> group_by(start_state) |>
        mutate(pct = round(100 * n / sum(n), 1)) |> ungroup() |>
        select(-n) |> pivot_wider(names_from = end_state, values_from = pct, names_prefix = "to_"))

cat("\nIntervals by sex:\n"); print(count(iv, female))
by_wave <- iv |> count(year, end_state) |> pivot_wider(names_from = end_state, values_from = n, values_fill = 0)
write_csv(by_wave, here("output/tables/03_end_state_by_year.csv"))
cat("\nIntervals by start year:\n"); print(count(iv, year) |> deframe())
cat("\nAge at start of interval:\n"); print(summary(iv$age))

dir.create(here("data/derived"), showWarnings = FALSE)
saveRDS(iv, here("data/derived/intervals.rds"))
