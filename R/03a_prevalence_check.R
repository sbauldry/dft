# 03a_prevalence_check.R -- prevalence of dual functional state by age (comparison with past studies)
# Input:  data/derived/person_wave.rds, person.rds
# Output: output/tables/03a_prevalence_by_age.csv, 03a_prevalence_by_age_sex.csv ; output/figures/03a_prevalence_dual_functional.{pdf,png}
# Dual functional = 0 ADLs (adl5a == 0) and Langa-Weir not demented (cogfunction != 3).
# Respondents (iwstat == 1) age 50+ with observed state; unweighted; pooled across waves.
suppressPackageStartupMessages({library(tidyverse); library(here)})

sex <- readRDS(here("data/derived/person.rds")) |>
  transmute(hhidpn, sex = factor(ragender, 1:2, c("Men", "Women")))
pw <- readRDS(here("data/derived/person_wave.rds")) |>
  left_join(sex, by = "hhidpn") |>
  filter(iwstat == 1, !is.na(state), age >= 50) |>
  mutate(age_i = pmin(floor(age), 100L))   # single years, 100 = 100+

prev <- function(d, lab) d |> group_by(age_i) |>
  summarise(n = n(), dual = mean(state == 1), phys = mean(state == 2),
            cog = mean(state == 3), both = mean(state == 4), .groups = "drop") |>
  mutate(sample = lab)

out <- bind_rows(prev(pw, "1998-2022"), prev(filter(pw, year >= 2000), "2000-2022"))
dir.create(here("output/figures"), recursive = TRUE, showWarnings = FALSE)
write_csv(out |> mutate(across(c(dual, phys, cog, both), ~ round(.x, 4))),
          here("output/tables/03a_prevalence_by_age.csv"))

cat("Person-waves:", nrow(pw), " (2000+:", sum(pw$year >= 2000), ")\n")
cat("\nDual functional prevalence, 5-year age groups (1998-2022 | 2000-2022):\n")
grp <- function(d) d |> mutate(g = cut(age, c(50,55,60,65,70,75,80,85,90,95,Inf), right = FALSE,
                                       labels = c("50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95+"))) |>
  group_by(g) |> summarise(n = n(), dual = round(100 * mean(state == 1), 1))
print(left_join(grp(pw), grp(filter(pw, year >= 2000)), by = "g", suffix = c("_98", "_00")), n = Inf)
cat("\nSelected single years (1998-2022):\n")
print(out |> filter(sample == "1998-2022", age_i %in% c(50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)) |>
        mutate(dual = round(100 * dual, 1)) |> select(age_i, n, dual), n = Inf)

p <- ggplot(out, aes(age_i, 100 * dual, colour = sample)) +
  geom_line() + geom_point(aes(size = n), alpha = .5) +
  scale_size_area(max_size = 3, guide = "none") +
  labs(x = "Age", y = "% dual functional (0 ADLs, no dementia)", colour = NULL) +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(here("output/figures/03a_prevalence_dual_functional.pdf"), p, width = 6.5, height = 4.5)
ggsave(here("output/figures/03a_prevalence_dual_functional.png"), p, width = 6.5, height = 4.5, dpi = 300)

# ---- By sex ----
sx <- bind_rows(prev(filter(pw, sex == "Men"), "Men"), prev(filter(pw, sex == "Women"), "Women"))
sx <- sx |> rename(sex = sample)
write_csv(sx |> mutate(across(c(dual, phys, cog, both), ~ round(.x, 4))),
          here("output/tables/03a_prevalence_by_age_sex.csv"))
cat("\nSex missing:", sum(is.na(pw$sex)), " Men:", sum(pw$sex == "Men", na.rm = TRUE),
    " Women:", sum(pw$sex == "Women", na.rm = TRUE), "\n")
cat("\nDual functional % by sex, 5-year age groups (1998-2022):\n")
gs <- pw |> filter(!is.na(sex)) |>
  mutate(g = cut(age, c(50,55,60,65,70,75,80,85,90,95,Inf), right = FALSE,
                 labels = c("50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95+"))) |>
  group_by(g, sex) |> summarise(n = n(), dual = round(100 * mean(state == 1), 1), .groups = "drop") |>
  pivot_wider(names_from = sex, values_from = c(n, dual))
print(gs, n = Inf)
cat("\nSelected single years:\n")
print(sx |> filter(age_i %in% c(50, 60, 70, 80, 90, 95, 100)) |>
        mutate(dual = round(100 * dual, 1)) |> select(sex, age_i, n, dual) |>
        pivot_wider(names_from = sex, values_from = c(n, dual)), n = Inf)
p2 <- ggplot(sx, aes(age_i, 100 * dual, colour = sex)) +
  geom_line() + geom_point(aes(size = n), alpha = .5) +
  scale_size_area(max_size = 3, guide = "none") +
  labs(x = "Age", y = "% dual functional (0 ADLs, no dementia)", colour = NULL) +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(here("output/figures/03a_prevalence_dual_functional_sex.pdf"), p2, width = 6.5, height = 4.5)
ggsave(here("output/figures/03a_prevalence_dual_functional_sex.png"), p2, width = 6.5, height = 4.5, dpi = 300)
