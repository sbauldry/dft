# 06_figures.R -- figures from the life table summaries (main prior, sd 5 by default).
# Usage: Rscript R/06_figures.R [prior_sd=5]
# Reads  output/tables/lifetable_summary_<sex>_psd<sd>.csv
# Writes output/figures/06_*.pdf and .png (overwrites figures from this script only).
# Colours: validated categorical palette slots 1-4 (adjacent-pair CVD safe), fixed state order.
suppressPackageStartupMessages({library(tidyverse); library(here)})
source(here("R/functions/constants.R"))

args <- commandArgs(trailingOnly = TRUE)
prior_sd <- if (length(args) >= 1) args[1] else "5"
max_age <- 100                                   # ages above 100 have little data support; table runs to 110+

tab <- map_dfr(c("men", "women"), \(s)
  read_csv(here("output/tables", sprintf("lifetable_summary_%s_psd%s.csv", s, prior_sd)), show_col_types = FALSE)) |>
  mutate(sex = factor(sex, c("men", "women"), c("Men", "Women")),
         state = factor(state, c(state_labels[1:4], "Total"))) |>
  filter(age <= max_age)

state_cols <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#e87ba4"), state_labels[1:4])
sex_cols   <- c(Men = "#2a78d6", Women = "#eb6834")
ink <- "#0b0b0b"; ink2 <- "#52514e"

theme_set(theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
        axis.text = element_text(colour = ink2), text = element_text(colour = ink),
        strip.text = element_text(face = "bold", hjust = 0), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12), plot.caption = element_text(colour = ink2, hjust = 0)))

save_fig <- function(p, name, w, h) {
  walk(c("pdf", "png"), \(ext) ggsave(here("output/figures", sprintf("06_%s.%s", name, ext)), p,
                                      width = w, height = h, dpi = 300, bg = "white"))
}
cap <- sprintf("Lines: posterior mean; bands: 84%% credible interval. Population-based radix at age 50; prior N(0, %s^2).", prior_sd)

# 1. State expectancies by age (rows: state, columns: sex)
d1 <- tab |> filter(quantity == "expectancy", radix == "population")
p1 <- ggplot(d1, aes(age, mean, colour = state, fill = state)) +
  geom_ribbon(aes(ymin = lo84, ymax = hi84), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 0.7) +
  facet_grid(state ~ sex, scales = "free_y") +
  scale_colour_manual(values = state_cols) + scale_fill_manual(values = state_cols) +
  scale_x_continuous(breaks = seq(50, 100, 10)) + expand_limits(y = 0) +
  labs(x = "Age", y = "Expected years in state", colour = NULL, fill = NULL,
       title = "Expected years in each state, by age", caption = cap) +
  guides(colour = "none", fill = "none")
save_fig(p1, "state_expectancy", 6.5, 8)

# 2. Total life expectancy, men vs women
d2 <- tab |> filter(quantity == "total_expectancy", radix == "population")
p2 <- ggplot(d2, aes(age, mean, colour = sex, fill = sex)) +
  geom_ribbon(aes(ymin = lo84, ymax = hi84), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = sex_cols) + scale_fill_manual(values = sex_cols) +
  scale_x_continuous(breaks = seq(50, 100, 10)) + expand_limits(y = 0) +
  labs(x = "Age", y = "Total life expectancy (years)", colour = NULL, fill = NULL,
       title = "Total life expectancy, by age and sex", caption = cap)
save_fig(p2, "total_expectancy", 6, 4)

# 3. Percent of remaining life in each state (stacked, posterior means)
d3 <- tab |> filter(quantity == "pct_remaining", radix == "population", state != "Total")
p3 <- ggplot(d3, aes(age, mean, fill = state)) +
  geom_area(colour = "white", linewidth = 0.3) +
  facet_wrap(~ sex) +
  scale_fill_manual(values = state_cols) +
  scale_x_continuous(breaks = seq(50, 100, 10)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Age", y = "% of remaining life", fill = NULL,
       title = "Share of remaining life spent in each state",
       caption = sprintf("Posterior means (computed draw by draw). Population-based radix; prior N(0, %s^2).", prior_sd))
save_fig(p3, "pct_remaining", 7, 4)

# 4. Status-based radix: expectancy at 50 by starting state
d4 <- tab |> filter(age == 50, radix %in% paste0("s", 1:4), state != "Total",
                    quantity %in% c("expectancy")) |>
  mutate(start = factor(radix, paste0("s", 1:4), paste("Start in", state_labels[1:4])))
p4 <- ggplot(d4, aes(mean, start, fill = state)) +
  geom_col(position = position_stack(reverse = TRUE), colour = "white", linewidth = 0.4, width = 0.7) +
  facet_wrap(~ sex, ncol = 1) +
  scale_fill_manual(values = state_cols) + scale_y_discrete(limits = rev) +
  labs(x = "Expected years from age 50", y = NULL, fill = NULL,
       title = "Expected years by state at age 50, by starting state",
       caption = "Status-based radix: everyone starts in the stated state at age 50. Posterior means.")
save_fig(p4, "status_radix_age50", 7, 5.5)
cat("saved figures to", here("output/figures"), "\n")
