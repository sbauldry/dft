# Shared constants for the pipeline
waves     <- 4:16                       # RAND waves used (1998-2022)
wave_year <- function(w) 1998L + 2L * (w - 4L)
state_labels <- c("1 Functional", "2 Physical only", "3 Cognitive only", "4 Both", "5 Dead")
