#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

year <- as.integer(Sys.getenv("STATE_POPULATION_YEAR", "2025"))
if (year != 2025L) {
  stop("This reproducible download currently supports STATE_POPULATION_YEAR=2025.")
}
url <- paste0(
  "https://www2.census.gov/programs-surveys/popest/datasets/",
  "2020-2025/state/totals/NST-EST2025-ALLDATA.csv"
)

raw <- read_csv(url, show_col_types = FALSE)
population <- raw %>%
  transmute(
    state = .data$NAME,
    population = .data$POPESTIMATE2025,
    census_state_fips = sprintf("%02d", as.integer(.data$STATE))
  )

population <- population %>%
  mutate(
    population = as.numeric(.data$population),
    source_year = year,
    source_dataset = "Census annual state population estimates, 2020-2025 vintage",
    accessed_date = as.character(Sys.Date())
  ) %>%
  filter(.data$state %in% c(state.name, "District of Columbia")) %>%
  arrange(.data$state)

if (nrow(population) != 51L || anyNA(population$population)) {
  stop("The Census response did not contain complete population data for 50 states and DC.")
}

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
write_csv(population, "data/processed/state_population.csv")
message("Wrote state population weights to data/processed/state_population.csv")
