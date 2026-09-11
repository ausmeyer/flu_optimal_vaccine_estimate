#!/usr/bin/env Rscript

# Cache Census cartographic state boundaries for manuscript maps.

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(tigris)
})

output_path <- "data/processed/us_states_2024_cartographic_boundaries.rds"
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

state_boundaries <- tigris::states(
  cb = TRUE,
  resolution = "20m",
  year = 2024,
  progress_bar = FALSE
) %>%
  filter(.data$STUSPS %in% c(state.abb, "DC")) %>%
  tigris::shift_geometry(position = "below") %>%
  select(
    state = NAME,
    state_abbreviation = STUSPS,
    GEOID,
    geometry
  )

if (nrow(state_boundaries) != 51L || any(!sf::st_is_valid(state_boundaries))) {
  stop("Expected valid boundaries for 50 states and the District of Columbia.")
}

saveRDS(state_boundaries, output_path)
message("Wrote ", output_path)
