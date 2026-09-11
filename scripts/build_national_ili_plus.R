#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("R/analysis_config.R")
source("R/ili_plus.R")

ilinet_path <- "data/processed/ilinet_national_all_age_primary_seasons.csv"
nrevss_path <- "data/processed/nrevss_clinical_labs_primary_seasons.csv"
output_path <- "data/processed/national_ili_plus_primary_seasons.csv"
diagnostic_path <- "outputs/diagnostics/data/national_ili_plus_coverage_by_season.csv"

ilinet <- read_csv(ilinet_path, show_col_types = FALSE) %>%
  filter_to_primary_seasons("national ILINet data")
nrevss <- read_csv(nrevss_path, show_col_types = FALSE) %>%
  filter_to_primary_seasons("national NREVSS data")

ili_plus <- prepare_national_ili_plus(ilinet, nrevss)
if (any(!is.finite(ili_plus$ili_plus_weighted)) ||
    any(!is.finite(ili_plus$ili_plus_unweighted)) ||
    any(ili_plus$ili_plus_weighted < 0 | ili_plus$ili_plus_weighted > 1) ||
    any(ili_plus$ili_plus_unweighted < 0 | ili_plus$ili_plus_unweighted > 1)) {
  stop("Constructed ILI+ values must be finite proportions in [0, 1].")
}

coverage <- ili_plus %>%
  summarise(
    n_weeks = n(),
    first_week_start = min(.data$week_start, na.rm = TRUE),
    last_week_start = max(.data$week_start, na.rm = TRUE),
    has_window_start = any(.data$mmwr_week == 36L),
    has_window_end = any(.data$mmwr_week == 22L),
    all_weeks_clinical = all(
      .data$nrevss_lab_series == "clinical_labs"
    ),
    nrevss_lab_series = paste(sort(unique(.data$nrevss_lab_series)), collapse = ","),
    .by = c("season", "season_start_year")
  ) %>%
  mutate(
    complete_burden_window =
      .data$has_window_start & .data$has_window_end,
    eligible_clinical_only =
      .data$complete_burden_window & .data$all_weeks_clinical
  ) %>%
  arrange(.data$season_start_year)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(diagnostic_path), recursive = TRUE, showWarnings = FALSE)
write_csv(ili_plus, output_path)
write_csv(coverage, diagnostic_path)

print(coverage, n = Inf)
message("Wrote national ILI+ input to ", output_path)
