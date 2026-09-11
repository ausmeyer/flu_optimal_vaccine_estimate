library(dplyr)
library(purrr)
library(readr)

source("R/analysis_config.R")
source("R/fluview_download.R")
source("R/lab_surveillance.R")

primary_season_start_years <- primary_season_start_years()
download_years <- fluview_download_seasons(primary_season_start_years)
nrevss_regions <- Sys.getenv("NREVSS_REGIONS", "national,hhs,state") %>%
  strsplit(",", fixed = TRUE) %>%
  unlist(use.names = FALSE) %>%
  trimws()

dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/diagnostics/data", showWarnings = FALSE, recursive = TRUE)

redownload <- tolower(Sys.getenv("REDOWNLOAD_NREVSS", "false")) %in% c("1", "true", "yes")
raw_path <- "data/raw/nrevss_who_raw.rds"

raw <- if (!redownload && file.exists(raw_path)) readRDS(raw_path) else NULL
raw_covers_window <- !is.null(raw) && all(nrevss_regions %in% names(raw)) &&
  all(vapply(raw[nrevss_regions], function(region_result) {
    table_names <- grep(
      "(clinical_labs|combined_prior_to_2015_16)$", names(region_result), value = TRUE
    )
    weeks <- bind_rows(lapply(region_result[table_names], function(x) {
      select(x, any_of(c("year", "week")))
    }))
    fluview_raw_covers_burden_window(weeks, primary_season_start_years)
  }, logical(1)))

if (!raw_covers_window) {
  message(
    "Downloading WHO/NREVSS seasons beginning in: ",
    paste(download_years, collapse = ", ")
  )
  raw <- setNames(
    lapply(
      nrevss_regions,
      function(region) fluview_download_who_nrevss(region = region, years = download_years)
    ),
    nrevss_regions
  )
  saveRDS(raw, raw_path)
}

clinical <- purrr::imap_dfr(raw, function(region_result, region_label) {
  clinical_name <- grep("clinical_labs$", names(region_result), value = TRUE)
  combined_name <- grep("combined_prior_to_2015_16$", names(region_result), value = TRUE)

  pieces <- list()
  if (length(combined_name) > 0L) {
    pieces[["combined_prior_to_2015_16"]] <- normalise_nrevss_clinical_labs(
      region_result[[combined_name[[1]]]],
      region_label = region_label
    ) %>%
      mutate(nrevss_lab_series = "combined_prior_to_2015_16")
  }
  if (length(clinical_name) > 0L) {
    pieces[["clinical_labs"]] <- normalise_nrevss_clinical_labs(
      region_result[[clinical_name[[1]]]],
      region_label = region_label
    ) %>%
      mutate(nrevss_lab_series = "clinical_labs")
  }

  if (length(pieces) == 0L) {
    return(tibble::tibble())
  }

  bind_rows(pieces)
}) %>%
  filter(
    .data$season_start_year %in% primary_season_start_years,
    .data$in_primary_burden_window
  ) %>%
  mutate(
    nrevss_lab_series_priority = if_else(.data$nrevss_lab_series == "clinical_labs", 1L, 2L)
  ) %>%
  arrange(
    .data$source_region_type,
    .data$source_region,
    .data$mmwr_year,
    .data$mmwr_week,
    .data$nrevss_lab_series_priority
  ) %>%
  distinct(
    .data$source_region_type,
    .data$source_region,
    .data$mmwr_year,
    .data$mmwr_week,
    .keep_all = TRUE
  ) %>%
  select(
    -"nrevss_lab_series_priority"
  )

if (nrow(clinical) == 0L) {
  stop("No WHO/NREVSS clinical or pre-2015 combined lab rows were normalized for the requested seasons.")
}

write_csv(clinical, "data/processed/nrevss_clinical_labs_primary_seasons.csv")

coverage <- clinical %>%
  group_by(.data$source_region_type, .data$season, .data$season_start_year) %>%
  summarise(
    n_locations = n_distinct(.data$state),
    first_week_start = suppressWarnings(min(.data$week_start, na.rm = TRUE)),
    last_week_start = suppressWarnings(max(.data$week_start, na.rm = TRUE)),
    observed_weeks = n_distinct(paste(.data$mmwr_year, .data$mmwr_week)),
    complete_burden_window = any(.data$mmwr_week == 36L) &&
      any(.data$mmwr_week == 22L) &&
      all(diff(sort(unique(.data$week_start))) == 7),
    .groups = "drop"
  )

write_csv(coverage, "outputs/diagnostics/data/nrevss_clinical_lab_coverage_by_season.csv")
print(coverage)
message("Wrote WHO/NREVSS clinical lab processed file.")
