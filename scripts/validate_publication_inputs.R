#!/usr/bin/env Rscript

# Fail before model fitting if freshly prepared public inputs do not match the
# prespecified publication population, seasons, and MMWR surveillance window.

suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
  library(readr)
})

source("R/package_requirements.R")
source("R/publication_spec.R")

expected_primary_seasons <- c(2010:2018, 2023:2025)
expected_excluded_seasons <- 2019:2022
expected_states <- sort(c(state.name, "District of Columbia"))
expected_age_groups <- sort(c("0-4", "5-24", "25-49", "50-64", "65+"))
expected_hhs_regions <- paste("Region", 1:10)
expected_clinical_only_seasons <- c(2016:2018, 2023:2025)

raw_paths <- c(
  "data/raw/ilinet_state_raw.csv",
  "data/raw/ilinet_hhs_raw.csv",
  "data/raw/ilinet_national_raw.csv",
  "data/raw/nrevss_who_raw.rds",
  "data/raw/cdc_ve_raw_tables.rds"
)
processed_paths <- c(
  state = "data/processed/ilinet_state_all_age_primary_seasons.csv",
  hhs = "data/processed/ilinet_hhs_age_primary_seasons.csv",
  national = "data/processed/ilinet_national_all_age_primary_seasons.csv",
  nrevss = "data/processed/nrevss_clinical_labs_primary_seasons.csv",
  population = "data/processed/state_population.csv",
  ve = "data/processed/cdc_ve_estimates.csv",
  boundaries = "data/processed/us_states_2024_cartographic_boundaries.rds"
)
input_paths <- c(raw_paths, unname(processed_paths))

missing_paths <- input_paths[!file.exists(input_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Publication input preparation did not create:\n",
    paste(missing_paths, collapse = "\n"),
    call. = FALSE
  )
}

assert_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(label, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

assert_exact_set <- function(observed, expected, label) {
  observed <- sort(unique(as.character(observed[!is.na(observed)])))
  expected <- sort(unique(as.character(expected)))
  if (!identical(observed, expected)) {
    stop(
      label, " does not match the publication specification. Observed: ",
      paste(observed, collapse = ", "), "; expected: ",
      paste(expected, collapse = ", "), ".",
      call. = FALSE
    )
  }
}

read_publication_csv <- function(path, required, label) {
  data <- read_csv(path, show_col_types = FALSE)
  if (nrow(data) == 0L) {
    stop(label, " is empty: ", path, call. = FALSE)
  }
  assert_columns(data, required, label)
  data
}

common_surveillance_columns <- c(
  "season", "season_start_year", "mmwr_year", "mmwr_week", "week_start"
)

validate_surveillance_scope <- function(data, label, delayed_start_seasons = integer()) {
  assert_columns(data, common_surveillance_columns, label)

  if (anyNA(data$season_start_year) || anyNA(data$mmwr_year) ||
      anyNA(data$mmwr_week)) {
    stop(label, " has missing season or MMWR identifiers.", call. = FALSE)
  }

  season_start_year <- as.integer(data$season_start_year)
  mmwr_year <- as.integer(data$mmwr_year)
  mmwr_week <- as.integer(data$mmwr_week)
  derived_start_year <- ifelse(mmwr_week >= 36L, mmwr_year, mmwr_year - 1L)
  expected_label <- paste0(
    season_start_year, "/", sprintf("%02d", (season_start_year + 1L) %% 100L)
  )

  assert_exact_set(
    season_start_year,
    expected_primary_seasons,
    paste0(label, " season start years")
  )
  if (any(season_start_year %in% expected_excluded_seasons)) {
    stop(label, " contains a prespecified excluded 2019/20-2022/23 season.", call. = FALSE)
  }
  if (any(!mmwr_week %in% c(36:53, 1:22))) {
    stop(label, " contains weeks outside the MMWR week 36-22 burden window.", call. = FALSE)
  }
  if (any(season_start_year != derived_start_year)) {
    stop(label, " has MMWR weeks assigned to the wrong influenza season.", call. = FALSE)
  }
  if (any(as.character(data$season) != expected_label)) {
    stop(label, " has an inconsistent season label.", call. = FALSE)
  }

  week_start <- as.Date(data$week_start)
  if (anyNA(week_start)) {
    stop(label, " has missing or invalid week_start dates.", call. = FALSE)
  }

  week_keys <- tibble::tibble(
    season_start_year = season_start_year,
    mmwr_year = mmwr_year,
    mmwr_week = mmwr_week,
    week_start = week_start
  ) %>%
    distinct()

  conflicting_dates <- week_keys %>%
    count(.data$season_start_year, .data$mmwr_year, .data$mmwr_week) %>%
    filter(.data$n != 1L)
  if (nrow(conflicting_dates) > 0L) {
    stop(label, " assigns more than one date to an MMWR week.", call. = FALSE)
  }

  coverage <- week_keys %>%
    arrange(.data$season_start_year, .data$week_start) %>%
    summarise(
      first_year = first(.data$mmwr_year),
      first_week = first(.data$mmwr_week),
      first_date = first(.data$week_start),
      last_year = last(.data$mmwr_year),
      last_week = last(.data$mmwr_week),
      last_date = last(.data$week_start),
      continuous = all(diff(.data$week_start) == 7),
      .by = "season_start_year"
    )

  permitted_first_week <- ifelse(
    coverage$season_start_year %in% delayed_start_seasons,
    coverage$first_week %in% c(36L, 40L),
    coverage$first_week == 36L
  )
  valid_end <- coverage$last_year == coverage$season_start_year + 1L &
    coverage$last_week == 22L
  valid_start_year <- coverage$first_year == coverage$season_start_year

  if (any(!permitted_first_week) || any(!valid_start_year) ||
      any(!valid_end) || any(!coverage$continuous)) {
    stop(
      label,
      " does not provide a continuous publication burden window. Seasons must ",
      "start at week 36 and end at week 22, except for explicitly permitted ",
      "week-40 starts in season start years: ",
      if (length(delayed_start_seasons)) paste(delayed_start_seasons, collapse = ", ") else "none",
      ".",
      call. = FALSE
    )
  }

  week_keys
}

assert_same_week_grid <- function(observed, expected, label) {
  by <- c("season_start_year", "mmwr_year", "mmwr_week", "week_start")
  missing <- anti_join(expected, observed, by = by)
  extra <- anti_join(observed, expected, by = by)
  if (nrow(missing) > 0L || nrow(extra) > 0L) {
    stop(label, " does not align exactly with the national ILINet week grid.", call. = FALSE)
  }
}

state_data <- read_publication_csv(
  processed_paths[["state"]],
  c(
    common_surveillance_columns, "state", "ilitotal", "total_patients",
    "unweighted_ili", "num_providers"
  ),
  "State ILINet input"
)
hhs_data <- read_publication_csv(
  processed_paths[["hhs"]],
  c(common_surveillance_columns, "hhs_region", "age_group", "ili_age_count", "total_patients"),
  "HHS age-specific ILINet input"
)
national_data <- read_publication_csv(
  processed_paths[["national"]],
  c(
    common_surveillance_columns, "state", "age_group", "ilitotal",
    "total_patients", "unweighted_ili", "weighted_ili"
  ),
  "National ILINet input"
)
nrevss_data <- read_publication_csv(
  processed_paths[["nrevss"]],
  c(
    common_surveillance_columns, "source_region_type", "source_region",
    "positive_tests", "total_tests", "nrevss_lab_series"
  ),
  "WHO/NREVSS input"
)

state_weeks <- validate_surveillance_scope(
  state_data, "State ILINet input", delayed_start_seasons = 2010L
)
hhs_weeks <- validate_surveillance_scope(hhs_data, "HHS age-specific ILINet input")
national_weeks <- validate_surveillance_scope(national_data, "National ILINet input")

national_endpoint <- national_weeks %>%
  arrange(.data$week_start) %>%
  slice_tail(n = 1L)
if (national_endpoint$mmwr_year[[1]] != 2026L ||
    national_endpoint$mmwr_week[[1]] != 22L ||
    national_endpoint$week_start[[1]] != as.Date("2026-05-31")) {
  stop(
    "National ILINet input must end at MMWR year 2026, week 22 ",
    "(week starting 2026-05-31).",
    call. = FALSE
  )
}

# Public state reporting began at week 40 of 2010. National and HHS series
# cover the September weeks as well.
expected_state_weeks <- national_weeks %>%
  filter(!(.data$season_start_year == 2010L & .data$mmwr_year == 2010L &
    .data$mmwr_week < 40L))
assert_same_week_grid(state_weeks, expected_state_weeks, "State ILINet input")
assert_same_week_grid(hhs_weeks, national_weeks, "HHS age-specific ILINet input")

assert_exact_set(state_data$state, expected_states, "State ILINet jurisdictions")
state_duplicates <- state_data %>%
  count(.data$season_start_year, .data$state, .data$mmwr_year, .data$mmwr_week) %>%
  filter(.data$n != 1L)
state_week_counts <- state_data %>%
  summarise(
    n_jurisdictions = n_distinct(.data$state),
    .by = c("season_start_year", "mmwr_year", "mmwr_week")
  )
if (nrow(state_duplicates) > 0L || any(state_week_counts$n_jurisdictions != 51L)) {
  stop("State ILINet must contain one row for each of 51 jurisdictions per week.", call. = FALSE)
}
state_missing_outpatient_counts <- is.na(state_data$ilitotal) &
  is.na(state_data$total_patients)
state_mismatched_outpatient_counts <- xor(
  is.na(state_data$ilitotal),
  is.na(state_data$total_patients)
)
state_observed_outpatient_counts <- !state_missing_outpatient_counts &
  !state_mismatched_outpatient_counts
state_invalid_outpatient_counts <- state_observed_outpatient_counts & (
  !is.finite(state_data$ilitotal) |
    !is.finite(state_data$total_patients) |
    state_data$ilitotal < 0 |
    state_data$total_patients <= 0 |
    state_data$ilitotal > state_data$total_patients
)
state_invalid_summaries <-
  (state_missing_outpatient_counts & (
    !is.na(state_data$unweighted_ili) |
      !is.na(state_data$num_providers)
  )) |
  (state_observed_outpatient_counts & (
    !is.finite(state_data$unweighted_ili) |
      !is.finite(state_data$num_providers) |
      state_data$num_providers <= 0
  ))
if (any(state_mismatched_outpatient_counts) ||
    any(state_invalid_outpatient_counts) ||
    any(state_invalid_summaries)) {
  stop(
    paste0(
      "State ILINet must have either paired missing outpatient counts and ",
      "summaries or valid reported counts, percentages, and provider totals."
    ),
    call. = FALSE
  )
}
state_missing_rows <- state_data[state_missing_outpatient_counts, , drop = FALSE]
state_missing_jurisdictions <- sort(unique(state_missing_rows$state))

assert_exact_set(hhs_data$age_group, expected_age_groups, "HHS ILINet age groups")
assert_exact_set(hhs_data$hhs_region, expected_hhs_regions, "HHS ILINet regions")
hhs_duplicates <- hhs_data %>%
  count(
    .data$season_start_year, .data$hhs_region, .data$age_group,
    .data$mmwr_year, .data$mmwr_week
  ) %>%
  filter(.data$n != 1L)
hhs_grid_columns <- c("season_start_year", "mmwr_year", "mmwr_week")
expected_hhs_grid <- tidyr::crossing(
  national_weeks %>% select(all_of(hhs_grid_columns)),
  hhs_region = expected_hhs_regions,
  age_group = expected_age_groups
)
observed_hhs_grid <- hhs_data %>%
  transmute(
    season_start_year = as.integer(.data$season_start_year),
    mmwr_year = as.integer(.data$mmwr_year),
    mmwr_week = as.integer(.data$mmwr_week),
    hhs_region = as.character(.data$hhs_region),
    age_group = as.character(.data$age_group)
  )
hhs_grid_keys <- c(hhs_grid_columns, "hhs_region", "age_group")
missing_hhs_grid <- anti_join(
  expected_hhs_grid,
  observed_hhs_grid,
  by = hhs_grid_keys
)
extra_hhs_grid <- anti_join(
  observed_hhs_grid,
  expected_hhs_grid,
  by = hhs_grid_keys
)
if (nrow(hhs_duplicates) > 0L ||
    nrow(missing_hhs_grid) > 0L ||
    nrow(extra_hhs_grid) > 0L) {
  stop(
    paste0(
      "HHS age-specific ILINet must contain the complete weekly cross-product ",
      "of 10 HHS regions and five age groups, with exactly one row per cell."
    ),
    call. = FALSE
  )
}
if (any(!is.finite(hhs_data$ili_age_count)) ||
    any(!is.finite(hhs_data$total_patients)) ||
    any(hhs_data$ili_age_count < 0) || any(hhs_data$total_patients <= 0) ||
    any(hhs_data$ili_age_count > hhs_data$total_patients)) {
  stop("HHS age-specific ILINet has invalid outpatient counts.", call. = FALSE)
}

if (!identical(unique(national_data$state), "US") ||
    !identical(unique(national_data$age_group), "all")) {
  stop("National ILINet must contain only the US, all-age series.", call. = FALSE)
}
national_duplicates <- national_data %>%
  count(.data$season_start_year, .data$mmwr_year, .data$mmwr_week) %>%
  filter(.data$n != 1L)
if (nrow(national_duplicates) > 0L ||
    any(!is.finite(national_data$ilitotal)) ||
    any(!is.finite(national_data$total_patients)) ||
    any(!is.finite(national_data$unweighted_ili)) ||
    any(!is.finite(national_data$weighted_ili)) ||
    any(national_data$ilitotal < 0) || any(national_data$total_patients <= 0) ||
    any(national_data$ilitotal > national_data$total_patients)) {
  stop("National ILINet has duplicate weeks or invalid values.", call. = FALSE)
}

national_nrevss <- nrevss_data %>%
  filter(.data$source_region_type == "national")
if (nrow(national_nrevss) == 0L) {
  stop("WHO/NREVSS input does not contain the required national series.", call. = FALSE)
}
nrevss_weeks <- validate_surveillance_scope(national_nrevss, "National WHO/NREVSS input")
assert_same_week_grid(nrevss_weeks, national_weeks, "National WHO/NREVSS input")
nrevss_duplicates <- national_nrevss %>%
  count(.data$season_start_year, .data$mmwr_year, .data$mmwr_week) %>%
  filter(.data$n != 1L)
if (nrow(nrevss_duplicates) > 0L ||
    any(!is.finite(national_nrevss$positive_tests)) ||
    any(!is.finite(national_nrevss$total_tests)) ||
    any(national_nrevss$positive_tests < 0) ||
    any(national_nrevss$total_tests <= 0) ||
    any(national_nrevss$positive_tests > national_nrevss$total_tests)) {
  stop("National WHO/NREVSS has duplicate weeks or invalid laboratory counts.", call. = FALSE)
}

clinical_coverage <- national_nrevss %>%
  arrange(.data$week_start) %>%
  summarise(
    all_clinical = all(.data$nrevss_lab_series == "clinical_labs"),
    starts_week_36 = first(.data$mmwr_week) == 36L,
    ends_week_22 = last(.data$mmwr_week) == 22L,
    .by = "season_start_year"
  ) %>%
  filter(.data$all_clinical, .data$starts_week_36, .data$ends_week_22)
assert_exact_set(
  clinical_coverage$season_start_year,
  expected_clinical_only_seasons,
  "Complete clinical-laboratory-only WHO/NREVSS seasons"
)

population <- read_publication_csv(
  processed_paths[["population"]],
  c("state", "population", "source_year"),
  "State population input"
)
assert_exact_set(population$state, expected_states, "Population jurisdictions")
if (nrow(population) != 51L || anyDuplicated(population$state) ||
    !identical(sort(unique(as.integer(population$source_year))), 2025L) ||
    any(!is.finite(population$population)) || any(population$population <= 0)) {
  stop("Population input must contain positive 2025 estimates for 50 states and DC.", call. = FALSE)
}

ve <- read_publication_csv(
  processed_paths[["ve"]],
  c(
    "season", "season_start_year", "ve_age_group", "ve", "lower_95", "upper_95",
    "log_or_mean", "log_or_sd", "source_id", "citation_key", "source_url",
    "source_age_group", "estimate_status", "source_table"
  ),
  "Vaccine-effectiveness input"
)
expected_ve_coverage <- publication_ve_pool_coverage()
expected_ve_seasons <- sort(unique(expected_ve_coverage$season_start_year))
assert_exact_set(ve$season_start_year, expected_ve_seasons, "Vaccine-effectiveness seasons")
assert_exact_set(
  ve$ve_age_group,
  unique(expected_ve_coverage$ve_age_group),
  "Vaccine-effectiveness age groups"
)
ve_missing_pairs <- anti_join(
  expected_ve_coverage, ve, by = c("season_start_year", "ve_age_group")
)
ve_extra_pairs <- anti_join(
  ve, expected_ve_coverage, by = c("season_start_year", "ve_age_group")
)
ve_duplicates <- ve %>%
  count(.data$season_start_year, .data$ve_age_group) %>%
  filter(.data$n != 1L)
if (nrow(ve) != nrow(expected_ve_coverage) ||
    nrow(ve_missing_pairs) > 0L || nrow(ve_extra_pairs) > 0L ||
    nrow(ve_duplicates) > 0L ||
    any(!is.finite(ve$ve)) || any(!is.finite(ve$lower_95)) ||
    any(!is.finite(ve$upper_95)) ||
    any(ve$lower_95 > ve$ve | ve$ve > ve$upper_95 | ve$upper_95 >= 1) ||
    any(!is.finite(ve$log_or_mean)) || any(!is.finite(ve$log_or_sd)) ||
    any(ve$log_or_sd <= 0)) {
  stop("Vaccine-effectiveness input is incomplete or contains invalid estimates.", call. = FALSE)
}
ve_provenance_columns <- c(
  "source_id", "citation_key", "source_url", "source_age_group",
  "estimate_status", "source_table"
)
if (any(vapply(ve[ve_provenance_columns], function(values) {
      anyNA(values) || any(!nzchar(trimws(values)))
    }, logical(1))) ||
    any(!grepl("^https://", ve$source_url)) ||
    any(ve$estimate_status != ifelse(
      ve$season_start_year == 2025L, "preliminary", "final"
    ))) {
  stop("Vaccine-effectiveness input has incomplete or incorrect source provenance.", call. = FALSE)
}

# The final 2024/25 US Flu VE Network report supersedes the interim input.
ve_2024_overall <- ve[ve$season_start_year == 2024L & ve$ve_age_group == "all", ]
if (nrow(ve_2024_overall) != 1L ||
    any(abs(c(ve_2024_overall$ve, ve_2024_overall$lower_95,
              ve_2024_overall$upper_95) - c(0.33, 0.24, 0.41)) > 1e-12) ||
    ve_2024_overall$citation_key != "chung2026VE") {
  stop("The 2024/25 overall VE input must use the final Chung report.", call. = FALSE)
}

boundaries <- readRDS(processed_paths[["boundaries"]])
if (!is.data.frame(boundaries)) {
  stop("State boundary input is not a data frame.", call. = FALSE)
}
assert_columns(boundaries, c("state", "state_abbreviation", "GEOID"), "State boundary input")
assert_exact_set(boundaries$state, expected_states, "Boundary jurisdictions")
if (nrow(boundaries) != 51L || anyDuplicated(boundaries$state) ||
    anyDuplicated(boundaries$GEOID)) {
  stop("State boundary input must contain one boundary for each of 50 states and DC.", call. = FALSE)
}

file_records <- lapply(input_paths, function(path) {
  info <- file.info(path)
  list(
    path = path,
    md5 = unname(tools::md5sum(path)),
    bytes = unname(info$size),
    local_modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE)
  )
})
installed_version_or_na <- function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    NA_character_
  }
}
installed_package_versions <- vapply(
  publication_required_packages,
  installed_version_or_na,
  character(1)
)
additional_model_package_versions <- vapply(
  c("Matrix", "nlme"),
  installed_version_or_na,
  character(1)
)
os_fields <- c("sysname", "release", "version", "machine")
operating_system <- as.list(as.character(Sys.info()[os_fields]))
names(operating_system) <- os_fields
rng_kind <- as.list(RNGkind())
names(rng_kind) <- c("kind", "normal_kind", "sample_kind")
session <- sessionInfo()
blas_path <- if (!is.null(session$BLAS)) session$BLAS else extSoftVersion()[["BLAS"]]
lapack_path <- if (!is.null(session$LAPACK)) session$LAPACK else La_library()
sf_external_versions <- if (requireNamespace("sf", quietly = TRUE)) {
  as.list(sf::sf_extSoftVersion())
} else {
  list()
}

validation_record <- list(
  validated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  software = list(
    r_version = R.version.string,
    r_version_numeric = as.character(getRversion()),
    r_platform = R.version$platform,
    os_type = .Platform$OS.type,
    operating_system = operating_system,
    locale = Sys.getlocale(),
    time_zone = Sys.timezone(),
    rng_kind = rng_kind,
    blas = list(library = basename(blas_path)),
    lapack = list(
      library = basename(lapack_path),
      version = as.character(La_version())
    ),
    cairo_available = unname(capabilities("cairo")),
    sf_external_versions = sf_external_versions,
    direct_package_versions = as.list(installed_package_versions),
    additional_model_package_versions = as.list(
      additional_model_package_versions
    )
  ),
  constraints = list(
    primary_season_start_years = expected_primary_seasons,
    excluded_season_start_years = expected_excluded_seasons,
    burden_window = "MMWR weeks 36 through 22, inclusive",
    delayed_start_exception = "State ILINet begins at week 40 in 2010/11; national and HHS series begin at week 36",
    observed_endpoint = list(mmwr_year = 2026L, mmwr_week = 22L),
    jurisdictions = 51L,
    hhs_regions = 10L,
    hhs_age_groups = expected_age_groups,
    ve_season_start_years = expected_ve_seasons,
    ve_season_age_coverage = expected_ve_coverage,
    state_population_year = 2025L,
    boundary_year = 2024L,
    clinical_only_nrevss_season_start_years = expected_clinical_only_seasons
  ),
  observed_missingness = list(
    state_missing_outpatient_count_rows = nrow(state_missing_rows),
    state_missing_outpatient_count_jurisdictions = state_missing_jurisdictions,
    national_zero_denominator_rows = sum(national_data$total_patients == 0)
  ),
  files = file_records
)
validation_path <- "outputs/diagnostics/data/publication_input_validation.json"
dir.create(dirname(validation_path), recursive = TRUE, showWarnings = FALSE)
write_json(
  validation_record,
  validation_path,
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null"
)

message(
  "Publication inputs validated: 12 seasons, MMWR weeks 36-22, ",
  "51 jurisdictions, endpoint 2026 week 22, and aligned national NREVSS coverage."
)
message("Wrote ", validation_path)
