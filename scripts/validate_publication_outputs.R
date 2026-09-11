#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
  library(readr)
})

source("R/analysis_config.R")
source("R/ili_plus.R")
source("R/package_requirements.R")
source("R/publication_spec.R")
source("R/ve_draws_from_estimates.R")

expected_draws <- 5000L
expected_primary_seasons <- c(2010:2018, 2023:2025)
validate_artifact_hygiene <- env_flag("VALIDATE_ARTIFACT_HYGIENE", FALSE)
primary_prefixes <- publication_primary_prefixes()
sensitivity_prefixes <- publication_sensitivity_prefixes()
loso_prefixes <- publication_loso_prefixes(expected_primary_seasons)

publication_raw_paths <- file.path(
  "data/raw",
  c(
    "cdc_ve_raw_tables.rds",
    "ilinet_hhs_raw.csv",
    "ilinet_national_raw.csv",
    "ilinet_state_raw.csv",
    "nrevss_who_raw.rds"
  )
)
publication_processed_paths <- file.path(
  "data/processed",
  c(
    "cdc_ve_estimates.csv",
    "ilinet_hhs_age_primary_seasons.csv",
    "ilinet_national_all_age_primary_seasons.csv",
    "ilinet_state_all_age_primary_seasons.csv",
    "national_ili_plus_primary_seasons.csv",
    "nrevss_clinical_labs_primary_seasons.csv",
    "state_population.csv",
    "us_states_2024_cartographic_boundaries.rds"
  )
)
diagnostic_paths <- c(
  "outputs/diagnostics/calendar/calendar_guardrail_week_wrap.csv",
  "outputs/diagnostics/data/cdc_ve_estimate_coverage.csv",
  "outputs/diagnostics/data/ilinet_state_coverage_by_season.csv",
  "outputs/diagnostics/data/national_ili_plus_coverage_by_season.csv",
  "outputs/diagnostics/data/nrevss_clinical_lab_coverage_by_season.csv",
  "outputs/diagnostics/data/publication_input_validation.json",
  "outputs/diagnostics/gam_fit/gam_fit_diagnostics_combined.csv",
  "outputs/diagnostics/gam_fit/gam_fit_diagnostics_flags.csv",
  "outputs/diagnostics/gam_fit/gam_fit_diagnostics_rollup.csv",
  "outputs/diagnostics/held_out_validation/held_out_season_decision_validation.csv",
  "outputs/diagnostics/held_out_validation/held_out_season_decision_validation_summary.csv",
  "outputs/diagnostics/monte_carlo/monte_carlo_exact_argmax_split_audit.csv",
  "outputs/diagnostics/monte_carlo/monte_carlo_exact_argmax_split_rollup.csv",
  "outputs/diagnostics/monte_carlo/monte_carlo_probability_precision.csv",
  "outputs/diagnostics/monte_carlo/monte_carlo_probability_precision_rollup.csv"
)

assert_no_unexpected_files <- function(root, allowed_paths, label) {
  observed_paths <- list.files(
    root,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )
  unexpected_paths <- setdiff(observed_paths, allowed_paths)
  if (length(unexpected_paths) > 0L) {
    stop(
      "Unexpected ", label, " artifact(s):\n",
      paste(unexpected_paths, collapse = "\n")
    )
  }
}

table_suffixes <- c(
  "_analysis_manifest.json",
  "_optimal_week_draws.csv",
  "_optimal_week_distribution.csv",
  "_optimal_week_intervals.csv",
  "_week_cutoff_probabilities.csv",
  "_regret_curve.csv",
  "_near_optimal_curve.csv",
  "_ve_waning_draws.csv"
)
primary_paths <- unlist(lapply(
  primary_prefixes,
  function(prefix) file.path(
    "outputs/primary/tables",
    paste0(prefix, table_suffixes)
  )
))
sensitivity_paths <- unlist(lapply(
  sensitivity_prefixes,
  function(prefix) file.path(
    "outputs/sensitivity",
    prefix,
    "tables",
    paste0(prefix, table_suffixes)
  )
))
loso_table_paths <- unlist(lapply(
  loso_prefixes,
  function(prefix) file.path(
    "outputs/sensitivity",
    prefix,
    "tables",
    paste0(prefix, table_suffixes)
  )
))
final_figure_paths <- file.path(
  "outputs/primary/final_figures",
  c(
    "figure_1_epidemic_timing_and_model_fit.pdf",
    "figure_1_epidemic_timing_and_model_fit.png",
    "figure_2_national_regret_curves_by_age.pdf",
    "figure_2_national_regret_curves_by_age.png",
    "figure_3_state_timing_and_optimality_probabilities.pdf",
    "figure_3_state_timing_and_optimality_probabilities.png",
    "figure_4_robustness_summary.pdf",
    "figure_4_robustness_summary.png",
    "figure_s1_gam_fit_diagnostics.pdf",
    "figure_s1_gam_fit_diagnostics.png",
    "figure_s2_alternative_outcome_relative_regret.pdf",
    "figure_s2_alternative_outcome_relative_regret.png",
    "figure_s3_national_surveillance_densities.pdf",
    "figure_s3_national_surveillance_densities.png",
    "figure_s4_vaccine_protection_density.pdf",
    "figure_s4_vaccine_protection_density.png"
  )
)
final_table_paths <- file.path(
  "outputs/primary/final_tables",
  c("table_s1_ve_inputs.tex", "table_s1_ve_inputs.csv", "ve_sources.bib",
    "table_s2_model_parameters.tex", "table_s2_model_parameters.csv")
)
figure_component_paths <- file.path(
  "outputs/primary/figure_components",
  c(
    "figure_3a_state_timing_choropleth.pdf",
    "figure_3a_state_timing_choropleth.png"
  )
)
figure_source_paths <- file.path(
  "outputs/primary/figure_source_data",
  c(
    "figure_s1_gam_fit_diagnostics_source_data.csv",
    "figure_s2_alternative_outcome_relative_regret_source_data.csv",
    "figure_s3_national_surveillance_densities_source_data.csv",
    "figure_s4_vaccine_protection_density_source_data.csv",
    "figure_s4_vaccine_protection_density_summary.csv"
  )
)
summary_paths <- c(
  "outputs/sensitivity/summary/tables/robustness_suite_us_summary.csv",
  "outputs/sensitivity/summary/tables/gam_leave_one_season_out_us_summary.csv"
)

required_paths <- c(
  publication_raw_paths,
  publication_processed_paths,
  primary_paths,
  sensitivity_paths,
  loso_table_paths,
  publication_manifest_paths(expected_primary_seasons),
  final_figure_paths,
  final_table_paths,
  figure_component_paths,
  figure_source_paths,
  summary_paths,
  diagnostic_paths
)
required_paths <- unique(required_paths)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing publication output(s):\n", paste(missing_paths, collapse = "\n"))
}
empty_paths <- required_paths[file.info(required_paths)$size <= 0]
if (length(empty_paths) > 0L) {
  stop("Empty publication output(s):\n", paste(empty_paths, collapse = "\n"))
}

# Supplementary input figures must retain complete, normalized source data.
national_density <- read_csv(
  "outputs/primary/figure_source_data/figure_s3_national_surveillance_densities_source_data.csv",
  show_col_types = FALSE
)
season_density <- filter(national_density, .data$curve_type == "season")
density_totals <- season_density %>%
  summarise(mass = sum(.data$normalized_weekly_mass), n_weeks = n(),
            .by = c("measure", "season"))
if (nrow(density_totals) != 36L || any(density_totals$n_weeks != 40L) ||
    any(abs(density_totals$mass - 1) > 1e-12)) {
  stop("National supplementary curves do not contain 36 normalized seasonal profiles.")
}
expected_density_mean <- season_density %>%
  summarise(expected = mean(.data$normalized_weekly_mass),
            .by = c("measure", "week"))
reported_density_mean <- national_density %>%
  filter(.data$curve_type == "equal_season_mean") %>%
  inner_join(expected_density_mean, by = c("measure", "week"),
             relationship = "one-to-one")
if (nrow(reported_density_mean) != 120L ||
    any(abs(reported_density_mean$normalized_weekly_mass - reported_density_mean$expected) > 1e-12)) {
  stop("Representative national curves are not equal-season averages.")
}
if (anyDuplicated(season_density[c("measure", "season", "week")]) ||
    !setequal(season_density$week, c(36L:53L, 1L:22L)) ||
    any(!is.finite(national_density$plot_week))) {
  stop("National supplementary curves do not align to the MMWR week calendar.")
}
source_hashes <- distinct(season_density, .data$source_path, .data$source_md5)
if (any(!file.exists(source_hashes$source_path)) ||
    any(unname(tools::md5sum(source_hashes$source_path)) != source_hashes$source_md5)) {
  stop("National supplementary curves have stale source data.")
}

protection_density <- read_csv(
  "outputs/primary/figure_source_data/figure_s4_vaccine_protection_density_source_data.csv",
  show_col_types = FALSE
)
protection_summary <- read_csv(
  "outputs/primary/figure_source_data/figure_s4_vaccine_protection_density_summary.csv",
  show_col_types = FALSE
)
protection_totals <- protection_density %>%
  summarise(mass = sum(.data$probability), draws = sum(.data$n_draws_in_bin),
            bins = n(), .by = c("analysis", "weeks_since_vaccination"))
if (nrow(protection_totals) != 1644L || any(protection_totals$bins != 200L) ||
    any(protection_totals$draws != expected_draws) ||
    any(abs(protection_totals$mass - 1) > 1e-12) ||
    !setequal(protection_totals$analysis, primary_prefixes) ||
    nrow(protection_summary) != 1644L ||
    !setequal(protection_summary$days_since_vaccination, 0:273) ||
    any(abs(protection_summary$weeks_since_vaccination * 7 -
              protection_summary$days_since_vaccination) > 1e-10) ||
    any(abs(protection_density$protection_bin_upper -
              protection_density$protection_bin_lower - 0.005) > 1e-12) ||
    any(protection_summary$median_protection < 0 | protection_summary$median_protection > 1)) {
  stop("Supplementary protection densities do not match the production draw coverage.")
}
draw_hashes <- distinct(protection_summary, .data$draw_source, .data$draw_source_md5)
if (any(!file.exists(draw_hashes$draw_source)) ||
    any(unname(tools::md5sum(draw_hashes$draw_source)) != draw_hashes$draw_source_md5)) {
  stop("Supplementary protection densities have stale simulation draws.")
}
parameter_table <- read_csv(
  "outputs/primary/final_tables/table_s2_model_parameters.csv", show_col_types = FALSE
)
required_parameter_ids <- c(
  "ve_source_mixture", "immune_lag", "ray_waning", "polynomial_waning",
  "outcome_date_delay", "geography_weights", "ve_calibration", "epidemic_shift"
)
if (anyDuplicated(parameter_table$parameter_id) ||
    !setequal(parameter_table$parameter_id, required_parameter_ids) ||
    anyNA(parameter_table[c("parameter", "specification", "source", "implementation")])) {
  stop("The supplementary empirical-input table is incomplete or lacks source provenance.")
}

input_validation_path <- "outputs/diagnostics/data/publication_input_validation.json"
input_validation <- read_json(input_validation_path, simplifyVector = FALSE)
software <- input_validation$software
required_software_fields <- c(
  "r_version",
  "r_version_numeric",
  "r_platform",
  "os_type",
  "operating_system",
  "locale",
  "time_zone",
  "rng_kind",
  "blas",
  "lapack",
  "cairo_available",
  "sf_external_versions",
  "direct_package_versions",
  "additional_model_package_versions"
)
if (!is.list(software) ||
    length(setdiff(required_software_fields, names(software))) > 0L) {
  stop("The publication input record lacks the required software provenance block.")
}

direct_package_versions <- unlist(
  software$direct_package_versions,
  use.names = TRUE
)
additional_model_package_versions <- unlist(
  software$additional_model_package_versions,
  use.names = TRUE
)
if (!setequal(names(direct_package_versions), publication_required_packages) ||
    length(direct_package_versions) != length(publication_required_packages) ||
    anyNA(direct_package_versions) || any(!nzchar(direct_package_versions))) {
  stop(
    "The publication input record does not contain exactly the direct required ",
    "package versions."
  )
}
if (!setequal(names(additional_model_package_versions), c("Matrix", "nlme")) ||
    length(additional_model_package_versions) != 2L ||
    anyNA(additional_model_package_versions) ||
    any(!nzchar(additional_model_package_versions))) {
  stop("The publication input record lacks Matrix or nlme version provenance.")
}

invalid_version_fields <- names(publication_minimum_package_versions)[vapply(
  names(publication_minimum_package_versions),
  function(package) {
    recorded <- tryCatch(
      numeric_version(direct_package_versions[[package]]),
      error = function(e) NULL
    )
    is.null(recorded) ||
      recorded < numeric_version(publication_minimum_package_versions[[package]])
  },
  logical(1)
)]
recorded_r_version <- tryCatch(
  numeric_version(software$r_version_numeric),
  error = function(e) NULL
)
if (is.null(recorded_r_version) || recorded_r_version < numeric_version("4.2.0")) {
  stop("The publication input record requires R 4.2.0 or later.")
}
if (length(invalid_version_fields) > 0L) {
  stop(
    "The publication input record has package versions below the declared ",
    "minimum: ", paste(invalid_version_fields, collapse = ", "), "."
  )
}

nonempty_scalar <- function(value) {
  length(value) == 1L && !is.na(value) && nzchar(as.character(value))
}
required_os_fields <- c("sysname", "release", "version", "machine")
required_rng_fields <- c("kind", "normal_kind", "sample_kind")
sf_versions <- unlist(software$sf_external_versions, use.names = TRUE)
if (!nonempty_scalar(software$r_version) ||
    !nonempty_scalar(software$r_platform) ||
    !nonempty_scalar(software$os_type) ||
    !setequal(names(software$operating_system), required_os_fields) ||
    any(!vapply(software$operating_system, nonempty_scalar, logical(1))) ||
    !nonempty_scalar(software$locale) ||
    !nonempty_scalar(software$time_zone) ||
    !setequal(names(software$rng_kind), required_rng_fields) ||
    any(!vapply(software$rng_kind, nonempty_scalar, logical(1))) ||
    !nonempty_scalar(software$blas$library) ||
    !nonempty_scalar(software$lapack$library) ||
    !nonempty_scalar(software$lapack$version) ||
    !isTRUE(software$cairo_available) ||
    !all(c("GDAL", "GEOS", "PROJ") %in% names(sf_versions)) ||
    any(!nzchar(sf_versions[c("GDAL", "GEOS", "PROJ")]))) {
  stop("The publication input record has incomplete platform or library provenance.")
}

manifest_paths <- publication_manifest_paths(expected_primary_seasons)
for (path in manifest_paths) {
  manifest <- read_json(path, simplifyVector = TRUE)
  if (manifest$configuration$n_draws != expected_draws) {
    stop(
      path, " records ", manifest$configuration$n_draws,
      " draws; expected ", expected_draws, "."
    )
  }
  if (!is.null(manifest$model_path) &&
      length(manifest$model_path) > 0L &&
      !is.na(manifest$model_path) &&
      nzchar(manifest$model_path) &&
      !file.exists(manifest$model_path)) {
    stop(path, " references a missing fitted model: ", manifest$model_path)
  }
  if (!is.data.frame(manifest$inputs) ||
      !all(c("path", "md5") %in% names(manifest$inputs))) {
    stop(path, " does not contain a tabular input-hash record.")
  }
  missing_manifest_inputs <- manifest$inputs$path[!file.exists(manifest$inputs$path)]
  if (length(missing_manifest_inputs) > 0L) {
    stop(
      path, " records missing input(s): ",
      paste(missing_manifest_inputs, collapse = ", ")
    )
  }
  current_md5 <- unname(tools::md5sum(manifest$inputs$path))
  stale_inputs <- manifest$inputs$path[current_md5 != manifest$inputs$md5]
  if (length(stale_inputs) > 0L) {
    stop(
      path, " was generated from an older version of: ",
      paste(stale_inputs, collapse = ", ")
    )
  }
}

manifest_path_for_prefix <- function(prefix) {
  if (prefix %in% primary_prefixes) {
    file.path(
      "outputs/primary/tables",
      paste0(prefix, "_analysis_manifest.json")
    )
  } else {
    file.path(
      "outputs/sensitivity",
      prefix,
      "tables",
      paste0(prefix, "_analysis_manifest.json")
    )
  }
}

configuration_value_matches <- function(observed, expected) {
  if (is.numeric(expected)) {
    return(identical(as.numeric(observed), as.numeric(expected)))
  }
  if (is.logical(expected)) {
    return(identical(as.logical(observed), as.logical(expected)))
  }
  identical(as.character(observed), as.character(expected))
}

assert_manifest_configuration <- function(
    prefix,
    expected,
    model_path_pattern = NULL) {
  path <- manifest_path_for_prefix(prefix)
  manifest <- read_json(path, simplifyVector = TRUE)
  configuration <- manifest$configuration
  missing_fields <- setdiff(names(expected), names(configuration))
  mismatched_fields <- names(expected)[vapply(
    names(expected),
    function(field) {
      field %in% names(configuration) &&
        !configuration_value_matches(configuration[[field]], expected[[field]])
    },
    logical(1)
  )]
  if (!identical(as.character(manifest$analysis), prefix) ||
      length(missing_fields) > 0L || length(mismatched_fields) > 0L) {
    stop(
      path, " does not match the prespecified configuration. Missing fields: ",
      paste(missing_fields, collapse = ", "), "; mismatched fields: ",
      paste(mismatched_fields, collapse = ", "), "."
    )
  }
  if (!is.null(model_path_pattern) &&
      (is.null(manifest$model_path) ||
        length(manifest$model_path) != 1L ||
        is.na(manifest$model_path) ||
        !grepl(model_path_pattern, basename(manifest$model_path), fixed = TRUE))) {
    stop(
      path, " does not reference the prespecified fitted-model configuration ",
      model_path_pattern, "."
    )
  }
  invisible(manifest)
}

baseline_configuration <- list(
  n_draws = expected_draws,
  primary_season_start_years = expected_primary_seasons,
  ve_reference_weeks = 8,
  waning_prior = "ray_2019",
  immune_lag_prior = "discrete_about_2w",
  timing_shift_sd_weeks = 0.75
)
baseline_gam_configuration <- utils::modifyList(
  baseline_configuration,
  list(
    outcome_date_prior = "none",
    gam_global_k = 14L,
    gam_season_k = 12L,
    gam_gamma = 1
  )
)

assert_manifest_configuration("gam_primary", baseline_gam_configuration)

age_configuration <- list(
  gam_age_0_4 = c(age_group = "0-4", ve_age_group = "6m-8"),
  gam_age_5_24 = c(age_group = "5-24", ve_age_group = "9-17"),
  gam_age_25_49 = c(age_group = "25-49", ve_age_group = "18-49"),
  gam_age_50_64 = c(age_group = "50-64", ve_age_group = "50-64"),
  gam_age_65_ = c(age_group = "65+", ve_age_group = "65+")
)
invisible(lapply(names(age_configuration), function(prefix) {
  assert_manifest_configuration(
    prefix,
    utils::modifyList(
      baseline_configuration,
      c(
        as.list(age_configuration[[prefix]]),
        list(gam_global_k = 14L, gam_season_k = 12L, gam_gamma = 1)
      )
    ),
    model_path_pattern = "_gk14_sk12_gamma1_"
  )
}))

assert_manifest_configuration(
  "empirical_primary",
  utils::modifyList(
    baseline_configuration,
    list(
      analysis_scope = "national all-age ILINet",
      empirical_curve_smoothing = "none",
      include_observation_noise = TRUE,
      season_sampling = "with replacement among eligible seasons extending through MMWR week 22",
      timing_shift_interpolation = "linear",
      timing_shift_boundary_rule = "hold the nearest observed boundary value constant"
    )
  )
)

gam_sensitivity_overrides <- list(
  gam_no_timing_shift_primary = list(timing_shift_sd_weeks = 0),
  gam_precovid_primary = list(primary_season_start_years = 2010:2018),
  gam_fixed_immune_lag_primary = list(immune_lag_prior = "fixed_2w"),
  gam_slow_waning_primary = list(waning_prior = "ray_2019_lower95"),
  gam_fast_waning_primary = list(waning_prior = "ray_2019_upper95"),
  gam_spencer_ferdinands_fast_waning = list(
    waning_prior = "spencer_ferdinands_fast"
  ),
  gam_spencer_ferdinands_onset_aligned = list(
    waning_prior = "spencer_ferdinands_fast",
    outcome_date_prior = "flu_ve_outpatient_onset_to_enrollment"
  ),
  gam_spencer_slow_waning = list(waning_prior = "spencer_slow"),
  gam_smooth_k16_k14_primary = list(
    gam_global_k = 16L,
    gam_season_k = 14L
  ),
  gam_ve_reference_4w = list(ve_reference_weeks = 4),
  gam_ve_reference_12w = list(ve_reference_weeks = 12)
)
invisible(lapply(names(gam_sensitivity_overrides), function(prefix) {
  assert_manifest_configuration(
    prefix,
    utils::modifyList(
      baseline_gam_configuration,
      gam_sensitivity_overrides[[prefix]]
    )
  )
}))

# The no-shift sensitivity must differ only in the timing perturbation.
primary_manifest <- read_json(
  manifest_path_for_prefix("gam_primary"), simplifyVector = TRUE
)
no_shift_manifest <- read_json(
  manifest_path_for_prefix("gam_no_timing_shift_primary"), simplifyVector = TRUE
)
no_shift_expected_configuration <- utils::modifyList(
  primary_manifest$configuration,
  list(timing_shift_sd_weeks = 0)
)
if (!isTRUE(all.equal(
      no_shift_manifest$configuration, no_shift_expected_configuration
    )) ||
    !identical(no_shift_manifest$model_path, primary_manifest$model_path) ||
    !identical(
      no_shift_manifest$inputs[c("path", "md5")],
      primary_manifest$inputs[c("path", "md5")]
    )) {
  stop("No-shift sensitivity differs from the primary analysis beyond timing shift.")
}
paired_ve_paths <- c(
  "outputs/primary/tables/gam_primary_ve_waning_draws.csv",
  paste0(
    "outputs/sensitivity/gam_no_timing_shift_primary/tables/",
    "gam_no_timing_shift_primary_ve_waning_draws.csv"
  )
)
if (length(unique(unname(tools::md5sum(paired_ve_paths)))) != 1L) {
  stop("No-shift sensitivity does not use the primary vaccine-protection draws.")
}

nrevss_configuration <- list(
  gam_nrevss_national = list(
    primary_season_start_years = expected_primary_seasons,
    analysis_region = "national",
    target_virus = "all"
  ),
  gam_nrevss_post2015_flu_a = list(
    primary_season_start_years = c(2015:2018, 2023:2025),
    analysis_region = "national",
    target_virus = "a"
  )
)
invisible(lapply(names(nrevss_configuration), function(prefix) {
  assert_manifest_configuration(
    prefix,
    utils::modifyList(
      baseline_configuration,
      nrevss_configuration[[prefix]]
    ),
    model_path_pattern = "_gk14_sk12_gamma1_"
  )
}))

ili_plus_configuration <- list(
  empirical_ili_plus_all_seasons_weighted = list(
    primary_season_start_years = expected_primary_seasons,
    available_season_start_years = expected_primary_seasons,
    ili_weighting = "cdc_population_weighted",
    nrevss_lab_scope = "all_available"
  ),
  empirical_ili_plus_all_seasons_unweighted = list(
    primary_season_start_years = expected_primary_seasons,
    available_season_start_years = expected_primary_seasons,
    ili_weighting = "national_unweighted",
    nrevss_lab_scope = "all_available"
  ),
  empirical_ili_plus_clinical_only_weighted = list(
    primary_season_start_years = c(2016:2018, 2023:2025),
    available_season_start_years = c(2016:2018, 2023:2025),
    ili_weighting = "cdc_population_weighted",
    nrevss_lab_scope = "clinical_only"
  )
)
invisible(lapply(names(ili_plus_configuration), function(prefix) {
  assert_manifest_configuration(
    prefix,
    utils::modifyList(
      baseline_configuration,
      utils::modifyList(
        list(smoothing_spar = 0.65),
        ili_plus_configuration[[prefix]]
      )
    )
  )
}))

invisible(lapply(expected_primary_seasons, function(omitted_year) {
  prefix <- paste0("gam_loso_drop_", omitted_year)
  assert_manifest_configuration(
    prefix,
    utils::modifyList(
      baseline_gam_configuration,
      list(
        primary_season_start_years = setdiff(
          expected_primary_seasons,
          omitted_year
        )
      )
    )
  )
}))

state_input <- read_csv(
  "data/processed/ilinet_state_all_age_primary_seasons.csv",
  show_col_types = FALSE
)
ve_input <- read_csv(
  "data/processed/cdc_ve_estimates.csv", show_col_types = FALSE
)
if (length(unique(unname(tools::md5sum(c(
  "data/processed/cdc_ve_estimates.csv",
  "outputs/primary/final_tables/table_s1_ve_inputs.csv"
))))) != 1L) {
  stop("The supplementary VE source table does not match the analysis input.")
}
ve_bibliography <- paste(readLines(
  "outputs/primary/final_tables/ve_sources.bib", warn = FALSE
), collapse = "\n")
missing_ve_citations <- unique(ve_input$citation_key)[!vapply(
  unique(ve_input$citation_key),
  function(key) grepl(paste0("{", key, ","), ve_bibliography, fixed = TRUE),
  logical(1)
)]
if (length(missing_ve_citations) > 0L) {
  stop("The supplementary VE bibliography lacks: ",
    paste(missing_ve_citations, collapse = ", "))
}
expected_state_names <- sort(unique(c(state_input$state, "US")))
age_group_by_prefix <- c(
  gam_age_0_4 = "0-4",
  gam_age_5_24 = "5-24",
  gam_age_25_49 = "25-49",
  gam_age_50_64 = "50-64",
  gam_age_65_ = "65+"
)

expected_groups_for_prefix <- function(prefix) {
  if (prefix %in% publication_us_only_prefixes()) {
    return(tibble::tibble(state = "US", age_group = "all"))
  }
  age_group <- if (prefix %in% names(age_group_by_prefix)) {
    unname(age_group_by_prefix[[prefix]])
  } else {
    "all"
  }
  tibble::tibble(
    state = expected_state_names,
    age_group = age_group
  )
}

assert_expected_groups <- function(data, prefix, table_label) {
  expected <- expected_groups_for_prefix(prefix)
  observed <- data %>%
    distinct(.data$state, .data$age_group)
  missing <- anti_join(
    expected,
    observed,
    by = c("state", "age_group")
  )
  extra <- anti_join(
    observed,
    expected,
    by = c("state", "age_group")
  )
  if (nrow(missing) > 0L || nrow(extra) > 0L) {
    stop(
      prefix, " has incorrect state-age coverage in ", table_label,
      ". Missing groups: ",
      paste(paste(missing$state, missing$age_group, sep = "/"), collapse = ", "),
      ". Unexpected groups: ",
      paste(paste(extra$state, extra$age_group, sep = "/"), collapse = ", "),
      "."
    )
  }
}

check_decision_tables <- function(prefix, table_dir) {
  draws <- read_csv(
    file.path(table_dir, paste0(prefix, "_optimal_week_draws.csv")),
    show_col_types = FALSE
  )
  cutoffs <- read_csv(
    file.path(table_dir, paste0(prefix, "_week_cutoff_probabilities.csv")),
    show_col_types = FALSE
  )
  regrets <- read_csv(
    file.path(table_dir, paste0(prefix, "_regret_curve.csv")),
    show_col_types = FALSE
  )
  intervals <- read_csv(
    file.path(table_dir, paste0(prefix, "_optimal_week_intervals.csv")),
    show_col_types = FALSE
  )
  distribution <- read_csv(
    file.path(table_dir, paste0(prefix, "_optimal_week_distribution.csv")),
    show_col_types = FALSE
  )
  near_optimal <- read_csv(
    file.path(table_dir, paste0(prefix, "_near_optimal_curve.csv")),
    show_col_types = FALSE
  )
  ve_draws <- read_csv(
    file.path(table_dir, paste0(prefix, "_ve_waning_draws.csv")),
    show_col_types = FALSE
  )

  assert_expected_groups(draws, prefix, "optimal-week draws")
  assert_expected_groups(cutoffs, prefix, "cutoff probabilities")
  assert_expected_groups(regrets, prefix, "regret curve")
  assert_expected_groups(intervals, prefix, "optimal-week intervals")
  assert_expected_groups(distribution, prefix, "optimal-week distribution")
  assert_expected_groups(near_optimal, prefix, "near-optimal curve")

  draw_checks <- draws %>%
    summarise(
      n_rows = n(),
      n_draws = n_distinct(.data$draw),
      valid_draw_ids = setequal(.data$draw, seq_len(expected_draws)),
      valid_weeks = all(.data$optimal_week %in% c(36:52, 1:12)),
      .by = c("state", "age_group")
    )
  if (any(
    draw_checks$n_rows != expected_draws |
      draw_checks$n_draws != expected_draws |
      !draw_checks$valid_draw_ids |
      !draw_checks$valid_weeks
  )) {
    stop(prefix, " has invalid or incomplete optimal-week draws.")
  }

  cutoff_checks <- cutoffs %>%
    arrange(.data$state, .data$age_group, .data$cutoff_week) %>%
    summarise(
      valid_cutoffs = setequal(
        .data$cutoff_week,
        c(39L, 44L, 48L, 52L)
      ),
      valid_probabilities = all(
        is.finite(.data$probability) &
          .data$probability >= 0 &
          .data$probability <= 1
      ),
      monotone = all(diff(.data$probability) >= -1e-12),
      .by = c("state", "age_group")
    )
  if (any(
    !cutoff_checks$valid_cutoffs |
      !cutoff_checks$valid_probabilities |
      !cutoff_checks$monotone
  ) ||
      any(!is.finite(cutoffs$probability)) ||
      any(cutoffs$probability < 0 | cutoffs$probability > 1)) {
    stop(prefix, " has invalid cutoff probabilities.")
  }

  regret_checks <- regrets %>%
    summarise(
      valid_weeks = setequal(
        .data$vaccination_week,
        c(36:52, 1:12)
      ),
      .by = c("state", "age_group")
    )
  if (any(!regret_checks$valid_weeks)) {
    stop(prefix, " has incomplete vaccination-week regret curves.")
  }
  numeric_regret <- regrets %>%
    select(
      "mean_regret", "lower_50", "upper_50", "lower_95", "upper_95"
    )
  if (any(!is.finite(as.matrix(numeric_regret))) ||
      any(as.matrix(numeric_regret) < -1e-12)) {
    stop(prefix, " has invalid regret summaries.")
  }

  interval_counts <- intervals %>%
    count(.data$state, .data$age_group, name = "n_rows")
  if (any(interval_counts$n_rows != 1L)) {
    stop(prefix, " must have one optimal-week interval row per group.")
  }

  distribution_checks <- distribution %>%
    summarise(
      complete_week_grid = setequal(.data$optimal_week, 1:52),
      zero_outside_candidate_window = all(
        .data$n[!.data$optimal_week %in% c(36:52, 1:12)] == 0L &
          .data$probability[!.data$optimal_week %in% c(36:52, 1:12)] == 0
      ),
      valid_draw_counts = all(
        .data$n >= 0L & .data$n_draws == expected_draws
      ) && sum(.data$n) == expected_draws,
      valid_probabilities = all(
        is.finite(.data$probability) &
          .data$probability >= 0 &
          .data$probability <= 1 &
          abs(.data$probability - .data$n / expected_draws) < 1e-12
      ) && abs(sum(.data$probability) - 1) < 1e-12,
      .by = c("state", "age_group")
    )
  if (any(
    !distribution_checks$complete_week_grid |
      !distribution_checks$zero_outside_candidate_window |
      !distribution_checks$valid_draw_counts |
      !distribution_checks$valid_probabilities
  )) {
    stop(prefix, " has an invalid optimal-week distribution.")
  }

  near_optimal_checks <- near_optimal %>%
    summarise(
      valid_weeks = setequal(.data$vaccination_week, c(36:52, 1:12)),
      valid_tolerances = setequal(
        .data$relative_regret_tolerance,
        c(0.01, 0.05, 0.10)
      ),
      complete_grid = n() == length(c(36:52, 1:12)) * 3L &&
        !anyDuplicated(paste(
          .data$vaccination_week,
          .data$relative_regret_tolerance
        )),
      valid_probabilities = all(
        is.finite(.data$probability_near_optimal) &
          .data$probability_near_optimal >= 0 &
          .data$probability_near_optimal <= 1
      ),
      valid_regret = all(
        is.finite(.data$mean_relative_regret) &
          .data$mean_relative_regret >= -1e-12
      ),
      .by = c("state", "age_group")
    )
  if (any(
    !near_optimal_checks$valid_weeks |
      !near_optimal_checks$valid_tolerances |
      !near_optimal_checks$complete_grid |
      !near_optimal_checks$valid_probabilities |
      !near_optimal_checks$valid_regret
  )) {
    stop(prefix, " has an invalid or incomplete near-optimal curve.")
  }

  expected_analysis_age_group <- unique(expected_groups_for_prefix(prefix)$age_group)
  expected_source_age_group <- if (prefix %in% names(age_configuration)) {
    unname(age_configuration[[prefix]][["ve_age_group"]])
  } else {
    "all"
  }
  expected_ve_pool <- ve_input %>%
    filter(.data$ve_age_group == expected_source_age_group)
  invalid_ve_sources <- ve_draws %>%
    anti_join(
      expected_ve_pool,
      by = c("ve_source_season" = "season", "ve_source" = "source_id")
    )
  if (nrow(invalid_ve_sources) > 0L ||
      !setequal(ve_draws$ve_source_season, expected_ve_pool$season)) {
    stop(prefix, " has vaccine-protection draws outside its specified age-season pool.")
  }
  valid_waning_parameters <- (
    ve_draws$protection_model == "exponential_effect_ratio" &
      is.finite(ve_draws$beta_wane_per_28d)
  ) | (
    ve_draws$protection_model %in% c(
      "spencer_fast_relative_ve",
      "spencer_slow_relative_ve"
    ) &
      is.finite(ve_draws$reference_retention)
  )
  if (nrow(ve_draws) != expected_draws ||
      !setequal(ve_draws$draw, seq_len(expected_draws)) ||
      !identical(sort(unique(ve_draws$age_group)), expected_analysis_age_group) ||
      any(!is.finite(ve_draws$initial_ve)) ||
      any(ve_draws$initial_ve < 0 | ve_draws$initial_ve > 1) ||
      any(!valid_waning_parameters) ||
      any(!is.finite(ve_draws$immune_lag_weeks))) {
    stop(prefix, " has invalid or incomplete vaccine-protection draws.")
  }

  # Replay only the reported-VE draws to verify age selection and CI sampling.
  # The all-age and five age-specific runners have fixed, distinct VE seeds.
  if (prefix %in% primary_prefixes) {
    ve_seed <- if (prefix == "gam_primary") {
      20260509L
    } else {
      20260520L + match(prefix, names(age_configuration))
    }
    expected_reported_draws <- draw_reported_ve_from_estimates(
      ve_input,
      n_draws = expected_draws,
      age_group = expected_analysis_age_group,
      ve_age_group = expected_source_age_group,
      seed = ve_seed
    )
    observed_reported_draws <- ve_draws %>%
      arrange(.data$draw) %>%
      select(all_of(names(expected_reported_draws)))
    if (!isTRUE(all.equal(
      observed_reported_draws, expected_reported_draws,
      tolerance = 1e-12, check.attributes = FALSE
    ))) {
      stop(prefix, " does not reproduce the specified age-specific VE sampling.")
    }
  }
}

invisible(lapply(
  primary_prefixes,
  check_decision_tables,
  table_dir = "outputs/primary/tables"
))
invisible(lapply(
  sensitivity_prefixes,
  function(prefix) check_decision_tables(
    prefix,
    file.path("outputs/sensitivity", prefix, "tables")
  )
))
invisible(lapply(
  loso_prefixes,
  function(prefix) check_decision_tables(
    prefix,
    file.path("outputs/sensitivity", prefix, "tables")
  )
))

held_out_path <- paste0(
  "outputs/diagnostics/held_out_validation/",
  "held_out_season_decision_validation.csv"
)
held_out_summary_path <- paste0(
  "outputs/diagnostics/held_out_validation/",
  "held_out_season_decision_validation_summary.csv"
)
held_out <- read_csv(held_out_path, show_col_types = FALSE)
held_out_summary <- read_csv(held_out_summary_path, show_col_types = FALSE)
held_out_probability_columns <- c(
  "held_out_optimal_in_training_50_interval",
  "held_out_optimal_in_training_95_interval",
  "p_by_week_39",
  "p_by_week_44",
  "p_by_week_48",
  "p_by_week_52"
)
held_out_numeric_columns <- c(
  "held_out_regret_at_training_week",
  "held_out_visits_lost_per_10000_at_training_week",
  "held_out_relative_regret_at_training_week",
  "optimal_week_log_score",
  "optimal_week_brier_score"
)
if (nrow(held_out) != length(expected_primary_seasons) ||
    anyDuplicated(held_out$held_out_season_start_year) ||
    !setequal(held_out$held_out_season_start_year, expected_primary_seasons) ||
    any(!held_out$training_regret_minimizing_week %in% c(36:52, 1:12)) ||
    any(!held_out$held_out_regret_minimizing_week %in% c(36:52, 1:12)) ||
    any(!held_out$held_out_median_optimal_week %in% c(36:52, 1:12)) ||
    any(!is.finite(as.matrix(held_out[held_out_probability_columns]))) ||
    any(as.matrix(held_out[held_out_probability_columns]) < 0) ||
    any(as.matrix(held_out[held_out_probability_columns]) > 1) ||
    any(!is.finite(as.matrix(held_out[held_out_numeric_columns]))) ||
    any(as.matrix(held_out[held_out_numeric_columns]) < 0) ||
    any(apply(
      held_out[c("p_by_week_39", "p_by_week_44", "p_by_week_48", "p_by_week_52")],
      1,
      function(probability) any(diff(probability) < -1e-12)
    ))) {
  stop("Held-out-season decision validation is incomplete or invalid.")
}

expected_held_out_summary <- tibble::tibble(
  n_held_out_seasons = nrow(held_out),
  median_visits_lost_per_10000 = stats::median(
    held_out$held_out_visits_lost_per_10000_at_training_week
  ),
  maximum_visits_lost_per_10000 = max(
    held_out$held_out_visits_lost_per_10000_at_training_week
  ),
  median_relative_regret = stats::median(
    held_out$held_out_relative_regret_at_training_week
  ),
  mean_50_interval_coverage = mean(
    held_out$held_out_optimal_in_training_50_interval
  ),
  mean_95_interval_coverage = mean(
    held_out$held_out_optimal_in_training_95_interval
  ),
  mean_optimal_week_log_score = mean(held_out$optimal_week_log_score),
  mean_optimal_week_brier_score = mean(held_out$optimal_week_brier_score)
)
if (nrow(held_out_summary) != 1L ||
    !isTRUE(all.equal(
      held_out_summary,
      expected_held_out_summary,
      tolerance = 1e-12,
      check.attributes = FALSE
    ))) {
  stop("Held-out-season summary does not reproduce the season-level records.")
}

ili_plus <- read_csv(
  "data/processed/national_ili_plus_primary_seasons.csv",
  show_col_types = FALSE
)
expected_season_start_years <- expected_primary_seasons
if (!setequal(
      unique(ili_plus$season_start_year),
      expected_season_start_years
    ) ||
    any(!is.finite(ili_plus$ili_plus_weighted)) ||
    any(!is.finite(ili_plus$ili_plus_unweighted))) {
  stop(
    "National ILI+ input failed season-set or finite-value validation. ",
    "Expected seasons: ",
    paste(expected_season_start_years, collapse = ", "),
    "."
  )
}

rebuilt_ili_plus <- prepare_national_ili_plus(
  read_csv(
    "data/processed/ilinet_national_all_age_primary_seasons.csv",
    show_col_types = FALSE
  ) %>%
    filter_to_primary_seasons("national ILINet validation data"),
  read_csv(
    "data/processed/nrevss_clinical_labs_primary_seasons.csv",
    show_col_types = FALSE
  ) %>%
    filter_to_primary_seasons("national NREVSS validation data")
)
comparison_columns <- c(
  "season", "season_start_year", "mmwr_year", "mmwr_week", "week_start",
  "ilinet_ili_visits", "ilinet_total_visits",
  "ilinet_unweighted_ili_prob", "ilinet_weighted_ili_prob",
  "nrevss_positive_tests", "nrevss_total_tests", "nrevss_positive_prob",
  "nrevss_lab_series", "ili_plus_weighted", "ili_plus_unweighted"
)
stored_comparison <- ili_plus %>%
  select(all_of(comparison_columns)) %>%
  arrange(.data$season_start_year, .data$mmwr_year, .data$mmwr_week)
rebuilt_comparison <- rebuilt_ili_plus %>%
  select(all_of(comparison_columns)) %>%
  arrange(.data$season_start_year, .data$mmwr_year, .data$mmwr_week)
if (!isTRUE(all.equal(
  stored_comparison,
  rebuilt_comparison,
  tolerance = 1e-12,
  check.attributes = FALSE
))) {
  stop(
    "Stored national ILI+ data do not match a fresh join of the current ",
    "ILINet and NREVSS inputs."
  )
}

clinical_ili_plus <- filter_ili_plus_lab_scope(
  ili_plus,
  "clinical_only"
)
clinical_season_start_years <- sort(unique(
  clinical_ili_plus$season_start_year
))
clinical_coverage <- clinical_ili_plus %>%
  summarise(
    all_weeks_clinical = all(
      .data$nrevss_lab_series == "clinical_labs"
    ),
    has_window_start = any(.data$mmwr_week == 36L),
    has_window_end = any(.data$mmwr_week == 22L),
    .by = c("season", "season_start_year")
  )
if (any(
  !clinical_coverage$all_weeks_clinical |
    !clinical_coverage$has_window_start |
    !clinical_coverage$has_window_end
)) {
  stop(
    "Clinical-laboratory-only ILI+ contains a mixed-series or incomplete season."
  )
}

if (validate_artifact_hygiene) {
  assert_no_unexpected_files(
    "data/raw",
    c("data/raw/.gitkeep", publication_raw_paths),
    "raw-data"
  )
  assert_no_unexpected_files(
    "data/processed",
    c("data/processed/.gitkeep", publication_processed_paths),
    "processed-data"
  )
  assert_no_unexpected_files(
    "outputs/diagnostics",
    diagnostic_paths,
    "diagnostic-output"
  )
}

clinical_manifest_path <- paste0(
  "outputs/sensitivity/empirical_ili_plus_clinical_only_weighted/tables/",
  "empirical_ili_plus_clinical_only_weighted_analysis_manifest.json"
)
clinical_manifest <- read_json(
  clinical_manifest_path,
  simplifyVector = TRUE
)
if (!setequal(
      as.integer(
        clinical_manifest$configuration$available_season_start_years
      ),
      clinical_season_start_years
    ) ||
    !setequal(
      as.integer(
        clinical_manifest$configuration$primary_season_start_years
      ),
      clinical_season_start_years
    )) {
  stop(
    "Clinical-laboratory-only ILI+ manifest does not match the complete, ",
    "single-series season set: ",
    paste(clinical_season_start_years, collapse = ", "),
    "."
  )
}

aligned_manifest_path <- paste0(
  "outputs/sensitivity/gam_spencer_ferdinands_onset_aligned/tables/",
  "gam_spencer_ferdinands_onset_aligned_analysis_manifest.json"
)
unaligned_manifest_path <- paste0(
  "outputs/sensitivity/gam_spencer_ferdinands_fast_waning/tables/",
  "gam_spencer_ferdinands_fast_waning_analysis_manifest.json"
)
aligned_manifest <- read_json(
  aligned_manifest_path,
  simplifyVector = TRUE
)
unaligned_manifest <- read_json(
  unaligned_manifest_path,
  simplifyVector = TRUE
)
if (aligned_manifest$configuration$waning_prior !=
      "spencer_ferdinands_fast" ||
    aligned_manifest$configuration$outcome_date_prior !=
      "flu_ve_outpatient_onset_to_enrollment" ||
    unaligned_manifest$configuration$waning_prior !=
      "spencer_ferdinands_fast" ||
    unaligned_manifest$configuration$outcome_date_prior != "none" ||
    aligned_manifest$model_path != unaligned_manifest$model_path) {
  stop("Onset-aligned sensitivity manifest has inconsistent time origins.")
}

paired_configuration_fields <- c(
  "n_draws",
  "primary_season_start_years",
  "burden_scale",
  "ve_reference_weeks",
  "waning_prior",
  "immune_lag_prior",
  "timing_shift_sd_weeks",
  "gam_global_k",
  "gam_season_k",
  "gam_gamma"
)
if (!isTRUE(all.equal(
  aligned_manifest$configuration[paired_configuration_fields],
  unaligned_manifest$configuration[paired_configuration_fields],
  check.attributes = FALSE
))) {
  stop(
    "The aligned and unaligned Spencer-Ferdinands analyses differ in ",
    "configuration beyond outcome-date alignment."
  )
}
if (aligned_manifest$configuration$outcome_date_category_counts !=
      "0-2 days:8520; 3-4 days:10754; 5-7 days:7861" ||
    aligned_manifest$configuration$outcome_date_seed != 20260531L) {
  stop("Onset-aligned manifest is missing the specified delay distribution.")
}

aligned_manifest_raw <- read_json(
  aligned_manifest_path,
  simplifyVector = FALSE
)
outcome_code_input <- Filter(
  function(input) identical(input$path, "R/outcome_time.R"),
  aligned_manifest_raw$inputs
)
if (length(outcome_code_input) != 1L ||
    !identical(
      outcome_code_input[[1]]$md5,
      unname(tools::md5sum("R/outcome_time.R"))
    )) {
  stop(
    "Onset-aligned manifest does not fingerprint the current ",
    "outcome-date implementation."
  )
}

if (validate_artifact_hygiene) {
  allowed_figure_roots <- c(
    normalizePath("outputs/primary/final_figures", mustWork = TRUE),
    normalizePath("outputs/primary/figure_components", mustWork = TRUE)
  )
  all_figure_paths <- list.files(
    c("outputs/primary", "outputs/sensitivity"),
    pattern = "[.](png|pdf)$",
    recursive = TRUE,
    full.names = TRUE
  )
  unexpected_figures <- all_figure_paths[
    !vapply(
      normalizePath(all_figure_paths, mustWork = TRUE),
      function(path) any(startsWith(path, paste0(allowed_figure_roots, "/"))),
      logical(1)
    )
  ]
  if (length(unexpected_figures) > 0L) {
    stop(
      "Unexpected exploratory figure artifact(s):\n",
      paste(unexpected_figures, collapse = "\n")
    )
  }

  state_gam_sensitivity_prefixes <- setdiff(
    sensitivity_prefixes,
    c(
      "empirical_primary",
      publication_us_only_prefixes()
    )
  )
  allowed_primary_files <- c(
    primary_paths,
    final_figure_paths,
    final_table_paths,
    figure_component_paths,
    figure_source_paths,
    file.path(
      "outputs/primary/tables",
      c(
        "gam_primary_ilinet_statewise_diagnostics.csv",
        "gam_age_ilinet_statewise_diagnostics.csv",
        "gam_age_specific_optimal_week_intervals.csv"
      )
    )
  )
  allowed_sensitivity_files <- c(
    sensitivity_paths,
    loso_table_paths,
    file.path(
      "outputs/sensitivity",
      c(state_gam_sensitivity_prefixes, loso_prefixes),
      "tables",
      paste0(
        c(state_gam_sensitivity_prefixes, loso_prefixes),
        "_ilinet_statewise_diagnostics.csv"
      )
    ),
    file.path(
      "outputs/sensitivity",
      c("gam_nrevss_national", "gam_nrevss_post2015_flu_a"),
      "tables",
      paste0(
        c("gam_nrevss_national", "gam_nrevss_post2015_flu_a"),
        "_diagnostics.csv"
      )
    ),
    summary_paths[startsWith(
      summary_paths,
      "outputs/sensitivity/"
    )]
  )
  actual_publication_files <- list.files(
    c("outputs/primary", "outputs/sensitivity"),
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = FALSE
  )
  unexpected_publication_files <- setdiff(
    actual_publication_files,
    unique(c(allowed_primary_files, allowed_sensitivity_files))
  )
  if (length(unexpected_publication_files) > 0L) {
    stop(
      "Unexpected primary or sensitivity output artifact(s):\n",
      paste(unexpected_publication_files, collapse = "\n")
    )
  }

  temporary_artifacts <- list.files(
    "outputs",
    pattern = "([.]DS_Store|Rplots[.]pdf|[.]tmp|[.]log)$",
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE
  )
  if (length(temporary_artifacts) > 0L) {
    stop(
      "Temporary output artifact(s) remain:\n",
      paste(temporary_artifacts, collapse = "\n")
    )
  }
}

message(
  "Publication output validation passed for ",
  length(primary_prefixes), " primary and ",
  length(sensitivity_prefixes), " manuscript sensitivity analyses, plus ",
  length(loso_prefixes), " leave-one-season-out analyses."
)
