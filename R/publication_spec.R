# Canonical publication artifact registry shared by validation and cache pruning.
# Keep these prefixes synchronized with run_manuscript_sensitivity_suite.sh.

publication_ve_pool_coverage <- function() {
  age_groups <- c("6m-8", "9-17", "18-49", "50-64", "65+")
  pool_seasons <- c(
    list(all = c(2010:2018, 2023:2025)),
    stats::setNames(rep(list(c(2010:2018, 2023L)), length(age_groups)), age_groups)
  )
  data.frame(
    ve_age_group = rep(names(pool_seasons), lengths(pool_seasons)),
    season_start_year = unlist(pool_seasons, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

publication_primary_prefixes <- function() {
  c(
    "gam_primary",
    "gam_age_0_4",
    "gam_age_5_24",
    "gam_age_25_49",
    "gam_age_50_64",
    "gam_age_65_"
  )
}

publication_sensitivity_prefixes <- function() {
  c(
    "empirical_primary",
    "gam_no_timing_shift_primary",
    "gam_precovid_primary",
    "gam_fixed_immune_lag_primary",
    "gam_slow_waning_primary",
    "gam_fast_waning_primary",
    "gam_spencer_ferdinands_fast_waning",
    "gam_spencer_ferdinands_onset_aligned",
    "gam_spencer_slow_waning",
    "gam_smooth_k16_k14_primary",
    "gam_ve_reference_4w",
    "gam_ve_reference_12w",
    "gam_nrevss_national",
    "gam_nrevss_post2015_flu_a",
    "empirical_ili_plus_all_seasons_weighted",
    "empirical_ili_plus_all_seasons_unweighted",
    "empirical_ili_plus_clinical_only_weighted"
  )
}

publication_loso_prefixes <- function(
    season_start_years = primary_season_start_years()) {
  paste0("gam_loso_drop_", sort(unique(as.integer(season_start_years))))
}

publication_us_only_prefixes <- function() {
  c(
    "empirical_primary",
    "gam_nrevss_national",
    "gam_nrevss_post2015_flu_a",
    "empirical_ili_plus_all_seasons_weighted",
    "empirical_ili_plus_all_seasons_unweighted",
    "empirical_ili_plus_clinical_only_weighted"
  )
}

publication_manifest_paths <- function(
    season_start_years = primary_season_start_years()) {
  primary <- file.path(
    "outputs/primary/tables",
    paste0(publication_primary_prefixes(), "_analysis_manifest.json")
  )
  sensitivity_prefixes <- c(
    publication_sensitivity_prefixes(),
    publication_loso_prefixes(season_start_years)
  )
  sensitivity <- file.path(
    "outputs/sensitivity",
    sensitivity_prefixes,
    "tables",
    paste0(sensitivity_prefixes, "_analysis_manifest.json")
  )
  c(primary, sensitivity)
}
