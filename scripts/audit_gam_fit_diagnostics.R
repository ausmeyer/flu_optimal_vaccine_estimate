#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
})

diagnostic_files <- list.files(
  c("outputs/primary", "outputs/sensitivity"),
  pattern = "^gam.*diagnostics\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)
diagnostic_files <- diagnostic_files[!grepl("gam_smoothing_tuning_", diagnostic_files)]

if (length(diagnostic_files) == 0L) {
  stop("No GAM diagnostic CSV files found under outputs/primary or outputs/sensitivity.")
}

diagnostics <- purrr::map_dfr(diagnostic_files, function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      diagnostic_file = basename(path),
      analysis = sub("_diagnostics\\.csv$", "", basename(path))
    )
})

diagnostics <- diagnostics %>%
  mutate(
    global_basis_warning = .data$global_k_index < 0.90 &
      .data$global_k_check_p_value < 0.05,
    season_basis_warning = .data$season_k_index < 0.90 &
      .data$season_k_check_p_value < 0.05,
    frequent_residual_autocorrelation =
      .data$share_seasons_abs_lag1_acf_gt_0_3 > 0.25
  )

rollup <- diagnostics %>%
  group_by(.data$analysis) %>%
  summarise(
    n_fits = n(),
    min_deviance_explained = min(.data$deviance_explained, na.rm = TRUE),
    median_deviance_explained = median(.data$deviance_explained, na.rm = TRUE),
    median_scale = median(.data$scale, na.rm = TRUE),
    max_scale = max(.data$scale, na.rm = TRUE),
    min_global_k_index = min(.data$global_k_index, na.rm = TRUE),
    min_season_k_index = min(.data$season_k_index, na.rm = TRUE),
    share_global_basis_warning = mean(.data$global_basis_warning, na.rm = TRUE),
    share_season_basis_warning = mean(.data$season_basis_warning, na.rm = TRUE),
    median_lag1_residual_acf = median(.data$median_lag1_residual_acf, na.rm = TRUE),
    share_fits_with_frequent_residual_autocorrelation = mean(
      .data$frequent_residual_autocorrelation,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(.data$analysis)

flags <- diagnostics %>%
  mutate(
    low_deviance_explained = .data$deviance_explained < 0.90,
    residual_autocorrelation_warning = .data$frequent_residual_autocorrelation
  ) %>%
  filter(
    .data$low_deviance_explained |
      .data$global_basis_warning |
      .data$season_basis_warning |
      .data$residual_autocorrelation_warning
  ) %>%
  arrange(.data$analysis, .data$deviance_explained, .data$global_k_index)

output_dir <- "outputs/diagnostics/gam_fit"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(diagnostics, file.path(output_dir, "gam_fit_diagnostics_combined.csv"))
write_csv(rollup, file.path(output_dir, "gam_fit_diagnostics_rollup.csv"))
write_csv(flags, file.path(output_dir, "gam_fit_diagnostics_flags.csv"))

cat("Wrote GAM fit diagnostic audit tables to ", output_dir, "\n", sep = "")
print(rollup, n = Inf)
