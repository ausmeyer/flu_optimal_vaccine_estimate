#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

source("R/calendar.R")

summary_path <- "outputs/sensitivity/summary/tables/robustness_suite_us_summary.csv"
final_figure_dir <- "outputs/primary/final_figures"
figure_names <- c(
  "figure_4_robustness_summary.pdf",
  "figure_4_robustness_summary.png"
)
final_paths <- file.path(final_figure_dir, figure_names)
dir.create(final_figure_dir, recursive = TRUE, showWarnings = FALSE)

staging_dir <- tempfile("robustness-figure-")
dir.create(staging_dir, recursive = TRUE)
on.exit(unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)
staged_paths <- file.path(staging_dir, figure_names)

if (!file.exists(summary_path)) {
  stop("Missing robustness summary: ", summary_path)
}

analysis_key <- tibble::tribble(
  ~analysis, ~analysis_label,
  "gam_primary", "Primary ILINet GAM",
  "gam_no_timing_shift_primary", "No epidemic timing shift",
  "empirical_primary", "Unsmoothed empirical ILINet",
  "gam_precovid_primary", "Pre-pandemic burden",
  "gam_fixed_immune_lag_primary", "Fixed 2-week immune lag",
  "gam_slow_waning_primary", "Ray slower waning",
  "gam_fast_waning_primary", "Ray faster waning",
  "gam_spencer_ferdinands_fast_waning", "Spencer-Ferdinands fast waning",
  "gam_spencer_ferdinands_onset_aligned",
    "Spencer-Ferdinands, onset aligned",
  "gam_spencer_slow_waning", "Spencer slower waning",
  "gam_smooth_k16_k14_primary", "Higher spline basis",
  "gam_ve_reference_4w", "VE reference at 4 weeks",
  "gam_ve_reference_12w", "VE reference at 12 weeks",
  "gam_nrevss_national", "ICLS/NREVSS all-influenza positivity",
  "gam_nrevss_post2015_flu_a", "ICLS/NREVSS influenza A, 2015/16 onward",
  "empirical_ili_plus_all_seasons_weighted",
    "ILI+, population-weighted ILI",
  "empirical_ili_plus_all_seasons_unweighted",
    "ILI+, unweighted ILI",
  "empirical_ili_plus_clinical_only_weighted",
    "ILI+, clinical labs only"
)

robustness <- read_csv(summary_path, show_col_types = FALSE) %>%
  filter(.data$state == "US", .data$age_group == "all") %>%
  inner_join(analysis_key, by = "analysis", relationship = "many-to-one")

missing <- setdiff(analysis_key$analysis, robustness$analysis)
if (length(missing) > 0L) {
  stop("Robustness figure missing analysis output(s): ", paste(missing, collapse = ", "))
}
if (anyNA(robustness$min_mean_regret_week) ||
    anyNA(robustness$lower_50_week) ||
    anyNA(robustness$upper_50_week) ||
    anyNA(robustness$lower_95_week) ||
    anyNA(robustness$upper_95_week)) {
  stop("Robustness figure inputs contain missing timing summaries.")
}

plot_data <- robustness %>%
  mutate(
    analysis_label = factor(
      .data$analysis_label,
      levels = rev(analysis_key$analysis_label)
    ),
    selected_index = week_to_season_index(.data$min_mean_regret_week),
    lower_50_index = week_to_season_index(.data$lower_50_week),
    upper_50_index = week_to_season_index(.data$upper_50_week),
    lower_95_index = week_to_season_index(.data$lower_95_week),
    upper_95_index = week_to_season_index(.data$upper_95_week)
  )

axis_weeks <- c(40L, 44L, 48L, 52L, 4L)
plot <- ggplot(plot_data, aes(y = .data$analysis_label)) +
  annotate(
    "rect",
    xmin = week_to_season_index(36L) - 0.5,
    xmax = week_to_season_index(44L) + 0.5,
    ymin = -Inf,
    ymax = Inf,
    fill = "#EEEEEE"
  ) +
  geom_segment(
    aes(
      x = .data$lower_95_index,
      xend = .data$upper_95_index,
      yend = .data$analysis_label
    ),
    color = "#9ECAE1",
    linewidth = 1.2,
    lineend = "round"
  ) +
  geom_segment(
    aes(
      x = .data$lower_50_index,
      xend = .data$upper_50_index,
      yend = .data$analysis_label
    ),
    color = "#3182BD",
    linewidth = 3.2,
    lineend = "round"
  ) +
  geom_point(
    aes(x = .data$selected_index),
    shape = 21,
    size = 2.8,
    stroke = 0.8,
    color = "#8C2D04",
    fill = "white"
  ) +
  scale_x_continuous(
    breaks = week_to_season_index(axis_weeks),
    labels = axis_weeks,
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = "Vaccination week (MMWR)",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "#252525")
  )

ggsave(
  staged_paths[[1]],
  plot,
  width = 9.2,
  height = 7.4,
  device = grDevices::cairo_pdf,
  bg = "white"
)
ggsave(
  staged_paths[[2]],
  plot,
  width = 9.2,
  height = 7.4,
  dpi = 300,
  device = "png",
  bg = "white"
)

if (!all(file.exists(staged_paths)) || any(file.info(staged_paths)$size <= 0L)) {
  stop("One or more staged Figure 4 assets are missing or empty.")
}

owned_existing_paths <- final_paths[file.exists(final_paths)]
if (length(owned_existing_paths) > 0L) {
  unlink(owned_existing_paths, force = TRUE)
}
if (!all(file.copy(staged_paths, final_paths, overwrite = TRUE))) {
  stop("Could not publish Figure 4 to ", final_figure_dir, ".")
}
if (!all(file.exists(final_paths)) || any(file.info(final_paths)$size <= 0L)) {
  stop("One or more published Figure 4 assets are missing or empty.")
}

message(
  "Wrote robustness figure:\n",
  paste(final_paths, collapse = "\n")
)
