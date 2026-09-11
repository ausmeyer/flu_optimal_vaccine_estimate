#!/usr/bin/env Rscript

# Build the supplementary manuscript figures from completed production outputs.

font_cache <- file.path(tempdir(), "fontconfig-cache")
dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(mgcv)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})

figure_dir <- "outputs/primary/final_figures"
source_data_dir <- "outputs/primary/figure_source_data"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)

figure_stems <- c(
  "figure_s1_gam_fit_diagnostics",
  "figure_s2_alternative_outcome_relative_regret"
)
figure_names <- c(paste0(figure_stems, ".pdf"), paste0(figure_stems, ".png"))
figure_paths <- file.path(figure_dir, figure_names)
source_names <- paste0(figure_stems, "_source_data.csv")
source_paths <- file.path(source_data_dir, source_names)

staging_dir <- tempfile("supplement-figures-")
dir.create(staging_dir, recursive = TRUE)
on.exit(unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)
staged_figure_paths <- file.path(staging_dir, figure_names)
staged_source_paths <- file.path(staging_dir, source_names)

palette <- list(
  ink = "#252525",
  gray = "#8C8C8C",
  light_gray = "#D9D9D9",
  blue = "#2A5C85",
  light_blue = "#A9C7DB",
  pale_blue = "#DCE9F1",
  orange = "#D88732"
)

week_axis <- function(week) {
  week <- as.integer(week)
  dplyr::case_when(
    week >= 36L & week <= 52L ~ as.numeric(week - 35L),
    week == 53L ~ 17.5,
    week >= 1L & week <= 22L ~ as.numeric(week + 17L),
    TRUE ~ NA_real_
  )
}

axis_break_weeks <- c(36L, 40L, 44L, 48L, 52L, 4L, 8L, 12L, 16L, 20L)
decision_break_weeks <- c(36L, 40L, 44L, 48L, 52L, 4L, 8L, 12L)

scale_x_flu_season <- function(decision_window = FALSE) {
  breaks <- if (decision_window) decision_break_weeks else axis_break_weeks
  ggplot2::scale_x_continuous(
    breaks = week_axis(breaks),
    labels = breaks,
    expand = expansion(mult = c(0.01, 0.01))
  )
}

theme_manuscript <- function(base_size = 9) {
  theme_light(base_family = "Helvetica", base_size = base_size) +
    theme(
      text = element_text(color = palette$ink),
      axis.text = element_text(color = palette$ink),
      axis.title = element_text(face = "plain"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#E6E6E6", linewidth = 0.25),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = palette$ink),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.key.height = grid::unit(7, "pt"),
      legend.key.width = grid::unit(17, "pt"),
      plot.margin = margin(5, 6, 5, 5)
    )
}

save_pdf <- function(plot, path, width, height) {
  ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )
}

save_png <- function(plot, path, width, height) {
  ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    device = "png",
    bg = "white",
    limitsize = FALSE
  )
}

require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required input(s):\n", paste(missing, collapse = "\n"))
  }
}


# Figure S1: targeted GAM-fit diagnostics -----------------------------------

primary_manifest_path <- "outputs/primary/tables/gam_primary_analysis_manifest.json"
higher_basis_manifest_path <- paste0(
  "outputs/sensitivity/gam_smooth_k16_k14_primary/tables/",
  "gam_smooth_k16_k14_primary_analysis_manifest.json"
)
state_population_path <- "data/processed/state_population.csv"
require_files(c(
  primary_manifest_path,
  higher_basis_manifest_path,
  state_population_path
))

primary_manifest <- jsonlite::read_json(primary_manifest_path, simplifyVector = TRUE)
higher_basis_manifest <- jsonlite::read_json(
  higher_basis_manifest_path,
  simplifyVector = TRUE
)
require_files(c(primary_manifest$model_path, higher_basis_manifest$model_path))

if (
  primary_manifest$configuration$gam_global_k != 14L ||
    primary_manifest$configuration$gam_season_k != 12L ||
    higher_basis_manifest$configuration$gam_global_k != 16L ||
    higher_basis_manifest$configuration$gam_season_k != 14L
) {
  stop("Unexpected spline-basis configuration for the Figure S1 inputs.")
}

target_seasons <- c("2018/19", "2024/25")
state_population <- read_csv(state_population_path, show_col_types = FALSE) %>%
  select("state", "population")

extract_state_fits <- function(models, fitted_name, include_observed = FALSE) {
  purrr::imap_dfr(models, function(model, state_name) {
    fitted_values <- as.numeric(stats::predict(
      model$fit,
      newdata = model$data,
      type = "response"
    ))
    if (length(fitted_values) != nrow(model$data)) {
      stop("Prediction length mismatch for ", state_name, ".")
    }

    output <- model$data %>%
      transmute(
        state = as.character(.data$state),
        season = as.character(.data$season),
        mmwr_year = as.integer(.data$mmwr_year),
        mmwr_week = as.integer(.data$mmwr_week),
        week_start = as.Date(.data$week_start),
        fitted_value = fitted_values
      )
    names(output)[names(output) == "fitted_value"] <- fitted_name

    if (include_observed) {
      output$observed_ili_proportion <- as.numeric(model$data$ili_prop)
    }
    output
  }) %>%
    filter(.data$season %in% target_seasons)
}

primary_state_fits <- extract_state_fits(
  readRDS(primary_manifest$model_path),
  fitted_name = "primary_gam_fit_ili_proportion",
  include_observed = TRUE
)
higher_basis_state_fits <- extract_state_fits(
  readRDS(higher_basis_manifest$model_path),
  fitted_name = "higher_basis_gam_fit_ili_proportion"
)

state_keys <- c("state", "season", "mmwr_year", "mmwr_week", "week_start")
if (
  anyDuplicated(primary_state_fits[state_keys]) ||
    anyDuplicated(higher_basis_state_fits[state_keys])
) {
  stop("Figure S1 model inputs contain duplicate state-season-week rows.")
}

fit_source_data <- primary_state_fits %>%
  inner_join(
    higher_basis_state_fits,
    by = state_keys,
    relationship = "one-to-one"
  ) %>%
  left_join(state_population, by = "state", relationship = "many-to-one")

if (
  n_distinct(fit_source_data$state) != 51L ||
    !setequal(unique(fit_source_data$season), target_seasons) ||
    anyNA(fit_source_data$population)
) {
  stop("Figure S1 inputs do not cover 51 jurisdictions and both target seasons.")
}

fit_source_data <- fit_source_data %>%
  summarise(
    observed_ili_proportion = weighted.mean(
      .data$observed_ili_proportion,
      .data$population,
      na.rm = TRUE
    ),
    primary_gam_fit_ili_proportion = weighted.mean(
      .data$primary_gam_fit_ili_proportion,
      .data$population,
      na.rm = TRUE
    ),
    higher_basis_gam_fit_ili_proportion = weighted.mean(
      .data$higher_basis_gam_fit_ili_proportion,
      .data$population,
      na.rm = TRUE
    ),
    .by = c("season", "mmwr_year", "mmwr_week", "week_start")
  ) %>%
  mutate(
    flu_season_week_index = week_axis(.data$mmwr_week),
    season = factor(.data$season, levels = target_seasons)
  ) %>%
  arrange(.data$season, .data$flu_season_week_index) %>%
  mutate(season = as.character(.data$season))

fit_values <- fit_source_data %>%
  select(
    "observed_ili_proportion",
    "primary_gam_fit_ili_proportion",
    "higher_basis_gam_fit_ili_proportion"
  )
if (
  nrow(fit_source_data) != length(target_seasons) * 39L ||
    any(!is.finite(as.matrix(fit_values))) ||
    any(as.matrix(fit_values) < 0 | as.matrix(fit_values) > 1)
) {
  stop("Figure S1 source data failed row-count or value-range validation.")
}

write_csv(fit_source_data, staged_source_paths[[1]])

fit_plot_data <- fit_source_data %>%
  pivot_longer(
    cols = c(
      "observed_ili_proportion",
      "primary_gam_fit_ili_proportion",
      "higher_basis_gam_fit_ili_proportion"
    ),
    names_to = "series",
    values_to = "ili_proportion"
  ) %>%
  mutate(
    series = factor(
      .data$series,
      levels = c(
        "observed_ili_proportion",
        "primary_gam_fit_ili_proportion",
        "higher_basis_gam_fit_ili_proportion"
      ),
      labels = c("Observed", "Primary GAM fit", "Higher-basis GAM fit")
    ),
    season = factor(.data$season, levels = target_seasons)
  )

s1_colors <- c(
  "Observed" = palette$gray,
  "Primary GAM fit" = palette$blue,
  "Higher-basis GAM fit" = palette$orange
)
s1_linetypes <- c(
  "Observed" = "solid",
  "Primary GAM fit" = "solid",
  "Higher-basis GAM fit" = "22"
)

figure_s1 <- ggplot(
  fit_plot_data,
  aes(
    x = .data$flu_season_week_index,
    y = .data$ili_proportion,
    color = .data$series,
    linetype = .data$series
  )
) +
  geom_line(
    aes(linewidth = .data$series),
    lineend = "round"
  ) +
  geom_point(
    data = fit_plot_data %>% filter(.data$series == "Observed"),
    size = 1.15,
    stroke = 0,
    show.legend = FALSE
  ) +
  facet_wrap(vars(.data$season), nrow = 1) +
  scale_color_manual(values = s1_colors, guide = guide_legend(title = NULL)) +
  scale_linetype_manual(values = s1_linetypes, guide = guide_legend(title = NULL)) +
  scale_linewidth_manual(
    values = c(
      "Observed" = 0.38,
      "Primary GAM fit" = 0.82,
      "Higher-basis GAM fit" = 0.76
    ),
    guide = "none"
  ) +
  scale_x_flu_season() +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    x = "MMWR week",
    y = "ILI visits (%)"
  ) +
  theme_manuscript(base_size = 9.2) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    panel.spacing = grid::unit(10, "pt")
  )

save_pdf(figure_s1, staged_figure_paths[[1]], width = 8.2, height = 3.5)
save_png(figure_s1, staged_figure_paths[[3]], width = 8.2, height = 3.5)


# Figure S2: relative regret under alternative burden curves ----------------

alternative_curve_key <- tribble(
  ~analysis, ~panel_label, ~near_optimal_path,
  "empirical_primary",
  "Unsmoothed empirical ILINet",
  paste0(
    "outputs/sensitivity/empirical_primary/tables/",
    "empirical_primary_near_optimal_curve.csv"
  ),
  "gam_nrevss_national",
  "ICLS/NREVSS all-influenza positivity",
  paste0(
    "outputs/sensitivity/gam_nrevss_national/tables/",
    "gam_nrevss_national_near_optimal_curve.csv"
  ),
  "empirical_ili_plus_all_seasons_weighted",
  "Population-weighted ILI+",
  paste0(
    "outputs/sensitivity/empirical_ili_plus_all_seasons_weighted/tables/",
    "empirical_ili_plus_all_seasons_weighted_near_optimal_curve.csv"
  )
)
empirical_manifest_path <- paste0(
  "outputs/sensitivity/empirical_primary/tables/",
  "empirical_primary_analysis_manifest.json"
)
require_files(c(alternative_curve_key$near_optimal_path, empirical_manifest_path))

empirical_manifest <- jsonlite::read_json(
  empirical_manifest_path,
  simplifyVector = TRUE
)
if (
  !identical(
    as.character(empirical_manifest$configuration$empirical_curve_smoothing),
    "none"
  ) ||
    !identical(
      as.character(empirical_manifest$configuration$timing_shift_interpolation),
      "linear"
    )
) {
  stop(
    "Figure S2 requires an unsmoothed empirical_primary analysis with ",
    "linear interpolation used only for continuous timing shifts."
  )
}

relative_regret_columns <- c(
  "mean_relative_regret",
  "median_relative_regret",
  "lower_50_relative_regret",
  "upper_50_relative_regret",
  "lower_95_relative_regret",
  "upper_95_relative_regret"
)

read_relative_regret <- function(analysis, panel_label, near_optimal_path) {
  curve <- read_csv(near_optimal_path, show_col_types = FALSE) %>%
    filter(.data$state == "US", .data$age_group == "all")

  required_columns <- c(
    "vaccination_week",
    "relative_regret_tolerance",
    relative_regret_columns
  )
  missing_columns <- setdiff(required_columns, names(curve))
  if (length(missing_columns) > 0L) {
    stop(
      near_optimal_path,
      " is missing column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  # Regret summaries are repeated for each near-optimality tolerance. Verify
  # that they are identical before retaining one row per vaccination week.
  consistency <- curve %>%
    summarise(
      across(all_of(relative_regret_columns), n_distinct),
      .by = "vaccination_week"
    )
  if (any(as.matrix(consistency[relative_regret_columns]) != 1L)) {
    stop("Relative-regret summaries vary by tolerance in ", near_optimal_path)
  }

  curve %>%
    select("vaccination_week", all_of(relative_regret_columns)) %>%
    distinct() %>%
    mutate(
      analysis = analysis,
      panel_label = panel_label,
      vaccination_week = as.integer(.data$vaccination_week),
      flu_season_week_index = week_axis(.data$vaccination_week)
    ) %>%
    select(
      "analysis",
      "panel_label",
      "vaccination_week",
      "flu_season_week_index",
      all_of(relative_regret_columns)
    )
}

relative_regret_source_data <- purrr::pmap_dfr(
  alternative_curve_key,
  read_relative_regret
) %>%
  mutate(
    panel_label = factor(
      .data$panel_label,
      levels = alternative_curve_key$panel_label
    )
  ) %>%
  arrange(.data$panel_label, .data$flu_season_week_index)

interval_values <- relative_regret_source_data %>%
  select(all_of(relative_regret_columns))
interval_order_valid <- with(
  relative_regret_source_data,
  lower_95_relative_regret <= lower_50_relative_regret &
    lower_50_relative_regret <= upper_50_relative_regret &
    upper_50_relative_regret <= upper_95_relative_regret
)
curve_counts <- relative_regret_source_data %>%
  summarise(
    n_weeks = n(),
    valid_weeks = setequal(.data$vaccination_week, c(36:52, 1:12)),
    .by = c("analysis", "panel_label")
  )
if (
  nrow(curve_counts) != 3L ||
    any(curve_counts$n_weeks != 29L) ||
    any(!curve_counts$valid_weeks) ||
    any(!is.finite(as.matrix(interval_values))) ||
    any(as.matrix(interval_values) < 0 | as.matrix(interval_values) > 1) ||
    any(!interval_order_valid)
) {
  stop("Figure S2 source data failed coverage, range, or interval validation.")
}

write_csv(
  relative_regret_source_data %>% mutate(panel_label = as.character(.data$panel_label)),
  staged_source_paths[[2]]
)

# Do not mark an "optimal" week on this figure. The minimum of mean relative
# regret need not equal the fixed decision selected by minimum mean absolute
# regret, which is the decision criterion reported in the manuscript.
figure_s2 <- ggplot(
  relative_regret_source_data,
  aes(x = .data$flu_season_week_index, y = .data$mean_relative_regret)
) +
  annotate(
    "rect",
    xmin = week_axis(36L) - 0.5,
    xmax = week_axis(44L) + 0.5,
    ymin = -Inf,
    ymax = Inf,
    fill = palette$light_gray,
    alpha = 0.45
  ) +
  geom_vline(
    xintercept = week_axis(c(44L, 48L, 52L)),
    color = "#C7C7C7",
    linewidth = 0.32,
    linetype = "22"
  ) +
  geom_ribbon(
    aes(
      ymin = .data$lower_95_relative_regret,
      ymax = .data$upper_95_relative_regret,
      fill = "Central 95% interval"
    ),
    color = NA,
    alpha = 0.58
  ) +
  geom_ribbon(
    aes(
      ymin = .data$lower_50_relative_regret,
      ymax = .data$upper_50_relative_regret,
      fill = "Central 50% interval"
    ),
    color = NA,
    alpha = 0.68
  ) +
  geom_line(
    aes(color = "Mean relative regret"),
    linewidth = 0.82,
    lineend = "round"
  ) +
  facet_wrap(vars(.data$panel_label), ncol = 1) +
  scale_fill_manual(
    values = c(
      "Central 95% interval" = palette$pale_blue,
      "Central 50% interval" = palette$light_blue
    ),
    breaks = c("Central 50% interval", "Central 95% interval")
  ) +
  scale_color_manual(values = c("Mean relative regret" = palette$blue)) +
  scale_x_flu_season(decision_window = TRUE) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 10),
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Vaccination week (MMWR)",
    y = "Achievable modeled protection left unused"
  ) +
  guides(
    color = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  ) +
  theme_manuscript(base_size = 9.2) +
  theme(
    legend.position = "top",
    legend.justification = "right",
    panel.spacing = grid::unit(7, "pt")
  )

save_pdf(figure_s2, staged_figure_paths[[2]], width = 7.4, height = 6.8)
save_png(figure_s2, staged_figure_paths[[4]], width = 7.4, height = 6.8)


# Publish only this script's outputs after every asset passes validation. ------

staged_paths <- c(staged_figure_paths, staged_source_paths)
if (!all(file.exists(staged_paths)) || any(file.info(staged_paths)$size <= 0L)) {
  stop("One or more staged supplementary figure assets are missing or empty.")
}

owned_paths <- c(figure_paths, source_paths)
owned_existing_paths <- owned_paths[file.exists(owned_paths)]
if (length(owned_existing_paths) > 0L) {
  unlink(owned_existing_paths, force = TRUE)
}
if (!all(file.copy(staged_figure_paths, figure_paths, overwrite = TRUE)) ||
    !all(file.copy(staged_source_paths, source_paths, overwrite = TRUE))) {
  stop("Could not publish one or more supplementary figure assets.")
}
if (!all(file.exists(owned_paths)) || any(file.info(owned_paths)$size <= 0L)) {
  stop("One or more published supplementary figure assets are missing or empty.")
}

message("Wrote supplementary manuscript figures:")
message(paste(figure_paths, collapse = "\n"))
message("Wrote supplementary figure source data:")
message(paste(source_paths, collapse = "\n"))
