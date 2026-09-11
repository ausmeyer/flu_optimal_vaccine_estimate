#!/usr/bin/env Rscript

# Build the final primary manuscript figures and retained figure component from
# completed analysis outputs.

font_cache <- file.path(tempdir(), "fontconfig-cache")
dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
  library(purrr)
  library(readr)
  library(sf)
  library(tibble)
  library(tidyr)
})

source("R/calendar.R")
source("R/hhs_regions.R")

table_dir <- "outputs/primary/tables"
figure_dir <- "outputs/primary/final_figures"
component_dir <- "outputs/primary/figure_components"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(component_dir, recursive = TRUE, showWarnings = FALSE)

publication_figure_stems <- c(
  "figure_1_epidemic_timing_and_model_fit",
  "figure_2_national_regret_curves_by_age",
  "figure_3_state_timing_and_optimality_probabilities"
)
component_figure_stems <- "figure_3a_state_timing_choropleth"
figure_stems <- c(publication_figure_stems, component_figure_stems)
pdf_names <- paste0(figure_stems, ".pdf")
png_names <- paste0(figure_stems, ".png")
figure_names <- c(pdf_names, png_names)
output_dirs <- c(
  rep(figure_dir, length(publication_figure_stems)),
  component_dir,
  rep(figure_dir, length(publication_figure_stems)),
  component_dir
)
figure_paths <- file.path(output_dirs, figure_names)
staging_dir <- tempfile("results-figures-")
dir.create(staging_dir, recursive = TRUE)
on.exit(unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)
staged_paths <- file.path(staging_dir, figure_names)
staged_pdf_paths <- file.path(staging_dir, pdf_names)
staged_png_paths <- file.path(staging_dir, png_names)

palette <- list(
  ink = "#252525",
  gray = "#9A9A9A",
  light_gray = "#D9D9D9",
  blue = "#2A5C85",
  light_blue = "#A9C7DB",
  pale_blue = "#DCE9F1",
  orange = "#D88732",
  green = "#438A72",
  vermillion = "#B6533C"
)

age_files <- tribble(
  ~age_group, ~age_label, ~prefix,
  "all", "All ages", "gam_primary",
  "0-4", "0-4 years", "gam_age_0_4",
  "5-24", "5-24 years", "gam_age_5_24",
  "25-49", "25-49 years", "gam_age_25_49",
  "50-64", "50-64 years", "gam_age_50_64",
  "65+", "65+ years", "gam_age_65_"
) %>%
  mutate(
    age_label = factor(
      .data$age_label,
      levels = c(
        "All ages", "0-4 years", "5-24 years", "25-49 years",
        "50-64 years", "65+ years"
      )
    )
  )

required_files <- c(
  file.path(table_dir, paste0(age_files$prefix, "_regret_curve.csv")),
  file.path(table_dir, paste0(age_files$prefix, "_optimal_week_draws.csv")),
  file.path(table_dir, "gam_primary_analysis_manifest.json"),
  "data/processed/state_population.csv",
  "data/processed/us_states_2024_cartographic_boundaries.rds"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required analysis output(s):\n", paste(missing_files, collapse = "\n"))
}

week_axis <- function(week) {
  week <- as.integer(week)
  case_when(
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
  theme_minimal(base_family = "Helvetica", base_size = base_size) +
    theme(
      text = element_text(color = palette$ink),
      axis.text = element_text(color = palette$ink),
      axis.text.x = element_text(angle = 0, vjust = 0.5),
      axis.title = element_text(face = "plain"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#E6E6E6", linewidth = 0.25),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = palette$ink),
      plot.title = element_text(face = "bold", size = rel(1), hjust = 0),
      plot.tag = element_text(face = "bold", size = rel(1.15)),
      plot.tag.position = c(0, 1),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_text(face = "plain"),
      legend.key.height = grid::unit(7, "pt"),
      legend.key.width = grid::unit(15, "pt"),
      plot.margin = margin(5, 6, 5, 5)
    )
}

# Keep panel letters at 12 points when the figure is placed at manuscript width.
panel_tag_theme <- function(figure_width) {
  theme(
    plot.tag = element_text(face = "bold", size = 12 * figure_width / 6.5,
                            hjust = 0, vjust = 1),
    plot.tag.position = c(0, 1)
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

state_regions <- hhs_region_crosswalk() %>%
  transmute(
    state = .data$state_or_territory,
    hhs_region = factor(
      paste("Region", .data$region_number),
      levels = paste("Region", 1:10)
    ),
    region_number = as.integer(.data$region_number)
  ) %>%
  filter(.data$state %in% c(state.name, "District of Columbia")) %>%
  distinct(.data$state, .keep_all = TRUE)

if (nrow(state_regions) != 51L) {
  stop("Expected HHS-region assignments for 50 states and the District of Columbia.")
}

state_positions <- state_regions %>%
  arrange(.data$region_number, .data$state) %>%
  mutate(state_y = rev(seq_len(n())))

region_positions <- state_positions %>%
  summarise(
    region_y = mean(.data$state_y),
    lower_edge = min(.data$state_y) - 0.5,
    .by = c("hhs_region", "region_number")
  ) %>%
  arrange(.data$region_number)

region_separators <- region_positions %>%
  filter(.data$region_number < max(.data$region_number)) %>%
  pull(.data$lower_edge)

region_separator_data <- tibble(separator_y = region_separators)


# Figure 1: historical epidemic timing and fitted curves ---------------------

manifest <- jsonlite::read_json(
  file.path(table_dir, "gam_primary_analysis_manifest.json"),
  simplifyVector = TRUE
)

if (!file.exists(manifest$model_path)) {
  stop("Primary model listed in the manifest does not exist: ", manifest$model_path)
}

primary_models <- readRDS(manifest$model_path)
state_population <- read_csv(
  "data/processed/state_population.csv",
  show_col_types = FALSE
) %>%
  select("state", "population")

fitted_state_curves <- purrr::imap_dfr(
  primary_models,
  function(model, state_name) {
    model$data %>%
      transmute(
        state = as.character(.data$state),
        season = as.character(.data$season),
        mmwr_week = as.integer(.data$mmwr_week),
        observed_ili = as.numeric(.data$ili_prop),
        fitted_ili = as.numeric(stats::fitted(model$fit))
      )
  }
) %>%
  left_join(state_population, by = "state") %>%
  mutate(week_axis = week_axis(.data$mmwr_week))

if (
  n_distinct(fitted_state_curves$state) != 51L ||
    n_distinct(fitted_state_curves$season) != 12L
) {
  stop("Figure 1 requires 51 jurisdictions and 12 influenza seasons.")
}

if (anyNA(fitted_state_curves$population)) {
  missing_states <- fitted_state_curves %>%
    filter(is.na(.data$population)) %>%
    distinct(.data$state) %>%
    pull(.data$state)
  stop("Missing population for: ", paste(missing_states, collapse = ", "))
}

national_curves <- fitted_state_curves %>%
  group_by(.data$season, .data$mmwr_week, .data$week_axis) %>%
  summarise(
    observed_ili = weighted.mean(.data$observed_ili, .data$population, na.rm = TRUE),
    fitted_ili = weighted.mean(.data$fitted_ili, .data$population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c("observed_ili", "fitted_ili"),
    names_to = "series",
    values_to = "ili_proportion"
  ) %>%
  mutate(
    series = factor(
      .data$series,
      levels = c("observed_ili", "fitted_ili"),
      labels = c("Observed", "GAM fit")
    ),
    season = factor(
      .data$season,
      levels = sprintf(
        "%d/%02d",
        manifest$configuration$primary_season_start_years,
        (manifest$configuration$primary_season_start_years + 1L) %% 100L
      )
    )
  )

season_levels <- levels(national_curves$season)
national_curves <- national_curves %>%
  mutate(season = factor(as.character(.data$season), levels = season_levels))

peak_timing <- fitted_state_curves %>%
  group_by(.data$state, .data$season) %>%
  slice_max(.data$fitted_ili, n = 1L, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(state_positions, by = "state") %>%
  mutate(
    season = factor(.data$season, levels = season_levels),
    season_x = match(.data$season, season_levels)
  )

timing_colors <- c("#2A4D69", "#4B86B4", "#8CB8CF", "#DDD6B5", "#D79B5B", "#A44A3F")
peak_timing_scale_weeks <- c(40L, 44L, 48L, 52L, 4L, 8L, 12L)
peak_timing_scale_limits <- week_axis(c(40L, 12L))
peak_timing_scale_breaks <- week_axis(peak_timing_scale_weeks)

# State-level vaccination decisions occupy a much narrower range than epidemic
# peak weeks. Use the full palette over weeks 40-52 so the meaningful variation
# within Figure 3A remains visible rather than being compressed into blue tones.
decision_timing_scale_weeks <- c(40L, 44L, 48L, 52L)
decision_timing_scale_limits <- week_axis(c(40L, 52L))
decision_timing_scale_breaks <- week_axis(decision_timing_scale_weeks)

p1_curves <- ggplot(
  national_curves,
  aes(x = .data$week_axis, y = .data$ili_proportion, color = .data$series)
) +
  geom_line(aes(linewidth = .data$series), lineend = "round") +
  facet_wrap(vars(.data$season), ncol = 4) +
  scale_color_manual(
    values = c("Observed" = palette$gray, "GAM fit" = palette$blue),
    guide = guide_legend(title = NULL)
  ) +
  scale_linewidth_manual(
    values = c("Observed" = 0.38, "GAM fit" = 0.72),
    guide = "none"
  ) +
  scale_x_flu_season() +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "MMWR week", y = "ILI visits (%)", tag = "A") +
  ggplot2::theme_light(base_family = "Helvetica", base_size = 8.4) +
  theme(
    text = element_text(color = palette$ink),
    axis.text = element_text(color = palette$ink),
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    axis.title = element_text(face = "plain"),
    panel.spacing = grid::unit(7, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = palette$ink),
    plot.tag = element_text(face = "bold", size = rel(1.15)),
    plot.tag.position = c(0, 1),
    legend.position = "top",
    legend.justification = "right",
    legend.title = element_text(face = "plain"),
    legend.key.height = grid::unit(7, "pt"),
    legend.key.width = grid::unit(15, "pt"),
    legend.margin = margin(0, 0, 2, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(5, 6, 5, 5)
  )

p1_timing_heatmap <- ggplot(
  peak_timing,
  aes(x = .data$season_x, y = .data$state_y, fill = .data$week_axis)
) +
  geom_tile(color = "white", linewidth = 0.12, width = 1, height = 1) +
  geom_text(
    data = state_positions,
    aes(x = 0.20, y = .data$state_y, label = .data$state),
    inherit.aes = FALSE,
    hjust = 1,
    size = 5.8 / ggplot2::.pt,
    family = "Helvetica",
    color = palette$ink
  ) +
  geom_text(
    data = region_positions,
    aes(x = -3.45, y = .data$region_y, label = .data$hhs_region),
    inherit.aes = FALSE,
    hjust = 0.5,
    fontface = "bold",
    size = 5.8 / ggplot2::.pt,
    family = "Helvetica",
    color = palette$ink
  ) +
  geom_rect(
    data = region_separator_data,
    aes(
      xmin = -4.5,
      xmax = 12.5,
      ymin = .data$separator_y - 0.09,
      ymax = .data$separator_y + 0.09
    ),
    inherit.aes = FALSE,
    fill = "white",
    color = NA
  ) +
  geom_hline(
    yintercept = region_separators,
    color = palette$ink,
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    colors = timing_colors,
    limits = peak_timing_scale_limits,
    breaks = peak_timing_scale_breaks,
    labels = peak_timing_scale_weeks,
    oob = scales::squish,
    name = "Peak week",
    guide = guide_colorbar(
      barwidth = grid::unit(40, "mm"),
      barheight = grid::unit(3, "mm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  scale_x_continuous(
    breaks = seq_along(season_levels),
    labels = season_levels,
    limits = c(-4.5, 12.5),
    expand = expansion(add = 0)
  ) +
  scale_y_continuous(
    limits = c(0.5, 51.5),
    expand = expansion(add = 0)
  ) +
  labs(x = "Season", y = NULL, tag = "B") +
  theme_manuscript(base_size = 7.1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(hjust = 0.65),
    panel.grid = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "top",
    legend.justification = "right",
    legend.key.width = grid::unit(24, "pt")
  )

figure_1 <- free(p1_curves, side = "t") +
  free(p1_timing_heatmap, side = "t") +
  plot_layout(widths = c(1.75, 1), guides = "keep")
figure_1 <- figure_1 & panel_tag_theme(11.0)

save_pdf(figure_1, staged_pdf_paths[[1]], width = 11.0, height = 8.7)
save_png(figure_1, staged_png_paths[[1]], width = 11.0, height = 8.7)


# Figure 2: national regret curves by age -----------------------------------

regret_data <- age_files %>%
  select("age_group", "age_label", "prefix") %>%
  pmap_dfr(function(age_group, age_label, prefix) {
    read_csv(
      file.path(table_dir, paste0(prefix, "_regret_curve.csv")),
      show_col_types = FALSE
    ) %>%
      filter(.data$state == "US") %>%
      transmute(
        age_group = age_group,
        age_label = age_label,
        vaccination_week = as.integer(.data$vaccination_week),
        week_axis = week_axis(.data$vaccination_week),
        mean_regret = 10000 * .data$mean_regret,
        lower_50 = 10000 * .data$lower_50,
        upper_50 = 10000 * .data$upper_50,
        lower_95 = 10000 * .data$lower_95,
        upper_95 = 10000 * .data$upper_95
      )
  }) %>%
  mutate(
    age_label = factor(as.character(.data$age_label), levels = levels(age_files$age_label))
  )

regret_minima <- regret_data %>%
  group_by(.data$age_group, .data$age_label) %>%
  slice_min(.data$mean_regret, n = 1L, with_ties = FALSE) %>%
  ungroup()

if (
  n_distinct(regret_data$age_group) != 6L ||
    nrow(regret_minima) != 6L ||
    any(count(regret_data, .data$age_group)$n != 29L)
) {
  stop("Figure 2 requires 29 candidate weeks for each of six age groups.")
}

regret_panel <- function(age, show_x_title = TRUE, show_y_title = TRUE) {
  data <- regret_data %>% filter(.data$age_group == age)
  minimum <- regret_minima %>% filter(.data$age_group == age)
  title <- as.character(unique(data$age_label))

  ggplot(data, aes(x = .data$week_axis, y = .data$mean_regret)) +
    annotate(
      "rect",
      xmin = week_axis(36L) - 0.5,
      xmax = week_axis(44L) + 0.5,
      ymin = -Inf,
      ymax = Inf,
      fill = palette$light_gray,
      alpha = 0.5
    ) +
    geom_vline(
      xintercept = week_axis(c(44L, 48L, 52L)),
      color = "#C7C7C7",
      linewidth = 0.32,
      linetype = "22"
    ) +
    geom_ribbon(
      aes(ymin = .data$lower_95, ymax = .data$upper_95),
      fill = palette$pale_blue,
      color = NA,
      alpha = 0.5
    ) +
    geom_ribbon(
      aes(ymin = .data$lower_50, ymax = .data$upper_50),
      fill = palette$light_blue,
      color = NA,
      alpha = 0.5
    ) +
    geom_line(
      color = palette$blue,
      linewidth = 0.72,
      lineend = "round",
      alpha = 0.5
    ) +
    geom_point(
      data = minimum,
      color = palette$vermillion,
      fill = "white",
      shape = 21,
      stroke = 0.75,
      size = 2.2
    ) +
    annotate(
      "text",
      x = week_axis(40L),
      y = Inf,
      label = "ACIP timing\nrecommendation",
      vjust = 1.1,
      size = 2.1,
      lineheight = 0.9,
      family = "Helvetica",
      color = "#666666"
    ) +
    geom_label(
      data = minimum,
      aes(label = .data$vaccination_week),
      nudge_y = max(data$upper_95, na.rm = TRUE) * 0.055,
      linewidth = 0,
      label.padding = grid::unit(1.2, "pt"),
      fill = scales::alpha("white", 0.85),
      color = palette$ink,
      size = 2.5,
      family = "Helvetica"
    ) +
    scale_x_flu_season(decision_window = TRUE) +
    scale_y_continuous(
      labels = scales::label_number(accuracy = 1, big.mark = ","),
      expand = expansion(mult = c(0, 0.10))
    ) +
    labs(
      x = if (show_x_title) "Vaccination week (MMWR)" else NULL,
      y = if (show_y_title) {
        paste(
          "ILI visits not averted vs best week",
          "(lower is better; per 10,000 weekly encounters)",
          sep = "\n"
        )
      } else {
        NULL
      },
      title = title
    ) +
    ggplot2::theme_light(base_family = "Helvetica", base_size = 8.5) +
    theme(
      text = element_text(color = palette$ink),
      axis.text = element_text(color = palette$ink),
      axis.text.x = element_text(angle = 0, vjust = 0.5),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.title.position = "panel",
      plot.margin = margin(5, 6, 5, 5)
    )
}

p2_all <- regret_panel("all", show_x_title = FALSE, show_y_title = TRUE)
p2_0_4 <- regret_panel("0-4", show_x_title = FALSE, show_y_title = FALSE)
p2_5_24 <- regret_panel("5-24", show_x_title = FALSE, show_y_title = FALSE)
p2_25_49 <- regret_panel("25-49", show_x_title = TRUE, show_y_title = TRUE)
p2_50_64 <- regret_panel("50-64", show_x_title = TRUE, show_y_title = FALSE)
p2_65 <- regret_panel("65+", show_x_title = TRUE, show_y_title = FALSE)

figure_2 <- wrap_plots(
  p2_all,
  p2_0_4,
  p2_5_24,
  p2_25_49,
  p2_50_64,
  p2_65,
  ncol = 3,
  nrow = 2
)

save_pdf(figure_2, staged_pdf_paths[[2]], width = 10.4, height = 8.2)
save_png(figure_2, staged_png_paths[[2]], width = 10.4, height = 8.2)


# Figure 3: state timing and cumulative probability -------------------------

state_timing <- regret_data %>%
  select("age_group", "age_label") %>%
  distinct() %>%
  left_join(age_files %>% select("age_group", "prefix"), by = "age_group") %>%
  pmap_dfr(function(age_group, age_label, prefix) {
    read_csv(
      file.path(table_dir, paste0(prefix, "_regret_curve.csv")),
      show_col_types = FALSE
    ) %>%
      filter(.data$state != "US") %>%
      group_by(.data$state) %>%
      slice_min(.data$mean_regret, n = 1L, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(
        state = .data$state,
        age_group = age_group,
        age_label = age_label,
        regret_minimizing_week = as.integer(.data$vaccination_week)
      )
  }) %>%
  left_join(state_positions, by = "state") %>%
  mutate(
    age_label = factor(as.character(.data$age_label), levels = levels(age_files$age_label)),
    age_x = match(.data$age_label, levels(age_files$age_label)),
    regret_minimizing_week_axis = week_axis(.data$regret_minimizing_week),
    text_color = if_else(
      .data$regret_minimizing_week <= 45L,
      "white",
      palette$ink
    )
  )

if (nrow(state_timing) != 51L * 6L || anyNA(state_timing$state_y)) {
  stop("Figure 3 requires one timing estimate per jurisdiction and age group.")
}
if (any(
  state_timing$regret_minimizing_week < min(decision_timing_scale_weeks) |
    state_timing$regret_minimizing_week > max(decision_timing_scale_weeks)
)) {
  stop("Figure 3 contains a selected week outside its week 40-52 color scale.")
}

state_boundaries <- readRDS(
  "data/processed/us_states_2024_cartographic_boundaries.rds"
)
if (!inherits(state_boundaries, "sf") || nrow(state_boundaries) != 51L) {
  stop("The state boundary file must contain 50 states and the District of Columbia.")
}

p3_timing_heatmap <- ggplot(
  state_timing,
  aes(
    x = .data$age_x,
    y = .data$state_y,
    fill = .data$regret_minimizing_week_axis
  )
) +
  geom_tile(color = "white", linewidth = 0.20, width = 1, height = 1) +
  geom_text(
    data = state_positions,
    aes(x = 0.20, y = .data$state_y, label = .data$state),
    inherit.aes = FALSE,
    hjust = 1,
    size = 6.2 / ggplot2::.pt,
    family = "Helvetica",
    color = palette$ink
  ) +
  geom_text(
    data = region_positions,
    aes(x = -3.55, y = .data$region_y, label = .data$hhs_region),
    inherit.aes = FALSE,
    hjust = 0.5,
    fontface = "bold",
    size = 6 / ggplot2::.pt,
    family = "Helvetica",
    color = palette$ink
  ) +
  geom_rect(
    data = region_separator_data,
    aes(
      xmin = -4.25,
      xmax = 6.5,
      ymin = .data$separator_y - 0.09,
      ymax = .data$separator_y + 0.09
    ),
    inherit.aes = FALSE,
    fill = "white",
    color = NA
  ) +
  geom_hline(
    yintercept = region_separators,
    color = palette$ink,
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = .data$regret_minimizing_week, color = .data$text_color),
    size = 1.95,
    family = "Helvetica"
  ) +
  scale_fill_gradientn(
    colors = timing_colors,
    limits = decision_timing_scale_limits,
    breaks = decision_timing_scale_breaks,
    labels = decision_timing_scale_weeks,
    oob = scales::squish,
    name = "MMWR week",
    guide = guide_colorbar(
      barwidth = grid::unit(40, "mm"),
      barheight = grid::unit(3, "mm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  scale_color_identity() +
  scale_x_continuous(
    breaks = seq_along(levels(age_files$age_label)),
    labels = levels(age_files$age_label),
    limits = c(-4.25, 6.5),
    position = "top",
    expand = expansion(add = 0)
  ) +
  scale_y_continuous(
    limits = c(0.5, 51.5),
    expand = expansion(add = 0)
  ) +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_manuscript(base_size = 7.5) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 7),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "top",
    legend.justification = "right",
    legend.key.width = grid::unit(30, "pt")
  )

cutoff_weeks <- tibble(
  cutoff_week = c(44L, 48L, 52L),
  cutoff_index = week_to_season_index(c(44L, 48L, 52L))
)

cutoff_probabilities <- age_files %>%
  select("age_group", "age_label", "prefix") %>%
  pmap_dfr(function(age_group, age_label, prefix) {
    draws <- read_csv(
      file.path(table_dir, paste0(prefix, "_optimal_week_draws.csv")),
      show_col_types = FALSE
    ) %>%
      filter(.data$state == "US")

    tidyr::crossing(draws, cutoff_weeks) %>%
      summarise(
        probability = mean(.data$optimal_season_index <= .data$cutoff_index),
        .by = c("cutoff_week")
      ) %>%
      mutate(age_group = age_group, age_label = age_label)
  }) %>%
  mutate(
    age_label = factor(
      as.character(.data$age_label),
      levels = rev(levels(age_files$age_label))
    ),
    cutoff_label = factor(
      paste("Week", .data$cutoff_week),
      levels = paste("Week", c(44L, 48L, 52L))
    ),
    probability_label = scales::percent(.data$probability, accuracy = 0.1),
    label_hjust = if_else(.data$probability >= 0.92, 1.12, -0.12),
    label_color = if_else(.data$probability >= 0.92, "white", palette$ink)
  )

cutoff_checks <- cutoff_probabilities %>%
  arrange(.data$age_group, .data$cutoff_week) %>%
  summarise(
    n_cutoffs = n(),
    nondecreasing = all(diff(.data$probability) >= 0),
    .by = "age_group"
  )
if (nrow(cutoff_probabilities) != 18L || any(cutoff_checks$n_cutoffs != 3L) ||
    !all(cutoff_checks$nondecreasing)) {
  stop("Cumulative week 44, 48, and 52 probabilities failed validation.")
}

cutoff_colors <- c(
  "Week 44" = "#377EB8",
  "Week 48" = "#D88732",
  "Week 52" = "#438A72"
)
cutoff_dodge <- position_dodge(width = 0.72, reverse = TRUE)

p3_probabilities <- ggplot(
  cutoff_probabilities,
  aes(
    x = .data$probability,
    y = .data$age_label,
    fill = .data$cutoff_label
  )
) +
  geom_col(
    position = cutoff_dodge,
    orientation = "y",
    width = 0.62,
    color = "white",
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      label = .data$probability_label,
      hjust = .data$label_hjust,
      color = .data$label_color
    ),
    position = cutoff_dodge,
    size = 2.35,
    family = "Helvetica",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = cutoff_colors, name = NULL) +
  scale_color_identity() +
  scale_x_continuous(
    labels = scales::label_percent(accuracy = 1),
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(x = "Probability the optimal week had occurred", y = NULL, tag = "B") +
  theme_manuscript(base_size = 8.5) +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.25),
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    legend.justification = "right"
  )

figure_3 <- free(p3_timing_heatmap, side = "t") +
  free(p3_probabilities, side = "t") +
  plot_layout(widths = c(1.55, 1), guides = "keep")
figure_3 <- figure_3 & panel_tag_theme(10.8)

save_pdf(figure_3, staged_pdf_paths[[3]], width = 10.8, height = 9.2)
save_png(figure_3, staged_png_paths[[3]], width = 10.8, height = 9.2)


# Standalone map of state timing estimates ----------------------------------

choropleth_data <- state_boundaries %>%
  filter(.data$state_abbreviation != "DC") %>%
  left_join(
    state_timing %>%
      filter(.data$state != "District of Columbia") %>%
      select(
        "state",
        "age_label",
        "regret_minimizing_week",
        "regret_minimizing_week_axis"
      ),
    by = "state",
    relationship = "one-to-many"
  )

if (nrow(choropleth_data) != 50L * 6L ||
    anyNA(choropleth_data$regret_minimizing_week)) {
  stop("The choropleth requires one timing estimate per state and age group.")
}

state_timing_choropleth <- ggplot(choropleth_data) +
  geom_sf(
    aes(fill = .data$regret_minimizing_week_axis),
    color = "white",
    linewidth = 0.24
  ) +
  facet_wrap(vars(.data$age_label), ncol = 3) +
  scale_fill_gradientn(
    colors = timing_colors,
    limits = decision_timing_scale_limits,
    breaks = decision_timing_scale_breaks,
    labels = decision_timing_scale_weeks,
    oob = scales::squish,
    name = "MMWR week",
    guide = guide_colorbar(
      barwidth = grid::unit(40, "mm"),
      barheight = grid::unit(3, "mm"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  coord_sf(datum = NA, clip = "off") +
  theme_void(base_family = "Helvetica", base_size = 9) +
  theme(
    text = element_text(color = palette$ink),
    strip.text = element_text(
      face = "bold",
      color = palette$ink,
      size = 9,
      lineheight = 1.2,
      margin = margin(10, 0, -10, 0)
    ),
    strip.clip = "off",
    panel.spacing = grid::unit(5, "pt"),
    legend.position = "top",
    legend.justification = "right",
    legend.title = element_text(face = "plain"),
    legend.key.width = grid::unit(30, "pt"),
    legend.box.margin = margin(-2, 0, 0, 0),
    plot.margin = margin(5, 6, 5, 5)
  )

save_pdf(
  state_timing_choropleth,
  staged_pdf_paths[[4]],
  width = 11.0,
  height = 6.2
)
save_png(
  state_timing_choropleth,
  staged_png_paths[[4]],
  width = 11.0,
  height = 6.2
)

# Replace only the figures owned by this script after every output builds
# successfully. Supplemental publication figures share figure_dir, while the
# unreferenced standalone map is retained separately as a figure component.
owned_existing_files <- figure_paths[file.exists(figure_paths)]
if (length(owned_existing_files) > 0L) {
  unlink(owned_existing_files, force = TRUE)
}
if (!all(file.copy(staged_paths, figure_paths, overwrite = TRUE))) {
  stop("Could not publish one or more primary figure assets.")
}

if (!all(file.exists(figure_paths)) || any(file.info(figure_paths)$size <= 0L)) {
  stop("One or more primary figure assets are missing or empty.")
}

message("Wrote primary manuscript figures and retained figure component:")
message(paste(figure_paths, collapse = "\n"))
