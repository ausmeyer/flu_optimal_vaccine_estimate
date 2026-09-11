#!/usr/bin/env Rscript

# Chart contract: show the timing and between-season shape of the three national
# surveillance measures, not their incomparable absolute magnitudes. Static
# three-panel line figure for the supplement; each of 12 seasonal profiles sums
# to one, and a blue equal-season mixture accompanies neutral seasonal lines.
# MMWR-week coordinates retain week 53 between weeks 52 and 1. No posterior draws, timing
# perturbations, uncertainty bands, or additional curve fitting are introduced.

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

source("R/burden_posterior.R")

# Use the same calendar display coordinates as the other manuscript figures.
week_axis <- function(week) {
  dplyr::case_when(
    week >= 36L & week <= 52L ~ as.numeric(week - 35L),
    week == 53L ~ 17.5,
    week >= 1L & week <= 22L ~ as.numeric(week + 17L),
    TRUE ~ NA_real_
  )
}
axis_break_weeks <- c(36L, 40L, 44L, 48L, 52L, 4L, 8L, 12L, 16L, 20L)

figure_stem <- "figure_s3_national_surveillance_densities"
figure_dir <- "outputs/primary/final_figures"
source_dir <- "outputs/primary/figure_source_data"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

read_manifest <- function(prefix, primary = FALSE) {
  path <- if (primary) {
    file.path("outputs/primary/tables", paste0(prefix, "_analysis_manifest.json"))
  } else {
    file.path("outputs/sensitivity", prefix, "tables",
              paste0(prefix, "_analysis_manifest.json"))
  }
  if (!file.exists(path)) stop("Missing completed analysis manifest: ", path)
  manifest <- read_json(path, simplifyVector = TRUE)
  manifest$manifest_path <- path
  manifest
}

verify_input <- function(manifest, path) {
  index <- which(manifest$inputs$path == path)
  if (length(index) != 1L || !file.exists(path) ||
      unname(tools::md5sum(path)) != manifest$inputs$md5[[index]]) {
    stop("Figure input does not match completed analysis: ", path)
  }
  invisible(path)
}

read_models <- function(manifest) {
  if (is.null(manifest$model_path) || !file.exists(manifest$model_path)) {
    stop("Missing manifest-referenced GAM cache: ", manifest$model_path)
  }
  readRDS(manifest$model_path)
}

primary <- read_manifest("gam_primary", primary = TRUE)
ili_plus <- read_manifest("empirical_ili_plus_all_seasons_weighted")
laboratory <- read_manifest("gam_nrevss_national")
expected_years <- c(2010:2018, 2023:2025)
expected_seasons <- paste0(expected_years, "/", substr(expected_years + 1L, 3, 4))
for (manifest in list(primary, ili_plus, laboratory)) {
  if (!identical(as.integer(manifest$configuration$primary_season_start_years),
                 as.integer(expected_years))) {
    stop("Unexpected season coverage for ", manifest$analysis, ".")
  }
}
if (ili_plus$configuration$ili_weighting != "cdc_population_weighted" ||
    ili_plus$configuration$nrevss_lab_scope != "all_available" ||
    laboratory$configuration$target_virus != "all" ||
    laboratory$configuration$analysis_region != "national") {
  stop("The national curve figure requires the prespecified all-season measures.")
}

population_path <- "data/processed/state_population.csv"
state_path <- "data/processed/ilinet_state_all_age_primary_seasons.csv"
ili_plus_path <- "data/processed/national_ili_plus_primary_seasons.csv"
laboratory_path <- "data/processed/nrevss_clinical_labs_primary_seasons.csv"
verify_input(primary, population_path)
verify_input(primary, state_path)
verify_input(primary, "R/latent_curve_gam.R")
verify_input(ili_plus, ili_plus_path)
verify_input(laboratory, laboratory_path)

population <- read_csv(population_path, show_col_types = FALSE) %>%
  select("state", "population")
if (nrow(population) != 51L || anyDuplicated(population$state) ||
    any(!is.finite(population$population) | population$population <= 0)) {
  stop("National ILI curves require one positive population weight per jurisdiction.")
}

state_models <- read_models(primary)
if (!setequal(names(state_models), population$state)) {
  stop("Primary GAM caches and population weights cover different jurisdictions.")
}
state_fits <- imap_dfr(state_models, function(model, state_name) {
  model$data %>%
    transmute(
      state = state_name,
      season = as.character(.data$season),
      week = as.integer(.data$mmwr_week),
      max_week = as.integer(.data$max_week),
      central_burden = as.numeric(predict(model$fit, newdata = model$data,
                                          type = "response")),
      observed_burden = as.numeric(.data$ili_prop),
      observed_denominator = as.numeric(.data$total_patients_model)
    )
}) %>%
  left_join(population, by = "state", relationship = "many-to-one")
if (anyDuplicated(state_fits[c("state", "season", "week")])) {
  stop("Primary GAM grids contain duplicate jurisdiction-season-week rows.")
}
ili_curves <- state_fits %>%
  summarise(
    central_burden = weighted.mean(.data$central_burden, .data$population),
    observed_burden = weighted.mean(.data$observed_burden, .data$population,
                                    na.rm = TRUE),
    observed_denominator = sum(.data$observed_denominator),
    n_jurisdictions = n_distinct(.data$state),
    n_observed_jurisdictions = sum(is.finite(.data$observed_burden)),
    .by = c("season", "week", "max_week")
  ) %>%
  mutate(
    measure = "ILI",
    analysis = primary$analysis,
    source_path = state_path,
    source_md5 = unname(tools::md5sum(state_path)),
    model_path = primary$model_path,
    method = "Population-weighted mean of state GAM response predictions"
  )
if (any(ili_curves$n_jurisdictions != 51L)) {
  stop("A primary national curve is missing a jurisdiction on its prediction grid.")
}

laboratory_models <- read_models(laboratory)
if (length(laboratory_models) != 1L) stop("Expected one national laboratory GAM.")
laboratory_model <- laboratory_models[[1]]
laboratory_curves <- laboratory_model$data %>%
  transmute(
    measure = "Laboratory positivity",
    season = as.character(.data$season),
    week = as.integer(.data$mmwr_week),
    max_week = as.integer(.data$max_week),
    central_burden = as.numeric(predict(laboratory_model$fit,
                                        newdata = laboratory_model$data,
                                        type = "response")),
    observed_burden = as.numeric(.data$ili_prop),
    observed_denominator = as.numeric(.data$total_patients_model),
    nrevss_positive_tests = as.numeric(.data$positive_tests),
    nrevss_total_tests = as.numeric(.data$total_tests),
    nrevss_lab_series = as.character(.data$nrevss_lab_series),
    analysis = laboratory$analysis,
    source_path = laboratory_path,
    source_md5 = unname(tools::md5sum(laboratory_path)),
    model_path = laboratory$model_path,
    method = "National all-influenza laboratory-positivity GAM response prediction"
  )

ili_plus_observed <- read_csv(ili_plus_path, show_col_types = FALSE)
ili_plus_curves <- ili_plus_observed %>%
  transmute(
    draw = 1L, state = "US", age_group = "all",
    season = as.character(.data$season),
    week = as.integer(.data$mmwr_week),
    burden = .data$ili_plus_weighted,
    observed_burden = .data$ili_plus_weighted,
    ilinet_weighted_ili_proportion = .data$ilinet_weighted_ili_prob,
    nrevss_positive_proportion = .data$nrevss_positive_prob,
    ilinet_total_visits = .data$ilinet_total_visits,
    nrevss_total_tests = .data$nrevss_total_tests,
    nrevss_lab_series = .data$nrevss_lab_series,
    timing_shift_weeks = 0
  ) %>%
  transform_empirical_burden_draws(
    timing_shift_sd = 0, curve_smoothing = "cubic_spline",
    spline_spar = ili_plus$configuration$smoothing_spar
  ) %>%
  rename(central_burden = "burden") %>%
  select(-"draw", -"state", -"age_group", -"timing_shift_weeks") %>%
  mutate(
    measure = "ILI+",
    analysis = ili_plus$analysis,
    source_path = ili_plus_path,
    source_md5 = unname(tools::md5sum(ili_plus_path)),
    model_path = NA_character_,
    method = paste0("Production log1p cubic-spline transformation of observed weighted ILI+; spar=",
                    ili_plus$configuration$smoothing_spar)
  )

curves <- bind_rows(ili_curves, ili_plus_curves, laboratory_curves) %>%
  mutate(elapsed_week = week_to_season_index(
    .data$week, start_week = 36L, end_week = 22L, max_week = .data$max_week
  ) - 1L)
if (anyDuplicated(curves[c("measure", "season", "week")]) ||
    any(!is.finite(curves$central_burden) | curves$central_burden < 0)) {
  stop("Seasonal central curves contain duplicate, invalid, or negative values.")
}
coverage <- curves %>%
  summarise(n_seasons = n_distinct(.data$season), .by = "measure")
if (nrow(coverage) != 3L || any(coverage$n_seasons != 12L) ||
    !setequal(curves$season, expected_seasons)) {
  stop("Each surveillance measure must contain the 12 included seasons.")
}
curves <- curves %>%
  mutate(
    season_total_burden = sum(.data$central_burden),
    n_available_weeks = n(),
    normalized_weekly_mass = .data$central_burden / .data$season_total_burden,
    .by = c("measure", "season")
  )
if (any(!is.finite(curves$normalized_weekly_mass))) {
  stop("Cannot normalize a seasonal curve with zero or invalid total burden.")
}

# The representative curve is an equal mixture of the normalized available
# seasonal profiles, aligned by MMWR week. Missing early weeks and week 53
# in 52-week years have no mass in this mixture. Seasonal lines omit absent
# weeks. Week 53 remains a separate point in seasons that contain it; the
# mean line connects weeks 52 and 1 without plotting a partial-season mean
# at week 53. The complete mass accounting remains in the source data.
metadata <- curves %>%
  distinct(.data$measure, .data$season, .data$max_week,
           .data$season_total_burden, .data$n_available_weeks,
           .data$analysis, .data$source_path, .data$source_md5,
           .data$model_path, .data$method)
source_data <- tidyr::expand_grid(
  measure = unique(curves$measure), season = expected_seasons,
  week = c(36L:53L, 1L:22L)
) %>%
  left_join(metadata, by = c("measure", "season"), relationship = "many-to-one") %>%
  left_join(
    curves %>% select(-any_of(setdiff(names(metadata), c("measure", "season")))),
    by = c("measure", "season", "week"), relationship = "one-to-one"
  ) %>%
  mutate(
    source_week_present = !is.na(.data$central_burden),
    normalized_weekly_mass = coalesce(.data$normalized_weekly_mass, 0),
    curve_type = "season",
    n_seasons = 12L,
    n_contributing_seasons = as.integer(.data$source_week_present),
    plot_weekly_mass = if_else(.data$source_week_present,
                               .data$normalized_weekly_mass, NA_real_)
  )
season_sums <- source_data %>%
  summarise(total = sum(.data$normalized_weekly_mass),
            .by = c("measure", "season"))
if (any(abs(season_sums$total - 1) > 1e-12)) {
  stop("A normalized seasonal curve does not sum to one.")
}
mean_curves <- source_data %>%
  summarise(
    normalized_weekly_mass = mean(.data$normalized_weekly_mass),
    plot_weekly_mass = mean(.data$normalized_weekly_mass),
    n_contributing_seasons = sum(.data$source_week_present),
    n_seasons = n(),
    .by = c("measure", "week")
  ) %>%
  mutate(
    curve_type = "equal_season_mean", season = "Equal-season mean",
    plot_weekly_mass = if_else(.data$week == 53L, NA_real_, .data$plot_weekly_mass)
  )
mean_sums <- mean_curves %>%
  summarise(total = sum(.data$normalized_weekly_mass), .by = "measure")
if (any(abs(mean_sums$total - 1) > 1e-12) || any(mean_curves$n_seasons != 12L)) {
  stop("A representative curve does not preserve equal-season unit mass.")
}
source_data <- bind_rows(source_data, mean_curves) %>%
  mutate(plot_week = week_axis(.data$week)) %>%
  arrange(.data$measure, .data$curve_type, .data$season, .data$plot_week)

panel_labels <- c("ILI" = "ILI",
                  "Laboratory positivity" = "Laboratory positivity",
                  "ILI+" = "ILI+")
plot_data <- source_data %>%
  mutate(measure = factor(.data$measure, levels = names(panel_labels)))
y_axis_upper <- ceiling(max(plot_data$plot_weekly_mass, na.rm = TRUE) / 0.02) * 0.02
plot <- ggplot(plot_data, aes(x = .data$plot_week, y = .data$plot_weekly_mass)) +
  geom_line(data = filter(plot_data, .data$curve_type == "season",
                          .data$source_week_present),
            aes(group = .data$season, color = "Individual season"),
            linewidth = 0.35, alpha = 0.75) +
  geom_line(data = filter(plot_data, .data$curve_type == "equal_season_mean",
                          is.finite(.data$plot_weekly_mass)),
            aes(color = "Equal-season mean"), linewidth = 0.9) +
  facet_wrap(~measure, ncol = 1, labeller = as_labeller(panel_labels)) +
  scale_color_manual(values = c("Individual season" = "#AAAAAA",
                                "Equal-season mean" = "#2A5C85"),
                     breaks = c("Individual season", "Equal-season mean")) +
  scale_x_continuous(breaks = week_axis(axis_break_weeks), labels = axis_break_weeks,
                     limits = week_axis(c(36L, 22L)),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     breaks = seq(0, y_axis_upper, 0.02), limits = c(0, y_axis_upper),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(x = "MMWR week", y = "Share of seasonal burden", color = NULL) +
  theme_light(base_family = "Helvetica", base_size = 9) +
  theme(
    text = element_text(color = "#252525"),
    axis.text = element_text(color = "#252525"),
    axis.title = element_text(face = "plain"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E6E6E6", linewidth = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "#252525", hjust = 0),
    legend.position = "top", legend.justification = "right",
    legend.key.height = grid::unit(7, "pt"),
    legend.key.width = grid::unit(17, "pt"),
    plot.margin = margin(5, 6, 5, 5)
  )

# Separate panel letters from the measure names and align them in the outer
# margin. Their size is 12 points at the manuscript's 6.5-inch figure width.
plot <- ggplotGrob(plot)
strip_rows <- sort(unique(plot$layout$t[grepl("^strip-t", plot$layout$name)]))
if (length(strip_rows) != length(panel_labels)) stop("Unexpected Figure S3 panel layout.")
panel_tag_size <- 12 * 7.2 / 6.5
plot <- gtable::gtable_add_cols(plot, grid::unit(panel_tag_size + 7, "pt"), pos = 0)
for (i in seq_along(strip_rows)) {
  plot <- gtable::gtable_add_grob(
    plot,
    grid::textGrob(LETTERS[[i]], x = grid::unit(5, "pt"), y = 0.5,
                   just = c("left", "centre"),
                   gp = grid::gpar(fontfamily = "Helvetica", fontface = "bold",
                                   fontsize = panel_tag_size, col = "#252525")),
    t = strip_rows[[i]], l = 1, clip = "off", name = paste0("panel-tag-", i)
  )
}

staging_dir <- tempfile("national-curve-figure-")
dir.create(staging_dir)
on.exit(unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)
source_name <- paste0(figure_stem, "_source_data.csv")
write_csv(source_data, file.path(staging_dir, source_name))
for (extension in c("pdf", "png")) {
  ggsave(file.path(staging_dir, paste0(figure_stem, ".", extension)),
         plot = plot, width = 7.2, height = 6.4, units = "in", bg = "white",
         device = if (extension == "pdf") grDevices::cairo_pdf else "png", dpi = 300)
}
if (!all(file.copy(file.path(staging_dir, paste0(figure_stem, c(".pdf", ".png"))),
                   figure_dir, overwrite = TRUE)) ||
    !file.copy(file.path(staging_dir, source_name), source_dir, overwrite = TRUE)) {
  stop("Could not publish all national surveillance figure artifacts.")
}
message("Built national surveillance profiles from completed production inputs.")
