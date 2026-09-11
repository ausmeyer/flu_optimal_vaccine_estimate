library(dplyr)
library(ggplot2)
library(scales)

source("R/calendar.R")

plot_week_axis <- function(start_week = 36L, end_week = 12L) {
  weeks <- c(seq.int(start_week, 52L), seq.int(1L, end_week))
  tibble::tibble(
    week = weeks,
    week_axis = seq_along(weeks)
  )
}

season_axis_breaks <- function(start_week = 36L, end_week = 12L) {
  axis <- plot_week_axis(start_week, end_week)
  center <- axis %>%
    filter(.data$week == 1L) %>%
    pull(.data$week_axis)
  wanted_axis <- sort(unique(c(
    seq.int(center, min(axis$week_axis), by = -5L),
    seq.int(center, max(axis$week_axis), by = 5L)
  )))
  axis %>%
    filter(.data$week_axis %in% wanted_axis) %>%
    arrange(.data$week_axis)
}

facet_panel_theme <- function() {
  theme(
    panel.border = element_rect(color = "#4D4D4D", fill = NA, linewidth = 0.35),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F2F2F2", color = "#4D4D4D", linewidth = 0.35),
    strip.text = element_text(face = "bold", color = "#222222"),
    panel.spacing = unit(0.7, "lines")
  )
}

regret_minimum_labels <- function(regret_curve, axis, group_vars = character(), x_offset = 0.35) {
  regret_curve %>%
    inner_join(axis, by = c("vaccination_week" = "week")) %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(max_regret_for_label = max(.data$upper_95, .data$mean_regret, na.rm = TRUE)) %>%
    slice_min(.data$mean_regret, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      week_axis_label = pmax(min(axis$week_axis), .data$week_axis - x_offset),
      regret_label_y = if_else(
        is.finite(.data$max_regret_for_label) & .data$max_regret_for_label > 0,
        0.92 * .data$max_regret_for_label,
        0
      ),
      regret_label = paste0("week ", .data$vaccination_week)
    )
}

standardized_visit_scale <- function(encounters_per_week) {
  scale_y_continuous(
    labels = label_number(
      scale = encounters_per_week,
      accuracy = 1,
      big.mark = ","
    )
  )
}

standardized_visit_caption <- function(encounters_per_week) {
  paste0(
    "Each simulation is compared with its optimal week.\n",
    "Standardized to ",
    label_comma()(encounters_per_week),
    " outpatient encounters per surveillance week."
  )
}

relative_regret_minimum_labels <- function(
    near_optimal_curve,
    axis,
    group_vars = character(),
    x_offset = 0.35) {
  near_optimal_curve %>%
    inner_join(axis, by = c("vaccination_week" = "week")) %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(max_regret_for_label = max(.data$mean_relative_regret, na.rm = TRUE)) %>%
    slice_min(.data$mean_relative_regret, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      week_axis_label = pmax(min(axis$week_axis), .data$week_axis - x_offset),
      regret_label_y = if_else(
        is.finite(.data$max_regret_for_label) & .data$max_regret_for_label > 0,
        0.92 * .data$max_regret_for_label,
        0
      ),
      regret_label = paste0("week ", .data$vaccination_week)
    )
}

add_relative_regret_intervals <- function(plot, plot_data) {
  interval_columns <- c(
    "lower_50_relative_regret",
    "upper_50_relative_regret",
    "lower_95_relative_regret",
    "upper_95_relative_regret"
  )
  if (!all(interval_columns %in% names(plot_data))) {
    return(plot)
  }

  plot +
    geom_ribbon(
      aes(
        ymin = .data$lower_95_relative_regret,
        ymax = .data$upper_95_relative_regret
      ),
      fill = "#9ECAE1",
      alpha = 0.35
    ) +
    geom_ribbon(
      aes(
        ymin = .data$lower_50_relative_regret,
        ymax = .data$upper_50_relative_regret
      ),
      fill = "#6BAED6",
      alpha = 0.45
    )
}

plot_us_optimal_week_distribution <- function(
    optimal_distribution,
    age_group = "all",
    start_week = 36L,
    end_week = 12L) {
  axis <- plot_week_axis(start_week, end_week)
  breaks <- season_axis_breaks(start_week, end_week)

  optimal_distribution %>%
    filter(.data$state == "US", .data$age_group == !!age_group) %>%
    inner_join(axis, by = c("optimal_week" = "week")) %>%
    ggplot(aes(x = .data$week_axis, y = .data$probability)) +
    geom_col(width = 0.85, fill = "#2C7FB8") +
    scale_x_continuous(
      breaks = breaks$week_axis,
      labels = breaks$week,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      x = "Optimal vaccination week (MMWR week)",
      y = "Posterior probability"
    ) +
    theme_minimal(base_size = 12)
}

plot_state_faceted_optimal_week_distribution <- function(
    optimal_distribution,
    age_group = "all",
    start_week = 36L,
    end_week = 12L) {
  axis <- plot_week_axis(start_week, end_week)
  breaks <- season_axis_breaks(start_week, end_week)

  optimal_distribution %>%
    filter(.data$state != "US", .data$age_group == !!age_group) %>%
    inner_join(axis, by = c("optimal_week" = "week")) %>%
    ggplot(aes(x = .data$week_axis, y = .data$probability)) +
    geom_col(width = 0.9, fill = "#2C7FB8") +
    facet_wrap(vars(.data$state), ncol = 7) +
    scale_x_continuous(
      breaks = breaks$week_axis,
      labels = breaks$week,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      x = "Optimal vaccination week (MMWR week)",
      y = "Posterior probability"
    ) +
    theme_minimal(base_size = 9) +
    facet_panel_theme()
}

plot_regret_curve <- function(
    regret_curve,
    state = "US",
    age_group = "all",
    start_week = 36L,
    end_week = 12L,
    standardized_encounters_per_week = 10000,
    y_label = "ILI visits potentially averted\nby optimizing timing",
    caption = NULL) {
  axis <- plot_week_axis(start_week, end_week)
  breaks <- season_axis_breaks(start_week, end_week)
  filtered_data <- regret_curve %>%
    filter(.data$state == !!state, .data$age_group == !!age_group) %>%
    arrange(.data$vaccination_week)
  plot_data <- filtered_data %>%
    inner_join(axis, by = c("vaccination_week" = "week"))
  minima <- regret_minimum_labels(filtered_data, axis, x_offset = 0.60)

  if (is.null(caption)) {
    caption <- standardized_visit_caption(standardized_encounters_per_week)
  }

  ggplot(plot_data, aes(x = .data$week_axis, y = .data$mean_regret)) +
    geom_ribbon(
      aes(ymin = .data$lower_95, ymax = .data$upper_95),
      fill = "#9ECAE1",
      alpha = 0.35
    ) +
    geom_ribbon(
      aes(ymin = .data$lower_50, ymax = .data$upper_50),
      fill = "#6BAED6",
      alpha = 0.45
    ) +
    geom_line(color = "#08519C", linewidth = 0.8) +
    geom_vline(
      data = minima,
      aes(xintercept = .data$week_axis),
      linetype = "dashed",
      color = "#4D4D4D",
      linewidth = 0.45
    ) +
    geom_text(
      data = minima,
      aes(
        x = .data$week_axis_label,
        y = .data$regret_label_y,
        label = .data$regret_label
      ),
      angle = 90,
      hjust = 1,
      vjust = 1,
      size = 3.8,
      color = "#333333",
      inherit.aes = FALSE
    ) +
    standardized_visit_scale(standardized_encounters_per_week) +
    scale_x_continuous(
      breaks = breaks$week_axis,
      labels = breaks$week,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      x = "Vaccination week (MMWR week)",
      y = y_label,
      caption = caption
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      plot.caption = element_text(color = "#4D4D4D", hjust = 0, size = 8.5)
    )
}

plot_state_faceted_regret_curve <- function(
    regret_curve,
    age_group = "all",
    start_week = 36L,
    end_week = 12L,
    free_y = TRUE,
    standardized_encounters_per_week = 10000,
    y_label = "ILI visits potentially averted\nby optimizing timing",
    caption = NULL) {
  axis <- plot_week_axis(start_week, end_week)
  breaks <- season_axis_breaks(start_week, end_week)
  filtered_data <- regret_curve %>%
    filter(.data$state != "US", .data$age_group == !!age_group) %>%
    arrange(.data$state, .data$vaccination_week)
  plot_data <- filtered_data %>%
    inner_join(axis, by = c("vaccination_week" = "week"))
  minima <- regret_minimum_labels(
    filtered_data,
    axis,
    group_vars = "state",
    x_offset = 1.55
  )

  if (is.null(caption)) {
    caption <- standardized_visit_caption(standardized_encounters_per_week)
  }

  ggplot(plot_data, aes(x = .data$week_axis, y = .data$mean_regret)) +
    geom_ribbon(
      aes(ymin = .data$lower_95, ymax = .data$upper_95),
      fill = "#9ECAE1",
      alpha = 0.30
    ) +
    geom_ribbon(
      aes(ymin = .data$lower_50, ymax = .data$upper_50),
      fill = "#6BAED6",
      alpha = 0.40
    ) +
    geom_line(color = "#08519C", linewidth = 0.45) +
    geom_vline(
      data = minima,
      aes(xintercept = .data$week_axis),
      linetype = "dashed",
      color = "#4D4D4D",
      linewidth = 0.30
    ) +
    geom_text(
      data = minima,
      aes(
        x = .data$week_axis_label,
        y = .data$regret_label_y,
        label = .data$regret_label
      ),
      angle = 90,
      hjust = 1,
      vjust = 1,
      size = 2.5,
      color = "#333333",
      inherit.aes = FALSE
    ) +
    facet_wrap(
      vars(.data$state),
      ncol = 7,
      scales = if (free_y) "free_y" else "fixed"
    ) +
    scale_x_continuous(
      breaks = breaks$week_axis,
      labels = breaks$week,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    standardized_visit_scale(standardized_encounters_per_week) +
    labs(
      x = "Vaccination week (MMWR week)",
      y = y_label,
      caption = caption
    ) +
    theme_minimal(base_size = 9) +
    facet_panel_theme() +
    theme(
      plot.caption = element_text(color = "#4D4D4D", hjust = 0, size = 7.5)
    )
}

plot_relative_regret_curve <- function(
    near_optimal_curve,
    state = "US",
    age_group = "all",
    relative_regret_tolerance = 0.05,
    start_week = 36L,
    end_week = 12L) {
  axis <- plot_week_axis(start_week, end_week)
  breaks <- season_axis_breaks(start_week, end_week)

  filtered_data <- near_optimal_curve %>%
    filter(
      .data$state == !!state,
      .data$age_group == !!age_group,
      abs(.data$relative_regret_tolerance - !!relative_regret_tolerance) < 1e-12
    )
  plot_data <- filtered_data %>%
    inner_join(axis, by = c("vaccination_week" = "week"))
  minima <- relative_regret_minimum_labels(filtered_data, axis, x_offset = 0.60)

  plot <- ggplot(plot_data, aes(x = .data$week_axis, y = .data$mean_relative_regret))
  plot <- add_relative_regret_intervals(plot, plot_data)
  plot +
    geom_line(color = "#08519C", linewidth = 0.8) +
    geom_vline(
      data = minima,
      aes(xintercept = .data$week_axis),
      linetype = "dashed",
      color = "#4D4D4D",
      linewidth = 0.45
    ) +
    geom_text(
      data = minima,
      aes(
        x = .data$week_axis_label,
        y = .data$regret_label_y,
        label = .data$regret_label
      ),
      angle = 90,
      hjust = 1,
      vjust = 1,
      size = 3.8,
      color = "#333333",
      inherit.aes = FALSE
    ) +
    scale_x_continuous(
      breaks = breaks$week_axis,
      labels = breaks$week,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(
      x = "Vaccination week (MMWR week)",
      y = "Preventable ILI burden lost\n(relative to optimal timing)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11)
    )
}

plot_state_faceted_relative_regret_curve <- function(
    near_optimal_curve,
    age_group = "all",
    relative_regret_tolerance = 0.05,
    start_week = 36L,
    end_week = 12L,
    free_y = TRUE) {
  axis <- plot_week_axis(start_week, end_week)
  breaks <- season_axis_breaks(start_week, end_week)

  filtered_data <- near_optimal_curve %>%
    filter(
      .data$state != "US",
      .data$age_group == !!age_group,
      abs(.data$relative_regret_tolerance - !!relative_regret_tolerance) < 1e-12
    )
  plot_data <- filtered_data %>%
    inner_join(axis, by = c("vaccination_week" = "week"))
  minima <- relative_regret_minimum_labels(
    filtered_data,
    axis,
    group_vars = "state",
    x_offset = 1.55
  )

  plot <- ggplot(plot_data, aes(x = .data$week_axis, y = .data$mean_relative_regret))
  plot <- add_relative_regret_intervals(plot, plot_data)
  plot +
    geom_line(color = "#08519C", linewidth = 0.45) +
    geom_vline(
      data = minima,
      aes(xintercept = .data$week_axis),
      linetype = "dashed",
      color = "#4D4D4D",
      linewidth = 0.30
    ) +
    geom_text(
      data = minima,
      aes(
        x = .data$week_axis_label,
        y = .data$regret_label_y,
        label = .data$regret_label
      ),
      angle = 90,
      hjust = 1,
      vjust = 1,
      size = 2.5,
      color = "#333333",
      inherit.aes = FALSE
    ) +
    facet_wrap(
      vars(.data$state),
      ncol = 7,
      scales = if (free_y) "free_y" else "fixed"
    ) +
    scale_x_continuous(
      breaks = breaks$week_axis,
      labels = breaks$week,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(
      x = "Vaccination week (MMWR week)",
      y = "Preventable ILI burden lost\n(relative to optimal timing)"
    ) +
    theme_minimal(base_size = 9) +
    facet_panel_theme()
}

plot_near_optimal_curve <- function(
    near_optimal_curve,
    state = "US",
    age_group = "all",
    relative_regret_tolerance = 0.05,
    start_week = 36L,
    end_week = 12L) {
  axis <- plot_week_axis(start_week, end_week)
  breaks <- season_axis_breaks(start_week, end_week)

  near_optimal_curve %>%
    filter(
      .data$state == !!state,
      .data$age_group == !!age_group,
      abs(.data$relative_regret_tolerance - !!relative_regret_tolerance) < 1e-12
    ) %>%
    inner_join(axis, by = c("vaccination_week" = "week")) %>%
    ggplot(aes(x = .data$week_axis, y = .data$probability_near_optimal)) +
    geom_col(width = 0.85, fill = "#2C7FB8") +
    scale_x_continuous(
      breaks = breaks$week_axis,
      labels = breaks$week,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
    labs(
      x = "Vaccination week (MMWR week)",
      y = paste0("Pr(within ", percent(relative_regret_tolerance, accuracy = 1), " of draw-best utility)")
    ) +
    theme_minimal(base_size = 12)
}

plot_state_faceted_near_optimal_curve <- function(
    near_optimal_curve,
    age_group = "all",
    relative_regret_tolerance = 0.05,
    start_week = 36L,
    end_week = 12L) {
  axis <- plot_week_axis(start_week, end_week)
  breaks <- season_axis_breaks(start_week, end_week)

  near_optimal_curve %>%
    filter(
      .data$state != "US",
      .data$age_group == !!age_group,
      abs(.data$relative_regret_tolerance - !!relative_regret_tolerance) < 1e-12
    ) %>%
    inner_join(axis, by = c("vaccination_week" = "week")) %>%
    ggplot(aes(x = .data$week_axis, y = .data$probability_near_optimal)) +
    geom_col(width = 0.9, fill = "#2C7FB8") +
    facet_wrap(vars(.data$state), ncol = 7) +
    scale_x_continuous(
      breaks = breaks$week_axis,
      labels = breaks$week,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
    labs(
      x = "Vaccination week (MMWR week)",
      y = paste0("Pr(within ", percent(relative_regret_tolerance, accuracy = 1), " of draw-best utility)")
    ) +
    theme_minimal(base_size = 9) +
    facet_panel_theme()
}
