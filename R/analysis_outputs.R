library(dplyr)

summarise_decision_outputs <- function(candidate_utilities) {
  optimal_weeks <- posterior_optimal_weeks(candidate_utilities)
  list(
    optimal_weeks = optimal_weeks,
    optimal_distribution = summarise_optimal_week_distribution(optimal_weeks),
    optimal_intervals = summarise_optimal_week_intervals(optimal_weeks),
    week_cutoff_probabilities = prob_optimal_by_week_cutoff(optimal_weeks),
    regret_curve = summarise_regret_curve(candidate_utilities),
    near_optimal_curve = summarise_near_optimal_curve(candidate_utilities)
  )
}

write_decision_outputs <- function(
    results,
    ve_draws,
    output_prefix,
    burden_draws = NULL,
    write_burden_draws = FALSE,
    output_dir = "outputs/primary/tables") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- function(suffix) file.path(output_dir, paste0(output_prefix, suffix))

  if (write_burden_draws && !is.null(burden_draws)) {
    readr::write_csv(burden_draws, path("_decision_burden_draws.csv.gz"))
  }
  readr::write_csv(results$optimal_weeks, path("_optimal_week_draws.csv"))
  readr::write_csv(results$optimal_distribution, path("_optimal_week_distribution.csv"))
  readr::write_csv(results$optimal_intervals, path("_optimal_week_intervals.csv"))
  readr::write_csv(
    results$week_cutoff_probabilities,
    path("_week_cutoff_probabilities.csv")
  )
  readr::write_csv(results$regret_curve, path("_regret_curve.csv"))
  readr::write_csv(results$near_optimal_curve, path("_near_optimal_curve.csv"))
  readr::write_csv(ve_draws, path("_ve_waning_draws.csv"))
}

save_decision_figures <- function(
    results,
    output_prefix,
    age_group,
    output_dir = "outputs/primary/figures",
    regret_y_label = "ILI visits potentially averted\nby optimizing timing",
    regret_caption = NULL) {
  write_figures <- tolower(Sys.getenv("WRITE_ANALYSIS_FIGURES", "true")) %in%
    c("1", "true", "yes", "y")
  if (!write_figures) {
    return(invisible(NULL))
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- function(suffix) file.path(output_dir, paste0(output_prefix, suffix))
  save_plot <- function(filename, plot, width, height) {
    ggplot2::ggsave(
      filename,
      plot,
      width = width,
      height = height,
      dpi = 300,
      bg = "white"
    )
  }

  save_plot(
    path("_us_optimal_week_distribution.png"),
    plot_us_optimal_week_distribution(results$optimal_distribution, age_group),
    7,
    4
  )
  save_plot(
    path("_us_regret_curve.png"),
    plot_regret_curve(
      results$regret_curve,
      state = "US",
      age_group = age_group,
      y_label = regret_y_label,
      caption = regret_caption
    ),
    7,
    4
  )
  save_plot(
    path("_us_relative_regret_curve.png"),
    plot_relative_regret_curve(
      results$near_optimal_curve,
      state = "US",
      age_group = age_group
    ),
    7,
    4
  )
  save_plot(
    path("_us_near_optimal_curve_5pct.png"),
    plot_near_optimal_curve(results$near_optimal_curve, state = "US", age_group = age_group),
    7,
    4
  )

  has_state_results <- any(results$optimal_distribution$state != "US")
  if (has_state_results) {
    save_plot(
      path("_state_optimal_week_distribution.png"),
      plot_state_faceted_optimal_week_distribution(
        results$optimal_distribution,
        age_group
      ),
      13,
      10
    )
    save_plot(
      path("_state_regret_curve.png"),
      plot_state_faceted_regret_curve(
        results$regret_curve,
        age_group,
        y_label = regret_y_label,
        caption = regret_caption
      ),
      13,
      10
    )
    save_plot(
      path("_state_relative_regret_curve.png"),
      plot_state_faceted_relative_regret_curve(
        results$near_optimal_curve,
        age_group
      ),
      13,
      10
    )
    save_plot(
      path("_state_near_optimal_curve_5pct.png"),
      plot_state_faceted_near_optimal_curve(
        results$near_optimal_curve,
        age_group
      ),
      13,
      10
    )
  }
}
