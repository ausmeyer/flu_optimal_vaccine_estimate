#!/usr/bin/env Rscript

# Supplementary protection distributions from the saved production draws.
# Chart contract: six age panels; daily probability histograms (0.5-percentage-
# point bins) and median curves. Linear blue intensity shows unconditional
# draw probability from 0% to 2%, with values above 2% sharing the darkest shade.
# All draws, including zero protection, remain in the histograms and medians.
# Daily evaluation is for display only; decision calculations remain weekly.
# Static PDF/PNG plus auditable source tables.
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(readr)
  library(tidyr)
})
source("R/decision_engine.R")
source("R/publication_spec.R")

prefixes <- publication_primary_prefixes()
age_labels <- c("All ages", "Ages 0–4", "Ages 5–24", "Ages 25–49",
                "Ages 50–64", "Ages 65+")
elapsed_days <- 0:(39L * 7L)
elapsed_weeks <- elapsed_days / 7
bin_width <- 0.005
n_bins <- as.integer(round(1 / bin_width))
density_rows <- list()
summary_rows <- list()

for (i in seq_along(prefixes)) {
  prefix <- prefixes[[i]]
  draw_path <- file.path("outputs/primary/tables", paste0(prefix, "_ve_waning_draws.csv"))
  manifest_path <- file.path("outputs/primary/tables", paste0(prefix, "_analysis_manifest.json"))
  manifest <- read_json(manifest_path, simplifyVector = TRUE)
  draws <- read_csv(draw_path, show_col_types = FALSE)
  stopifnot(nrow(draws) == manifest$configuration$n_draws,
            !anyDuplicated(draws$draw), n_distinct(draws$age_group) == 1L,
            all(draws$protection_model == "exponential_effect_ratio"),
            all(draws$ve_reference_weeks == manifest$configuration$ve_reference_weeks))

  # Evaluate the exact bounded exponential effect-ratio transformation used
  # in the decision engine. Keep the saved lag and back-calibrated VE together.
  since_response <- outer(draws$immune_lag_weeks, elapsed_weeks, function(lag, t) t - lag)
  protection <- pmin(pmax(
    1 - (1 - draws$initial_ve) *
      exp(draws$beta_wane_per_28d * pmax(since_response, 0) / 4), 0), 1)
  protection[since_response < 0] <- 0
  stopifnot(all(is.finite(protection)), all(protection >= 0 & protection <= 1))

  # Check the plotted transformation against the production utility engine:
  # unit burden at one week makes utility equal protection at that week.
  check_ids <- unique(c(seq_len(min(20L, nrow(draws))),
                        which(draws$initial_ve_bounded)[seq_len(min(5L, sum(draws$initial_ve_bounded)))]))
  check_ids <- check_ids[!is.na(check_ids)]
  check_burden <- crossing(draw = draws$draw[check_ids], elapsed = 0:39) %>%
    mutate(state = sprintf("t%02d", .data$elapsed), age_group = draws$age_group[[1]],
           week = if_else(.data$elapsed <= 17L, 36L + .data$elapsed, .data$elapsed - 17L),
           max_week = 53L, burden = 1)
  checked <- evaluate_candidate_utilities(check_burden, draws[check_ids, ], candidate_weeks = 36L) %>%
    mutate(elapsed = as.integer(sub("t", "", .data$state)),
           expected = protection[cbind(match(.data$draw, draws$draw), .data$elapsed * 7L + 1L)])
  stopifnot(nrow(checked) == length(check_ids) * 40L,
            max(abs(checked$utility - checked$expected)) < 1e-12)

  for (j in seq_along(elapsed_weeks)) {
    values <- protection[, j]
    bin <- pmin(floor(values / bin_width) + 1L, n_bins)
    counts <- tabulate(bin, nbins = n_bins)
    quantiles <- quantile(values, c(.025, .25, .5, .75, .975), names = FALSE)
    density_rows[[length(density_rows) + 1L]] <- tibble(
      analysis = prefix, age_group = draws$age_group[[1]], age_label = age_labels[[i]],
      days_since_vaccination = elapsed_days[[j]],
      weeks_since_vaccination = elapsed_weeks[[j]],
      protection_bin_lower = (0:(n_bins - 1L)) * bin_width,
      protection_bin_upper = (1:n_bins) * bin_width,
      n_draws_in_bin = counts, n_draws = nrow(draws),
      probability = counts / nrow(draws), draw_source = draw_path)
    summary_rows[[length(summary_rows) + 1L]] <- tibble(
      analysis = prefix, age_group = draws$age_group[[1]], age_label = age_labels[[i]],
      days_since_vaccination = elapsed_days[[j]],
      weeks_since_vaccination = elapsed_weeks[[j]], n_draws = nrow(draws),
      mean_protection = mean(values), lower_95 = quantiles[[1]], lower_50 = quantiles[[2]],
      median_protection = quantiles[[3]], upper_50 = quantiles[[4]], upper_95 = quantiles[[5]],
      probability_zero = mean(values == 0), probability_one = mean(values == 1),
      draw_source = draw_path, draw_source_md5 = unname(tools::md5sum(draw_path)))
  }
}

density_data <- bind_rows(density_rows)
summary_data <- bind_rows(summary_rows)
mass <- density_data %>% summarise(mass = sum(.data$probability),
                                  .by = c("analysis", "weeks_since_vaccination"))
stopifnot(max(abs(mass$mass - 1)) < 1e-12)

plot <- ggplot(density_data, aes(x = .data$days_since_vaccination)) +
  geom_raster(aes(y = 100 * (.data$protection_bin_lower + .data$protection_bin_upper) / 2,
                  fill = if_else(.data$probability > 0, 100 * .data$probability, NA_real_)),
              interpolate = FALSE) +
  geom_line(data = summary_data, aes(y = 100 * .data$median_protection),
            color = "#252525", linewidth = 0.4, linetype = "longdash") +
  facet_wrap(~factor(age_label, levels = age_labels), ncol = 2) +
  scale_fill_gradientn(colours = c("#EEF5FA", "#A9C7DB", "#568BB1", "#2A5C85", "#102A43"),
                       limits = c(0, 2), oob = scales::squish,
                       breaks = c(0, 0.5, 1, 1.5, 2),
                       labels = c("0", "0.5", "1", "1.5", "\u22652"),
                       na.value = "white", name = "Draws per 0.5-point bin (%)") +
  scale_x_continuous(breaks = seq(0, 36, 6) * 7, labels = seq(0, 36, 6),
                     limits = c(-.5, max(elapsed_days) + .5), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100),
                     expand = expansion(add = c(0.6, 0.6))) +
  labs(x = "Weeks since vaccination", y = "Modeled vaccine protection (%)") +
  theme_light(base_family = "Helvetica", base_size = 9) +
  theme(text = element_text(color = "#252525"),
        axis.text = element_text(color = "#252525"), panel.grid = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = "bold", color = "#252525"),
        legend.position = "bottom", legend.title.position = "top",
        legend.key.width = grid::unit(25, "pt"),
        plot.margin = margin(6, 8, 6, 6))

figure_dir <- "outputs/primary/final_figures"
source_dir <- "outputs/primary/figure_source_data"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
stem <- "figure_s4_vaccine_protection_density"
write_csv(density_data, file.path(source_dir, paste0(stem, "_source_data.csv")))
write_csv(summary_data, file.path(source_dir, paste0(stem, "_summary.csv")))
ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = 7, height = 7,
       device = grDevices::cairo_pdf, bg = "white")
ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = 7, height = 7,
       dpi = 300, bg = "white")
message("Generated protection-density figure from saved draws; production-equation checks passed.")
