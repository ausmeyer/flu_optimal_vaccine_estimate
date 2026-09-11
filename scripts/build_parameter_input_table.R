#!/usr/bin/env Rscript

# Summarize empirical inputs and key assumptions from the completed run.
# Fitting and computational settings remain in the Methods, manifests, and code.
# This builder does not regenerate draws or compile LaTeX.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})
source("R/outcome_time.R")

output_dir <- "outputs/primary/final_tables"
ve_path <- "data/processed/cdc_ve_estimates.csv"
population_path <- "data/processed/state_population.csv"

read_manifest <- function(prefix, primary = FALSE) {
  directory <- if (primary) "outputs/primary/tables" else file.path("outputs/sensitivity", prefix, "tables")
  path <- file.path(directory, paste0(prefix, "_analysis_manifest.json"))
  if (!file.exists(path)) stop("Complete the analysis before building the parameter table: ", path)
  jsonlite::read_json(path, simplifyVector = TRUE)
}

primary <- read_manifest("gam_primary", primary = TRUE)
config <- primary$configuration
age_prefixes <- c("gam_age_0_4", "gam_age_5_24", "gam_age_25_49", "gam_age_50_64", "gam_age_65_")
age_manifests <- lapply(age_prefixes, read_manifest, primary = TRUE)
sensitivity_prefixes <- c(
  "gam_no_timing_shift_primary", "gam_fixed_immune_lag_primary",
  "gam_slow_waning_primary", "gam_fast_waning_primary",
  "gam_spencer_ferdinands_fast_waning", "gam_spencer_ferdinands_onset_aligned",
  "gam_spencer_slow_waning", "gam_ve_reference_4w", "gam_ve_reference_12w"
)
sensitivities <- setNames(lapply(sensitivity_prefixes, read_manifest), sensitivity_prefixes)
sc <- function(prefix) sensitivities[[prefix]]$configuration
ve <- read_csv(ve_path, show_col_types = FALSE)
population <- read_csv(population_path, show_col_types = FALSE)
ve_hash <- unname(tools::md5sum(ve_path))
for (manifest in c(list(primary), age_manifests, sensitivities)) {
  input_hash <- manifest$inputs$md5[manifest$inputs$path == ve_path]
  if (length(input_hash) != 1L || input_hash != ve_hash) {
    stop("Parameter table cannot mix VE input versions: ", manifest$analysis)
  }
}
stopifnot(
  config$waning_prior == "ray_2019",
  config$immune_lag_prior == "discrete_about_2w",
  config$outcome_date_prior == "none",
  sc("gam_slow_waning_primary")$waning_prior == "ray_2019_lower95",
  sc("gam_fast_waning_primary")$waning_prior == "ray_2019_upper95",
  sc("gam_spencer_ferdinands_fast_waning")$waning_prior == "spencer_ferdinands_fast",
  sc("gam_spencer_slow_waning")$waning_prior == "spencer_slow",
  sc("gam_fixed_immune_lag_primary")$immune_lag_prior == "fixed_2w"
)
format_number <- function(x) format(x, trim = TRUE, scientific = FALSE)
age_pool_sizes <- ve %>% count(.data$ve_age_group)
all_age_size <- age_pool_sizes$n[age_pool_sizes$ve_age_group == "all"]
other_pool_sizes <- unique(age_pool_sizes$n[age_pool_sizes$ve_age_group != "all"])
population_year <- unique(population$source_year)
stopifnot(length(all_age_size) == 1L, length(other_pool_sizes) == 1L, length(population_year) == 1L)
outcome_distribution <- flu_ve_onset_to_enrollment_distribution()

rows <- list()
add_row <- function(section, id, parameter, specification, source,
                    implementation, citations = "", specification_latex = specification) {
  rows[[length(rows) + 1L]] <<- tibble(
    section = section, parameter_id = id, parameter = parameter,
    specification = specification, source = source, citation_keys = citations,
    implementation = implementation, specification_latex = specification_latex
  )
}

section <- "External inputs and published scenarios"
add_row(section, "ve_source_mixture", "Vaccine effectiveness (VE)",
  sprintf("%d all-age estimates and %d per age-specific pool, with 95%% confidence intervals (Table S1). Sample source seasons equally and independently of epidemic season; sample normal log odds ratios using each estimate and interval.", all_age_size, other_pool_sizes),
  "Cited in the VE input table; sampling choices are ours.",
  "data/processed/cdc_ve_estimates.csv; R/ve_draws_from_estimates.R:draw_reported_ve_from_estimates",
  specification_latex = sprintf("%d all-age estimates and %d per age-specific pool, with 95\\%% confidence intervals (Table~\\ref{tab:ve-inputs}). Sample source seasons equally and independently of epidemic season; sample normal log odds ratios using each estimate and interval.", all_age_size, other_pool_sizes))
add_row(section, "immune_lag", "Time to full immune response",
  "About 2 weeks. Use 1, 2, or 3 weeks with probabilities 0.15, 0.70, and 0.15; fixed 2 weeks in sensitivity analysis.",
  "CDC supports the 2-week center; distribution and sensitivity are our choices.",
  "R/ve_draws_from_estimates.R:draw_immune_response_lag", "cdcKeyFacts")
add_row(section, "ray_waning", "Primary waning",
  "Odds ratio 1.16 (95% CI, 1.13-1.20) per 28 days. beta ~ Normal(log(1.16),s^2), where s={log(1.20)-log(1.13)}/{2 Phi^(-1)(0.975)}. Fixed endpoint sensitivities use beta=log(1.13) or log(1.20).",
  "Ray et al.; normal reconstruction and fixed-endpoint sensitivities are our choices.",
  "R/ve_draws_from_estimates.R:draw_primary_waning", "ray2019",
  "Odds ratio 1.16 (95\\% CI, 1.13--1.20) per 28 days. $\\beta\\sim N\\{\\log(1.16),s^2\\}$, where\\newline $s=\\dfrac{\\log(1.20)-\\log(1.13)}{2\\Phi^{-1}(0.975)}$.\\newline Fixed endpoint sensitivities use $\\beta=\\log(1.13)$ or $\\log(1.20)$.")
add_row(section, "polynomial_waning", "Alternative waning scenarios",
  "For H weeks after full response, let b=H/2: fast r(H)=[(55-1.37b+0.18b^2-0.03b^3)/55]_0^1; slow r(H)=[(55-0.50b+0.05b^2-0.01b^3)/55]_0^1. The cubic coefficients use two-week units; division by 55 gives relative retention.",
  "Spencer et al. fast and slow scenarios; the fast curve draws on Ferdinands et al., whereas the slow curve is hypothetical.",
  "R/protection_models.R:spencer_relative_ve_retention; R/ve_draws_from_estimates.R:calibrate_initial_ve",
  "spencer2026,ferdinands2020",
  "For $H$ weeks after full response, let $b=H/2$:\\newline Fast $r(H)=\\left[\\dfrac{55-1.37b+0.18b^2-0.03b^3}{55}\\right]_0^1$.\\newline Slow $r(H)=\\left[\\dfrac{55-0.50b+0.05b^2-0.01b^3}{55}\\right]_0^1$.\\newline The cubic coefficients use two-week units; division by 55 gives relative retention.")
delay_categories <- paste0(
  outcome_distribution$lower_day, "-", outcome_distribution$upper_day, " days: ",
  format(outcome_distribution$n_observed, big.mark = ",", trim = TRUE),
  collapse = "; "
)
add_row(section, "outcome_date_delay", "Onset-to-enrollment delay",
  paste0(delay_categories, " (n=", format(sum(outcome_distribution$n_observed), big.mark = ",", trim = TRUE),
    "). For onset alignment, sample categories in proportion to these counts and an integer day uniformly within the selected category. Shift burden earlier by that delay; primary delay is zero."),
  "Balasubramani et al. counts summed across sites; within-category distribution and alignment are our choices.",
  "R/outcome_time.R:flu_ve_onset_to_enrollment_distribution;draw_outcome_date_offsets;apply_outcome_date_alignment",
  "balasubramani2020")
add_row(section, "geography_weights", "National population weights",
  sprintf("July 1, %d populations for %d jurisdictions (50 states and DC), held fixed across seasons and age groups.", population_year, nrow(population)),
  "U.S. Census Bureau; fixed-year weighting is our choice.",
  "scripts/download_state_population.R; R/latent_curve_gam.R:add_population_weighted_us",
  "censusStatePopulation2025")

section <- "Key author-specified assumptions"
add_row(section, "ve_calibration", "VE reference time",
  sprintf("%s weeks after full response; %s and %s weeks in sensitivity analyses. Estimate initial VE using VE at the reference time and the sampled waning curve.",
    format_number(config$ve_reference_weeks), format_number(sc("gam_ve_reference_4w")$ve_reference_weeks), format_number(sc("gam_ve_reference_12w")$ve_reference_weeks)),
  "Modeling assumption; no common reference time is reported by the source studies.",
  "R/ve_draws_from_estimates.R:calibrate_initial_ve")
add_row(section, "epidemic_shift", "Additional epidemic timing shift",
  sprintf("Normal(0, %s^2) weeks, shared across states; SD=%s in the no-shift sensitivity.",
    format_number(config$timing_shift_sd_weeks), format_number(sc("gam_no_timing_shift_primary")$timing_shift_sd_weeks)),
  "Modeling assumption.", "R/latent_curve_gam.R:sample_decision_scenarios;shift_burden_curves",
  specification_latex = sprintf("$N(0,%s^2)$ weeks, shared across states; SD=%s in the no-shift sensitivity.",
    format_number(config$timing_shift_sd_weeks), format_number(sc("gam_no_timing_shift_primary")$timing_shift_sd_weeks)))

parameters <- bind_rows(rows)
stopifnot(!anyDuplicated(parameters$parameter_id), !anyNA(parameters))
escape_tex <- function(text) {
  for (symbol in c("%", "&", "_", "#")) text <- gsub(symbol, paste0("\\", symbol), text, fixed = TRUE)
  text
}
table_rows <- unlist(lapply(unique(parameters$section), function(section) {
  block <- parameters[parameters$section == section, ]
  c(
    paste0("\\multicolumn{3}{l}{\\textbf{", escape_tex(section), "}} \\\\*"),
    vapply(seq_len(nrow(block)), function(i) {
      row <- block[i, ]
      specification <- if (row$specification_latex == row$specification) escape_tex(row$specification) else row$specification_latex
      source <- escape_tex(row$source)
      if (nzchar(row$citation_keys)) source <- paste0(source, " \\citep{", row$citation_keys, "}")
      paste0(escape_tex(row$parameter), " & ", specification, " & ", source, " \\\\[5pt]")
    }, character(1)),
    "\\addlinespace[3pt]"
  )
}), use.names = FALSE)

table_lines <- c(
  "% Generated by scripts/build_parameter_input_table.R; do not edit values by hand.",
  "\\clearpage", "\\begingroup", "\\singlespacing", "\\small",
  "\\setlength{\\tabcolsep}{4pt}", "\\setlength{\\LTcapwidth}{\\linewidth}",
  "\\renewcommand{\\arraystretch}{1.12}",
  "\\begin{longtable}{@{}>{\\raggedright\\arraybackslash}p{0.19\\linewidth}>{\\raggedright\\arraybackslash}p{0.50\\linewidth}>{\\raggedright\\arraybackslash}p{0.25\\linewidth}@{}}",
  "\\caption{Empirical inputs, published waning scenarios, and key parameter distributions. Table~\\ref{tab:ve-inputs} lists individual VE estimates and their sources.}\\label{tab:model-parameters}\\\\",
  "\\toprule", "Input & Value or distribution & Source and adaptation \\\\", "\\midrule", "\\endfirsthead",
  "\\multicolumn{3}{l}{Table~\\thetable. Continued} \\\\",
  "\\toprule", "Input & Value or distribution & Source and adaptation \\\\", "\\midrule", "\\endhead",
  "\\midrule", "\\multicolumn{3}{r}{Continued on next page} \\\\", "\\endfoot",
  "\\bottomrule", "\\endlastfoot", table_rows, "\\end{longtable}",
  "\\noindent $[x]_0^1=\\min\\{\\max(x,0),1\\}$; $\\Phi^{-1}$ is the standard-normal quantile. CDC indicates Centers for Disease Control and Prevention; CI, confidence interval; DC, District of Columbia; SD, standard deviation. Model fitting, protection equations, and other analysis settings are described in the Methods and reproducible code.",
  "\\endgroup"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(table_lines, file.path(output_dir, "table_s2_model_parameters.tex"))
write_csv(parameters, file.path(output_dir, "table_s2_model_parameters.csv"))
message("Generated concise parameter table with ", nrow(parameters), " rows from the completed publication configuration.")
