#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("R/analysis_config.R")

n_draws <- Sys.getenv("N_DRAWS", "5000")
base_years <- primary_season_start_years()
if (length(base_years) < 3L) {
  stop("Leave-one-season-out sensitivity requires at least 3 seasons.")
}

run_one <- function(omitted_year) {
  years <- setdiff(base_years, omitted_year)
  prefix <- paste0("gam_loso_drop_", omitted_year)
  manifest <- file.path(
    "outputs", "sensitivity", prefix, "tables",
    paste0(prefix, "_analysis_manifest.json")
  )
  if (env_flag("SKIP_EXISTING", FALSE) && file.exists(manifest)) {
    message("Skipping leave-one-season-out sensitivity for ", omitted_year)
    return(invisible(NULL))
  }
  message("Running leave-one-season-out sensitivity without ", omitted_year, " -> ", prefix)
  status <- system2(
    command = "Rscript",
    args = "scripts/run_gam_primary_analysis.R",
    env = c(
      paste0("N_DRAWS=", n_draws),
      paste0("PRIMARY_SEASON_START_YEARS=", paste(years, collapse = ",")),
      paste0("OUTPUT_PREFIX=", prefix),
      paste0("OUTPUT_CLASS=", Sys.getenv("OUTPUT_CLASS", "sensitivity")),
      "REUSE_GAM_MODELS=true"
    )
  )
  if (!identical(status, 0L)) {
    stop("Leave-one-season-out run failed for omitted season_start_year=", omitted_year)
  }
}

invisible(lapply(base_years, run_one))
message("GAM leave-one-season-out sensitivity complete.")
