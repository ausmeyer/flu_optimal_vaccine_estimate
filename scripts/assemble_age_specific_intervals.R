#!/usr/bin/env Rscript

# Reassemble the combined interval table after the five separate age runs.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Expected at most one argument: the table directory.")
table_dir <- if (length(args)) args[[1]] else "outputs/primary/tables"
age_groups <- c("0-4", "5-24", "25-49", "50-64", "65+")
prefixes <- paste0("gam_age_", gsub("[^A-Za-z0-9]+", "_", age_groups))
paths <- file.path(table_dir, paste0(prefixes, "_optimal_week_intervals.csv"))
missing <- paths[!file.exists(paths)]
if (length(missing)) {
  stop("Missing age-specific interval inputs:\n", paste(missing, collapse = "\n"))
}

intervals <- bind_rows(Map(function(path, age_group) {
  rows <- read_csv(path, show_col_types = FALSE)
  if (!"age_group" %in% names(rows) || nrow(rows) == 0L ||
      anyNA(rows$age_group) || !identical(unique(rows$age_group), age_group)) {
    stop("Age label does not match ", age_group, ": ", path)
  }
  rows
}, paths, age_groups))

output_path <- file.path(table_dir, "gam_age_specific_optimal_week_intervals.csv")
write_csv(intervals, output_path)
message("Assembled age-specific intervals in canonical age order: ", output_path)
