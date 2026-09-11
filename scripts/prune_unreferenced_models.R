#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

source("R/analysis_config.R")
source("R/publication_spec.R")

model_dir <- normalizePath("models", mustWork = TRUE)
manifest_paths <- publication_manifest_paths(c(2010:2018, 2023:2025))
missing_manifests <- manifest_paths[!file.exists(manifest_paths)]
if (length(missing_manifests) > 0L) {
  stop(
    "The publication analysis is incomplete; refusing to prune model cache. ",
    "Missing manifest(s):\n",
    paste(missing_manifests, collapse = "\n")
  )
}

referenced <- unique(unlist(lapply(manifest_paths, function(path) {
  manifest <- read_json(path, simplifyVector = TRUE)
  model_path <- manifest$model_path
  if (is.null(model_path) || length(model_path) == 0L || is.na(model_path) ||
      !nzchar(model_path)) {
    return(character())
  }
  if (!file.exists(model_path)) {
    stop("Manifest references a missing model: ", model_path)
  }
  normalizePath(model_path, mustWork = TRUE)
}), use.names = FALSE))

cached <- list.files(
  model_dir,
  pattern = "[.]rds$",
  full.names = TRUE
)
if (length(cached) == 0L) {
  message("No cached RDS models are present.")
  quit(save = "no", status = 0)
}
cached <- normalizePath(cached, mustWork = TRUE)
if (any(dirname(cached) != model_dir)) {
  stop("Model-cache pruning resolved a path outside the models directory.")
}

unreferenced <- setdiff(cached, referenced)
if (length(unreferenced) == 0L) {
  message("All cached models are referenced by completed analysis manifests.")
  quit(save = "no", status = 0)
}

removed <- file.remove(unreferenced)
if (!all(removed)) {
  stop(
    "Could not remove unreferenced model(s): ",
    paste(unreferenced[!removed], collapse = ", ")
  )
}
message(
  "Removed ", length(unreferenced),
  " unreferenced cached model(s):\n",
  paste(basename(unreferenced), collapse = "\n")
)
