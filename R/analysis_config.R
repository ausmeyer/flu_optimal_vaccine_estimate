parse_year_spec <- function(value, default = integer()) {
  if (is.null(value) || !nzchar(trimws(value))) {
    return(as.integer(default))
  }

  pieces <- unlist(strsplit(value, ",", fixed = TRUE), use.names = FALSE)
  years <- unlist(lapply(trimws(pieces), function(piece) {
    if (!nzchar(piece)) {
      return(integer())
    }
    if (grepl("^[0-9]{4}\\s*[:\\-]\\s*[0-9]{4}$", piece)) {
      bounds <- as.integer(unlist(strsplit(piece, "\\s*[:\\-]\\s*"), use.names = FALSE))
      return(seq.int(bounds[[1]], bounds[[2]]))
    }
    as.integer(piece)
  }), use.names = FALSE)

  years <- years[!is.na(years)]
  sort(unique(as.integer(years)))
}

collapse_year_ranges <- function(years) {
  years <- sort(unique(as.integer(years)))
  if (length(years) == 0L) {
    return("none")
  }

  starts <- years[c(TRUE, diff(years) > 1L)]
  ends <- c(years[which(diff(years) > 1L)], years[[length(years)]])
  paste0(
    ifelse(starts == ends, starts, paste0(starts, "-", ends)),
    collapse = "_"
  )
}

primary_season_start_years <- function() {
  parse_year_spec(
    Sys.getenv("PRIMARY_SEASON_START_YEARS", ""),
    default = c(2010:2018, 2023:2025)
  )
}

excluded_season_start_years <- function() {
  parse_year_spec(
    Sys.getenv("EXCLUDED_SEASON_START_YEARS", ""),
    default = 2019:2022
  )
}

season_set_slug <- function(years, prefix = "seasons") {
  paste0(prefix, "_", collapse_year_ranges(years))
}

file_fingerprint <- function(paths) {
  paths <- sort(unique(paths[nzchar(paths)]))
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Cannot fingerprint missing file(s): ", paste(missing, collapse = ", "))
  }
  if (length(paths) == 0L) {
    return("no_inputs")
  }

  hashes <- unname(tools::md5sum(paths))
  digest::digest(paste(paths, hashes, collapse = "|"), algo = "xxhash64")
}

latest_observed_week <- function(model_data) {
  required <- c("mmwr_year", "mmwr_week")
  missing <- setdiff(required, names(model_data))
  if (length(missing) > 0L) {
    stop("model_data missing: ", paste(missing, collapse = ", "))
  }

  model_data %>%
    dplyr::filter(!is.na(.data$mmwr_year), !is.na(.data$mmwr_week)) %>%
    dplyr::distinct(.data$mmwr_year, .data$mmwr_week) %>%
    dplyr::arrange(.data$mmwr_year, .data$mmwr_week) %>%
    dplyr::slice_tail(n = 1L)
}

model_path_with_seasons <- function(
    base_name,
    model_data,
    provenance_paths = character(),
    suffix = ".rds") {
  years <- sort(unique(as.integer(model_data$season_start_year)))
  observed <- latest_observed_week(model_data)
  observed_tag <- paste0(
    "through_",
    observed$mmwr_year[[1]],
    "w",
    sprintf("%02d", observed$mmwr_week[[1]])
  )
  provenance_tag <- if (length(provenance_paths) > 0L) {
    paste0("_src", substr(file_fingerprint(provenance_paths), 1L, 10L))
  } else {
    ""
  }
  file.path(
    "models",
    paste0(base_name, "_", season_set_slug(years), "_", observed_tag, provenance_tag, suffix)
  )
}

filter_to_primary_seasons <- function(data, label = deparse(substitute(data))) {
  if (!"season_start_year" %in% names(data)) {
    stop(label, " does not include season_start_year.")
  }

  years <- primary_season_start_years()
  out <- data %>%
    dplyr::filter(.data$season_start_year %in% years)

  if (nrow(out) == 0L) {
    stop(
      "No rows remain in ", label, " after filtering to PRIMARY_SEASON_START_YEARS=",
      paste(years, collapse = ",")
    )
  }

  message(
    "Using ", label, " seasons: ",
    paste(sort(unique(out$season_start_year)), collapse = ", ")
  )
  out
}

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, if (default) "true" else "false")
  tolower(value) %in% c("1", "true", "yes", "y")
}

env_number <- function(name, default, lower = -Inf, upper = Inf) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (length(value) != 1L || is.na(value) || value < lower || value > upper) {
    stop(name, " must be a number between ", lower, " and ", upper, ".")
  }
  value
}

timing_shift_sd <- function() {
  value <- Sys.getenv("TIMING_SHIFT_SD", "0.75")
  parsed <- suppressWarnings(as.numeric(value))
  if (length(parsed) != 1L || is.na(parsed) || parsed < 0) {
    stop("TIMING_SHIFT_SD must be a nonnegative number.")
  }
  parsed
}

ve_reference_weeks <- function() {
  env_number("VE_REFERENCE_WEEKS", default = 8, lower = 0)
}

read_state_population <- function(
    path = "data/processed/state_population.csv",
    expected_states = NULL) {
  if (!file.exists(path)) {
    stop(
      "State population data are missing. Run scripts/download_state_population.R before the analysis."
    )
  }

  population <- readr::read_csv(path, show_col_types = FALSE)
  required <- c("state", "population")
  missing <- setdiff(required, names(population))
  if (length(missing) > 0L) {
    stop(path, " missing: ", paste(missing, collapse = ", "))
  }
  if (anyNA(population$population) || any(population$population <= 0)) {
    stop(path, " contains missing or nonpositive population values.")
  }
  if (anyDuplicated(population$state)) {
    stop(path, " contains duplicate states.")
  }

  if (!is.null(expected_states)) {
    missing_states <- setdiff(sort(unique(expected_states)), population$state)
    if (length(missing_states) > 0L) {
      stop(path, " does not cover: ", paste(missing_states, collapse = ", "))
    }
  }
  population
}

write_analysis_manifest <- function(
    output_prefix,
    input_paths,
    configuration,
    model_path = NULL,
    output_dir = "outputs/primary/tables") {
  input_paths <- sort(unique(input_paths))
  missing <- input_paths[!file.exists(input_paths)]
  if (length(missing) > 0L) {
    stop("Manifest input file(s) missing: ", paste(missing, collapse = ", "))
  }

  inputs <- lapply(input_paths, function(path) {
    info <- file.info(path)
    list(
      path = path,
      md5 = unname(tools::md5sum(path)),
      bytes = unname(info$size),
      modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE)
    )
  })

  manifest <- list(
    analysis = output_prefix,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    r_version = R.version.string,
    model_path = model_path,
    inputs = inputs,
    configuration = configuration
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(output_dir, paste0(output_prefix, "_analysis_manifest.json"))
  jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(path)
}

analysis_output_dirs <- function(
    output_prefix,
    default_class = c("primary", "sensitivity")) {
  default_class <- match.arg(default_class)
  output_class <- Sys.getenv("OUTPUT_CLASS", default_class)
  if (!output_class %in% c("primary", "sensitivity")) {
    stop("OUTPUT_CLASS must be primary or sensitivity.")
  }

  base_dir <- if (output_class == "primary") {
    file.path("outputs", "primary")
  } else {
    file.path("outputs", "sensitivity", output_prefix)
  }
  paths <- list(
    class = output_class,
    base = base_dir,
    tables = file.path(base_dir, "tables"),
    figures = file.path(base_dir, "figures")
  )
  dir.create(paths$tables, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$figures, recursive = TRUE, showWarnings = FALSE)
  paths
}

gam_smoothing_slug <- function(
    global_k = Sys.getenv("GAM_GLOBAL_K", "14"),
    season_k = Sys.getenv("GAM_SEASON_K", "12"),
    gamma = Sys.getenv("GAM_GAMMA", "1")) {
  gamma_slug <- gsub("[^A-Za-z0-9]+", "p", as.character(gamma))
  paste0("_gk", global_k, "_sk", season_k, "_gamma", gamma_slug)
}
