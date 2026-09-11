# Local client for the public CDC FluView download service.
#
# These functions use the current official /flu2 endpoints directly and
# preserve the returned schemas expected by this project. Keeping the client
# here also makes endpoint, schema, and region validation part of this
# repository's reproducible input pipeline.

.fluview_region_map <- c(
  national = 3L,
  hhs = 1L,
  census = 2L,
  state = 5L
)

.fluview_metadata_cache <- new.env(parent = emptyenv())

# CDC download seasons begin at week 40. An analytic season beginning at
# week 36 also needs weeks 36-39 from the preceding download season.
fluview_download_seasons <- function(season_start_years) {
  sort(unique(c(season_start_years - 1L, season_start_years)))
}

fluview_raw_covers_burden_window <- function(
    data, season_start_years, first_available_date = NULL) {
  if (!all(c("year", "week") %in% names(data)) || nrow(data) == 0L) {
    return(FALSE)
  }
  observed <- MMWRweek::MMWRweek2Date(data$year, data$week)
  all(vapply(season_start_years, function(year) {
    expected <- seq(
      MMWRweek::MMWRweek2Date(year, 36L),
      MMWRweek::MMWRweek2Date(year + 1L, 22L),
      by = "week"
    )
    if (!is.null(first_available_date)) {
      expected <- expected[expected >= as.Date(first_available_date)]
    }
    length(expected) > 0L && all(expected %in% observed)
  }, logical(1)))
}

fluview_api_base <- function() {
  sub(
    "/+$",
    "",
    Sys.getenv("CDC_FLUVIEW_API_BASE", "https://gis.cdc.gov/flu2")
  )
}

fluview_timeout_seconds <- function() {
  value <- suppressWarnings(as.numeric(Sys.getenv("CDC_FLUVIEW_TIMEOUT_SECONDS", "120")))
  if (length(value) != 1L || is.na(value) || value <= 0) {
    stop("CDC_FLUVIEW_TIMEOUT_SECONDS must be a positive number.", call. = FALSE)
  }
  value
}

fluview_user_agent <- function() {
  paste0(
    "flu-optimal-vaccine-estimate/1.0 ",
    "(+https://github.com/ausmeyer/flu_optimal_vaccine_estimate)"
  )
}

fluview_clean_names <- function(data) {
  cleaned <- tolower(names(data))
  cleaned <- gsub("[[:punct:][:space:]]+", "_", cleaned)
  cleaned <- gsub("_+", "_", cleaned)
  cleaned <- gsub("(^_|_$)", "", cleaned)
  cleaned <- gsub("^x_", "", cleaned)
  names(data) <- make.unique(cleaned, sep = "_")
  data
}

fluview_to_numeric <- function(x) {
  x <- gsub("%", "", x, fixed = TRUE)
  x <- gsub(">", "", x, fixed = TRUE)
  x <- gsub("<", "", x, fixed = TRUE)
  x <- gsub(",", "", x, fixed = TRUE)
  x <- gsub(" ", "", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

fluview_response_detail <- function(response, body_path = NULL) {
  detail <- ""
  if (!is.null(body_path) && file.exists(body_path)) {
    size <- min(file.info(body_path)$size, 1000)
    if (is.finite(size) && size > 0) {
      detail <- tryCatch(
        readChar(body_path, nchars = size, useBytes = TRUE),
        error = function(e) ""
      )
    }
  } else {
    detail <- tryCatch(
      httr::content(response, as = "text", encoding = "UTF-8"),
      error = function(e) ""
    )
  }
  detail <- gsub("[[:space:]]+", " ", detail)
  substr(trimws(detail), 1L, 500L)
}

fluview_assert_success <- function(response, label, body_path = NULL) {
  status <- httr::status_code(response)
  if (status >= 200L && status < 300L) {
    return(invisible(response))
  }
  detail <- fluview_response_detail(response, body_path)
  stop(
    label, " failed with HTTP status ", status,
    if (nzchar(detail)) paste0(": ", detail) else ".",
    call. = FALSE
  )
}

fluview_get_metadata <- function(refresh = FALSE) {
  base <- fluview_api_base()
  if (!refresh && exists(base, envir = .fluview_metadata_cache, inherits = FALSE)) {
    return(get(base, envir = .fluview_metadata_cache, inherits = FALSE))
  }

  url <- paste0(base, "/GetPhase02InitApp?appVersion=Public")
  response <- httr::GET(
    url,
    httr::user_agent(fluview_user_agent()),
    httr::add_headers(
      Accept = "application/json, text/plain, */*",
      Referer = "https://gis.cdc.gov/fluview/"
    ),
    httr::timeout(fluview_timeout_seconds())
  )
  fluview_assert_success(response, "CDC FluView metadata request")

  body <- httr::content(response, as = "text", encoding = "UTF-8")
  metadata <- tryCatch(
    jsonlite::fromJSON(body),
    error = function(e) {
      stop("CDC FluView metadata was not valid JSON: ", conditionMessage(e), call. = FALSE)
    }
  )
  required <- c("seasonid", "description")
  if (!is.data.frame(metadata$seasons) ||
      !all(required %in% names(metadata$seasons)) ||
      nrow(metadata$seasons) == 0L) {
    stop("CDC FluView metadata did not contain the expected seasons schema.", call. = FALSE)
  }

  assign(base, metadata, envir = .fluview_metadata_cache)
  metadata
}

fluview_season_ids <- function(years, metadata) {
  available <- sort(unique(as.integer(metadata$seasons$seasonid)))
  available <- available[!is.na(available)]
  if (length(available) == 0L) {
    stop("CDC FluView metadata did not contain usable season IDs.", call. = FALSE)
  }
  if (is.null(years)) {
    return(available)
  }

  requested <- suppressWarnings(as.numeric(years))
  if (length(requested) == 0L || anyNA(requested)) {
    stop("FluView years must be numeric season start years or season IDs.", call. = FALSE)
  }
  requested_ids <- ifelse(requested > 1996, requested - 1960, requested)
  requested_ids <- sort(unique(as.integer(requested_ids)))
  unavailable <- setdiff(requested_ids, available)
  if (length(unavailable) > 0L) {
    message(
      "Ignoring unavailable CDC FluView season ID(s): ",
      paste(unavailable, collapse = ", ")
    )
  }
  selected <- intersect(requested_ids, available)
  if (length(selected) == 0L) {
    stop(
      "None of the requested seasons are available from CDC FluView. ",
      "Available season IDs: ", paste(available, collapse = ", "), ".",
      call. = FALSE
    )
  }
  selected
}

fluview_subregions <- function(region) {
  switch(
    region,
    national = list(list(ID = 0L, Name = "")),
    hhs = lapply(1:10, function(id) list(ID = id, Name = as.character(id))),
    census = lapply(1:9, function(id) list(ID = id, Name = as.character(id))),
    state = lapply(1:59, function(id) list(ID = id, Name = as.character(id)))
  )
}

fluview_download_parameters <- function(dataset, region, years, metadata) {
  list(
    AppVersion = "Public",
    DatasourceDT = list(list(ID = 1L, Name = dataset)),
    RegionTypeId = unname(.fluview_region_map[[region]]),
    SubRegionsDT = fluview_subregions(region),
    SeasonsDT = lapply(
      fluview_season_ids(years, metadata),
      function(id) list(ID = id, Name = as.character(id))
    )
  )
}

fluview_download_csvs <- function(parameters, label) {
  archive_path <- tempfile("fluview_", fileext = ".zip")
  extract_dir <- tempfile("fluview_extract_")
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(
    unlink(c(archive_path, extract_dir), recursive = TRUE, force = TRUE),
    add = TRUE
  )

  response <- httr::POST(
    paste0(fluview_api_base(), "/PostPhase02DataDownload"),
    httr::user_agent(fluview_user_agent()),
    httr::add_headers(
      Origin = "https://gis.cdc.gov",
      Accept = "application/json, text/plain, */*",
      Referer = "https://gis.cdc.gov/fluview/"
    ),
    encode = "json",
    body = parameters,
    httr::timeout(fluview_timeout_seconds()),
    httr::write_disk(archive_path, overwrite = TRUE)
  )
  fluview_assert_success(response, paste0("CDC FluView ", label, " download"), archive_path)

  archive_size <- file.info(archive_path)$size
  magic <- if (is.finite(archive_size) && archive_size >= 4L) {
    readBin(archive_path, what = "raw", n = 4L)
  } else {
    raw()
  }
  if (length(magic) < 2L ||
      !identical(as.integer(magic[1:2]), c(0x50L, 0x4bL))) {
    stop(
      "CDC FluView ", label,
      " response was not a nonempty ZIP archive. Response content type: ",
      httr::headers(response)[["content-type"]] %||% "unknown", ".",
      call. = FALSE
    )
  }

  listing <- tryCatch(
    utils::unzip(archive_path, list = TRUE),
    error = function(e) {
      stop("Could not inspect CDC FluView ", label, " ZIP: ", conditionMessage(e), call. = FALSE)
    }
  )
  member_names <- as.character(listing$Name)
  member_names <- member_names[grepl("[.]csv$", member_names, ignore.case = TRUE)]
  if (length(member_names) == 0L) {
    stop("CDC FluView ", label, " ZIP contained no CSV files.", call. = FALSE)
  }
  unsafe <- grepl("(^/|(^|/)[.][.](/|$))", member_names)
  if (any(unsafe) || anyDuplicated(basename(member_names))) {
    stop("CDC FluView ", label, " ZIP contained unsafe or duplicate CSV paths.", call. = FALSE)
  }

  paths <- tryCatch(
    utils::unzip(
      archive_path,
      files = member_names,
      exdir = extract_dir,
      junkpaths = TRUE,
      overwrite = TRUE
    ),
    error = function(e) {
      stop("Could not extract CDC FluView ", label, " ZIP: ", conditionMessage(e), call. = FALSE)
    }
  )
  if (length(paths) != length(member_names) || any(!file.exists(paths))) {
    stop("CDC FluView ", label, " ZIP extraction was incomplete.", call. = FALSE)
  }

  tables <- lapply(paths, function(path) {
    table <- tryCatch(
      utils::read.csv(path, skip = 1L, stringsAsFactors = FALSE),
      error = function(e) {
        stop("Could not parse ", basename(path), ": ", conditionMessage(e), call. = FALSE)
      }
    )
    if (ncol(table) == 0L) {
      stop("CDC FluView returned an empty-schema CSV: ", basename(path), ".", call. = FALSE)
    }
    fluview_clean_names(table)
  })
  names(tables) <- basename(paths)
  tables
}

fluview_assert_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(label, " is missing column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
}

fluview_validate_week_columns <- function(data, label) {
  if (anyNA(data$year) || anyNA(data$week) ||
      any(data$week < 1L | data$week > 53L)) {
    stop(label, " contains missing or invalid MMWR year/week identifiers.", call. = FALSE)
  }
}

fluview_download_ilinet <- function(
    region = c("national", "hhs", "census", "state"),
    years = NULL) {
  region <- match.arg(tolower(region), names(.fluview_region_map))
  metadata <- fluview_get_metadata()
  parameters <- fluview_download_parameters("ILINet", region, years, metadata)
  tables <- fluview_download_csvs(parameters, paste0("ILINet ", region))
  if (length(tables) != 1L) {
    stop(
      "CDC FluView ILINet download returned ", length(tables),
      " CSV files; exactly one was expected.",
      call. = FALSE
    )
  }
  data <- tables[[1]]

  required <- c(
    "region_type", "region", "year", "week", "weighted_ili",
    "unweighted_ili", "age_0_4", "age_25_49", "age_25_64",
    "age_5_24", "age_50_64", "age_65", "ilitotal",
    "num_of_providers", "total_patients"
  )
  fluview_assert_columns(data, required, "CDC FluView ILINet response")
  if (nrow(data) == 0L) {
    stop("CDC FluView ILINet response contained no rows.", call. = FALSE)
  }

  numeric_columns <- c(
    "weighted_ili", "unweighted_ili", "age_0_4", "age_25_49",
    "age_25_64", "age_5_24", "age_50_64", "age_65", "ilitotal",
    "num_of_providers", "total_patients"
  )
  data[numeric_columns] <- lapply(data[numeric_columns], fluview_to_numeric)
  data$year <- suppressWarnings(as.numeric(data$year))
  data$week <- suppressWarnings(as.numeric(data$week))
  fluview_validate_week_columns(data, "CDC FluView ILINet response")
  data$week_start <- MMWRweek::MMWRweek2Date(data$year, data$week)

  if (region == "national") {
    data$region <- "National"
  }
  if (region == "hhs") {
    data$region <- factor(data$region, levels = sprintf("Region %s", 1:10))
    if (!setequal(as.character(stats::na.omit(unique(data$region))), sprintf("Region %s", 1:10))) {
      stop("CDC FluView HHS ILINet response did not contain all 10 HHS regions.", call. = FALSE)
    }
  }
  if (region == "national" && !identical(unique(data$region), "National")) {
    stop("CDC FluView national ILINet response contained unexpected regions.", call. = FALSE)
  }
  if (region == "state" && anyNA(data$region)) {
    stop("CDC FluView state ILINet response contained a missing region name.", call. = FALSE)
  }

  class(data) <- c("tbl_df", "tbl", "data.frame")
  data <- suppressMessages(readr::type_convert(data))
  dplyr::arrange(data, .data$week_start)
}

fluview_mmwr_week_to_date <- function(year, week, day = NULL) {
  year <- as.numeric(year)
  week <- as.numeric(week)
  day <- if (is.null(day)) rep(1, length(week)) else as.numeric(day)
  week <- ifelse(week > 0 & week < 54, week, NA)
  as.Date(
    ifelse(is.na(week), NA, MMWRweek::MMWRweek2Date(year, week, day)),
    origin = "1970-01-01"
  )
}

fluview_download_who_nrevss <- function(
    region = c("national", "hhs", "census", "state"),
    years = NULL) {
  region <- match.arg(tolower(region), names(.fluview_region_map))
  metadata <- fluview_get_metadata()
  parameters <- fluview_download_parameters("WHO_NREVSS", region, years, metadata)
  tables <- fluview_download_csvs(parameters, paste0("WHO/NREVSS ", region))

  output_names <- sub(
    "who_nrevss_",
    "",
    tools::file_path_sans_ext(tolower(names(tables)))
  )
  if (anyDuplicated(output_names)) {
    stop("CDC FluView WHO/NREVSS response produced duplicate table names.", call. = FALSE)
  }
  names(tables) <- output_names

  tables <- lapply(names(tables), function(table_name) {
    data <- tables[[table_name]]
    data[data == "X"] <- NA
    data[data == "XX"] <- NA
    fluview_assert_columns(
      data,
      c("region_type", "region"),
      paste0("CDC FluView WHO/NREVSS table ", table_name)
    )
    if (all(c("year", "week") %in% names(data))) {
      data$year <- suppressWarnings(as.numeric(data$year))
      data$week <- suppressWarnings(as.numeric(data$week))
      fluview_validate_week_columns(
        data,
        paste0("CDC FluView WHO/NREVSS table ", table_name)
      )
      data$wk_date <- suppressWarnings(
        fluview_mmwr_week_to_date(data$year, data$week)
      )
    } else {
      data$wk_date <- as.Date(NA)
    }
    if (region == "national") {
      data$region <- "National"
    }
    class(data) <- c("tbl_df", "tbl", "data.frame")
    data
  })
  names(tables) <- output_names

  if (!any(grepl("(clinical_labs|combined_prior_to_2015_16)$", names(tables)))) {
    stop(
      "CDC FluView WHO/NREVSS response did not include a clinical-laboratory table.",
      call. = FALSE
    )
  }
  tables
}

# `%||%` is kept local to avoid adding rlang as a direct dependency.
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1]])) y else x
}
