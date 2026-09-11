#!/usr/bin/env Rscript

# Install packages required to download inputs and run the publication analysis.
# Usage (from the repository root):
#   Rscript requirements.R

source("R/package_requirements.R")

pkgs <- publication_required_packages
minimum_versions <- publication_minimum_package_versions
minimum_r_version <- numeric_version("4.2.0")

if (getRversion() < minimum_r_version) {
  stop(
    "R 4.2.0 or later is required; found ", getRversion(), ".",
    call. = FALSE
  )
}
if (!isTRUE(capabilities("cairo"))) {
  stop(
    "This pipeline requires an R build with Cairo graphics support.",
    call. = FALSE
  )
}

installed <- vapply(
  pkgs,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)
outdated <- names(minimum_versions)[vapply(
  names(minimum_versions),
  function(package) {
    installed[[package]] &&
      utils::packageVersion(package) <
        base::package_version(minimum_versions[[package]])
  },
  logical(1)
)]
to_install <- union(pkgs[!installed], outdated)

if (length(to_install) == 0L) {
  message("All required packages are already installed.")
  quit(save = "no", status = 0)
}

message(
  "Installing missing or outdated packages: ",
  paste(to_install, collapse = ", ")
)
install.packages(to_install, repos = "https://cloud.r-project.org")

still_missing <- pkgs[!vapply(
  pkgs,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(still_missing) > 0L) {
  stop(
    "Failed to install: ", paste(still_missing, collapse = ", "),
    call. = FALSE
  )
}
still_outdated <- names(minimum_versions)[vapply(
  names(minimum_versions),
  function(package) {
    utils::packageVersion(package) <
      base::package_version(minimum_versions[[package]])
  },
  logical(1)
)]
if (length(still_outdated) > 0L) {
  stop(
    "Installed package version is too old: ",
    paste(
      paste0(
        still_outdated,
        " < ",
        minimum_versions[still_outdated]
      ),
      collapse = ", "
    ),
    call. = FALSE
  )
}

message("Package installation complete.")
