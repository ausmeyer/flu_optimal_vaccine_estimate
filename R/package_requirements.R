# Direct R package requirements for the publication pipeline. The project
# records installed versions for provenance but intentionally does not freeze a
# complete package environment.

publication_required_packages <- c(
  "dplyr",
  "tidyr",
  "tibble",
  "purrr",
  "readr",
  "stringr",
  "ggplot2",
  "gtable",
  "scales",
  "mgcv",
  "MASS",
  "digest",
  "httr",
  "jsonlite",
  "MMWRweek",
  "rvest",
  "patchwork",
  "sf",
  "tigris"
)

publication_minimum_package_versions <- c(
  dplyr = "1.1.1",
  tidyr = "1.0.0",
  readr = "2.0.0",
  stringr = "1.3.0",
  ggplot2 = "4.0.0",
  scales = "1.1.0",
  jsonlite = "1.2",
  MMWRweek = "0.1.3",
  patchwork = "1.3.0",
  tigris = "1.4.0"
)
