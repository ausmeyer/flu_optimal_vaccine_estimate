library(dplyr)
library(purrr)
library(readr)
library(rvest)
library(stringr)
library(tibble)

# Download archived CDC season pages to a raw table collection. The normalization
# script validates the expected coverage and combines the scraped 2011/12-2018/19
# estimates with cited 2010/11 and 2023/24-2025/26 published estimates encoded
# there. Recent estimates have explicit source and final/preliminary status.

season_pages <- tribble(
  ~season, ~url,
  "2010/11", "https://archive.cdc.gov/www_cdc_gov/flu/vaccines-work/2010-2011.html",
  "2011/12", "https://archive.cdc.gov/www_cdc_gov/flu/vaccines-work/2011-2012.html",
  "2012/13", "https://archive.cdc.gov/www_cdc_gov/flu/vaccines-work/2012-2013.html",
  "2013/14", "https://archive.cdc.gov/www_cdc_gov/flu/vaccines-work/2013-2014.html",
  "2014/15", "https://archive.cdc.gov/www_cdc_gov/flu/vaccines-work/2014-2015.html",
  "2015/16", "https://archive.cdc.gov/www_cdc_gov/flu/vaccines-work/2015-2016.html",
  "2016/17", "https://archive.cdc.gov/www_cdc_gov/flu/vaccines-work/2016-2017.html",
  "2017/18", "https://archive.cdc.gov/www_cdc_gov/flu/vaccines-work/2017-2018.html",
  "2018/19", "https://archive.cdc.gov/www_cdc_gov/flu/vaccines-work/2018-2019.html"
)

clean_ci <- function(x) {
  out <- str_match(x, "(-?\\d+)\\s*(?:,|to)\\s*(-?\\d+)")
  tibble(lower_95 = as.numeric(out[, 2]), upper_95 = as.numeric(out[, 3]))
}

extract_page_tables <- function(season, url) {
  page <- tryCatch(
    read_html(url),
    error = function(e) {
      warning("Could not read ", season, " page: ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(page)) {
    return(tibble(
      season = season,
      source_url = url,
      table_index = NA_integer_,
      raw_table = list(tibble()),
      scrape_status = "failed"
    ))
  }
  tables <- html_table(page, fill = TRUE)

  imap_dfr(tables, function(tbl, table_index) {
    names(tbl) <- make.names(names(tbl), unique = TRUE)
    tbl <- as_tibble(tbl)

    if (!any(str_detect(names(tbl), regex("VE|effectiveness", ignore_case = TRUE)))) {
      return(tibble())
    }

    tibble(
      season = season,
      source_url = url,
      table_index = table_index,
      raw_table = list(tbl),
      scrape_status = "ok"
    )
  })
}

raw_tables <- pmap_dfr(season_pages, extract_page_tables)

dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
saveRDS(raw_tables, "data/raw/cdc_ve_raw_tables.rds")

message("Saved raw CDC VE tables to data/raw/cdc_ve_raw_tables.rds")
message("Next step: normalize and validate the publication VE input.")
