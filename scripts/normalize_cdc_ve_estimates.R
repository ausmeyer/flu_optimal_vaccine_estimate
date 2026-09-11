library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

raw_tables <- readRDS("data/raw/cdc_ve_raw_tables.rds")

parse_percent <- function(x) {
  as.numeric(str_remove_all(as.character(x), "[^0-9\\-\\.]+")) / 100
}

parse_ci <- function(x) {
  mat <- str_match(as.character(x), "\\(?\\s*(-?\\d+)\\s*(?:,|to|–|-)\\s*(-?\\d+)\\s*\\)?")
  tibble(
    lower_95 = as.numeric(mat[, 2]) / 100,
    upper_95 = as.numeric(mat[, 3]) / 100
  )
}

standardize_age <- function(age) {
  age <- str_squish(as.character(age))
  case_when(
    str_detect(age, regex("^All", ignore_case = TRUE)) ~ "all",
    str_detect(age, "6\\s*mos|6\\s*months") & str_detect(age, "8") ~ "6m-8",
    str_detect(age, "^9") & str_detect(age, "17") ~ "9-17",
    str_detect(age, "^18") & str_detect(age, "49") ~ "18-49",
    str_detect(age, "^50") & str_detect(age, "64") ~ "50-64",
    str_detect(age, "65") ~ "65+",
    TRUE ~ age
  )
}

extract_first_age_table <- function(season) {
  season_tables <- raw_tables %>%
    filter(.data$season == !!season, .data$scrape_status == "ok")

  if (nrow(season_tables) == 0) {
    return(tibble())
  }

  tbl <- season_tables$raw_table[[1]]
  names(tbl) <- make.names(names(tbl), unique = TRUE)

  age_col <- names(tbl)[str_detect(names(tbl), regex("Age.group", ignore_case = TRUE))][1]
  ve_col <- names(tbl)[str_detect(names(tbl), regex("Adjusted.VE|^VE$", ignore_case = TRUE))][1]
  ci_col <- names(tbl)[str_detect(names(tbl), regex("Adjusted.95|X95|95", ignore_case = TRUE))][1]

  if (any(is.na(c(age_col, ve_col, ci_col)))) {
    warning("Could not identify VE columns for ", season)
    return(tibble())
  }

  ci <- parse_ci(tbl[[ci_col]])

  tibble(
    season = season,
    ve_age_group = standardize_age(tbl[[age_col]]),
    ve = parse_percent(tbl[[ve_col]]),
    lower_95 = ci$lower_95,
    upper_95 = ci$upper_95,
    outcome = "medically_attended_influenza_A_or_B",
    source_id = "cdc_archived_season_page",
    citation_key = paste0("cdcVE", str_sub(season, 1, 4)),
    source_status = "scraped",
    source_url = season_tables$source_url[[1]],
    source_age_group = str_squish(as.character(tbl[[age_col]])),
    source_table = paste0(
      "First age-specific VE table (HTML table ",
      season_tables$table_index[[1]], ")"
    ),
    estimate_status = "final"
  ) %>%
    filter(!is.na(.data$ve), !is.na(.data$lower_95), !is.na(.data$upper_95))
}

scraped <- map_dfr(unique(raw_tables$season), extract_first_age_table)

# The archived 2015/16 CDC table reports an overall lower confidence limit
# of 43%. Jackson et al. report 48% (95% CI, 41%-55%) in the primary paper's
# Abstract and Results. Retain the raw scrape and explicitly prefer the
# primary source for this one confidence interval.
manual_2015_overall <- scraped %>%
  filter(.data$season == "2015/16", .data$ve_age_group == "all") %>%
  mutate(
    ve = 0.48,
    lower_95 = 0.41,
    upper_95 = 0.55,
    source_id = "jackson_2017_primary_overall",
    citation_key = "jackson2017VE",
    source_status = "manual_source_correction",
    source_url = "https://pmc.ncbi.nlm.nih.gov/articles/PMC5727917/",
    source_age_group = "All ages (6 months and older)",
    source_table = "Abstract and Results, overall vaccine effectiveness",
    estimate_status = "final"
  )

manual_2010 <- tribble(
  ~season, ~ve_age_group, ~ve, ~lower_95, ~upper_95, ~outcome, ~source_id, ~citation_key, ~source_status,
  "2010/11", "all",   0.60,  0.54, 0.66, "medically_attended_influenza_A_or_B", "treanor_2012_table3", "treanor2012", "manual",
  "2010/11", "6m-8",  0.63,  0.52, 0.72, "medically_attended_influenza_A_or_B", "treanor_2012_text_combined_6m_8", "treanor2012", "manual",
  "2010/11", "9-17",  0.51,  0.36, 0.62, "medically_attended_influenza_A_or_B", "treanor_2012_table3_9_49_proxy", "treanor2012", "manual_proxy",
  "2010/11", "18-49", 0.51,  0.36, 0.62, "medically_attended_influenza_A_or_B", "treanor_2012_table3_9_49_proxy", "treanor2012", "manual_proxy",
  "2010/11", "50-64", 0.51,  0.25, 0.68, "medically_attended_influenza_A_or_B", "treanor_2012_table3", "treanor2012", "manual",
  "2010/11", "65+",   0.36, -0.22, 0.66, "medically_attended_influenza_A_or_B", "treanor_2012_table3", "treanor2012", "manual"
) %>%
  mutate(
    source_url = "https://pmc.ncbi.nlm.nih.gov/articles/PMC3657521/",
    source_age_group = case_when(
      .data$ve_age_group == "all" ~ "All ages",
      .data$ve_age_group == "6m-8" ~ "6 months-8 years",
      .data$ve_age_group %in% c("9-17", "18-49") ~ "9-49 years",
      .data$ve_age_group == "50-64" ~ "50-64 years",
      .data$ve_age_group == "65+" ~ "65 years and older"
    ),
    source_table = if_else(
      .data$ve_age_group == "6m-8",
      "Results text, combined 6 months-8 years",
      "Table 3, all influenza, adjusted VE"
    ),
    estimate_status = "final"
  )

# Published U.S. Flu VE Network estimates, transcribed from the source
# locations recorded below. Use one all-influenza estimate per season and VE age pool. The
# 2023/24 youngest stratum starts at 8 months and is an explicit proxy for
# the historical 6m-8 pool. The 2024/25 final estimate replaces the interim
# estimate without expanding the existing age-specific season pools.
# The 2025/26 interim estimate also enters only the all-age pool. VISION estimates are not
# substituted for the U.S. Flu VE Network estimates.
manual_recent <- tribble(
  ~season, ~ve_age_group, ~ve, ~lower_95, ~upper_95, ~source_age_group,
  "2023/24", "all",   0.44,  0.36, 0.51, "All ages (8 months and older)",
  "2023/24", "6m-8",  0.68,  0.51, 0.79, "8 months-8 years",
  "2023/24", "9-17",  0.59,  0.35, 0.75, "9-17 years",
  "2023/24", "18-49", 0.38,  0.24, 0.50, "18-49 years",
  "2023/24", "50-64", 0.16, -0.11, 0.41, "50-64 years",
  "2023/24", "65+",   0.37,  0.05, 0.58, "65 years and older",
  "2024/25", "all",   0.33,  0.24, 0.41, "All ages (8 months and older)",
  "2025/26", "all",   0.24,  0.08, 0.38, "All ages (8 months and older)"
) %>%
  mutate(
    outcome = "medically_attended_influenza_A_or_B",
    source_id = case_when(
      .data$season == "2023/24" ~ "chung_2025_table2",
      .data$season == "2024/25" ~ "chung_2026_abstract_results",
      .data$season == "2025/26" ~ "maloney_2026_table2_us_flu_ve"
    ),
    citation_key = case_when(
      .data$season == "2023/24" ~ "chung2025VE",
      .data$season == "2024/25" ~ "chung2026VE",
      .data$season == "2025/26" ~ "maloney2026VE"
    ),
    source_status = if_else(
      .data$season == "2023/24" & .data$ve_age_group == "6m-8",
      "manual_proxy", "manual"
    ),
    source_url = case_when(
      .data$season == "2023/24" ~
        "https://academic.oup.com/cid/article/81/4/e184/7943530",
      .data$season == "2024/25" ~
        "https://academic.oup.com/cid/advance-article/doi/10.1093/cid/ciag437/8733560",
      .data$season == "2025/26" ~
        "https://www.cdc.gov/mmwr/volumes/75/wr/mm7509a2.htm"
    ),
    source_table = case_when(
      .data$season == "2023/24" ~ "Table 2, all influenza viruses",
      .data$season == "2024/25" ~ "Abstract, Results",
      .data$season == "2025/26" ~
        "Table 2, all ages, any influenza, U.S. Flu VE outpatient"
    ),
    estimate_status = if_else(.data$season == "2025/26", "preliminary", "final")
  )

ve_estimates <- bind_rows(
  scraped %>% filter(
    .data$season != "2010/11",
    !(.data$season == "2015/16" & .data$ve_age_group == "all")
  ),
  manual_2010,
  manual_2015_overall,
  manual_recent
) %>%
  mutate(
    season_start_year = as.integer(str_sub(.data$season, 1, 4)),
    effect_ratio = 1 - .data$ve,
    effect_ratio_lower = 1 - .data$upper_95,
    effect_ratio_upper = 1 - .data$lower_95,
    log_or_mean = log(pmax(.data$effect_ratio, 1e-4)),
    log_or_sd = (log(pmax(.data$effect_ratio_upper, 1e-4)) -
      log(pmax(.data$effect_ratio_lower, 1e-4))) /
      (2 * qnorm(0.975))
  ) %>%
  arrange(.data$season_start_year, .data$ve_age_group)

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_csv(ve_estimates, "data/processed/cdc_ve_estimates.csv")

coverage <- ve_estimates %>%
  count(.data$season, name = "n_age_rows")
dir.create("outputs/diagnostics/data", recursive = TRUE, showWarnings = FALSE)
write_csv(coverage, "outputs/diagnostics/data/cdc_ve_estimate_coverage.csv")

print(coverage)
message("Wrote normalized CDC VE estimates to data/processed/cdc_ve_estimates.csv")
