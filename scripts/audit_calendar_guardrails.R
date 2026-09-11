library(dplyr)
library(readr)

source("R/calendar.R")

ilinet_national <- read_csv(
  "data/processed/ilinet_national_all_age_primary_seasons.csv",
  show_col_types = FALSE
)

season_calendar <- ilinet_national %>%
  distinct(.data$season, .data$season_start_year, .data$mmwr_year, .data$mmwr_week, .data$week_start) %>%
  group_by(.data$season) %>%
  mutate(
    max_week = if_else(any(.data$mmwr_week == 53L), 53L, 52L),
    season_index = week_to_season_index(.data$mmwr_week, max_week = .data$max_week)
  ) %>%
  arrange(.data$season, .data$season_index) %>%
  ungroup()

guardrail_rows <- season_calendar %>%
  filter(.data$mmwr_week %in% c(52L, 53L, 1L, 2L)) %>%
  arrange(.data$season, .data$season_index)

dir.create("outputs/diagnostics/calendar", recursive = TRUE, showWarnings = FALSE)
write_csv(guardrail_rows, "outputs/diagnostics/calendar/calendar_guardrail_week_wrap.csv")

ordinary <- guardrail_rows %>% filter(.data$season == "2011/12")
week53 <- guardrail_rows %>% filter(.data$season == "2014/15")

stopifnot(
  ordinary$mmwr_week[order(ordinary$season_index)] |> identical(c(52, 1, 2)),
  week53$mmwr_week[order(week53$season_index)] |> identical(c(52, 53, 1, 2)),
  weeks_since_full_response(52, 1, max_week = 52) == -1,
  weeks_since_full_response(52, 2, max_week = 52) == 0,
  weeks_since_full_response(52, 1, max_week = 53) == 0
)

print(guardrail_rows %>% filter(.data$season %in% c("2011/12", "2014/15")))
message("Calendar guardrail audit passed.")
