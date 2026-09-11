# Repository-local HHS region crosswalk for the 50 states and District of
# Columbia. Assignments follow the US Department of Health and Human Services
# regional-office listing:
# https://www.hhs.gov/about/agencies/iea/regional-offices/index.html
# The page reports a review date of August 14, 2024; assignments were checked
# against it on August 22, 2026.

hhs_region_crosswalk <- function() {
  crosswalk <- tibble::tribble(
    ~state_or_territory, ~region_number,
    "Connecticut", 1L,
    "Maine", 1L,
    "Massachusetts", 1L,
    "New Hampshire", 1L,
    "Rhode Island", 1L,
    "Vermont", 1L,
    "New Jersey", 2L,
    "New York", 2L,
    "Delaware", 3L,
    "District of Columbia", 3L,
    "Maryland", 3L,
    "Pennsylvania", 3L,
    "Virginia", 3L,
    "West Virginia", 3L,
    "Alabama", 4L,
    "Florida", 4L,
    "Georgia", 4L,
    "Kentucky", 4L,
    "Mississippi", 4L,
    "North Carolina", 4L,
    "South Carolina", 4L,
    "Tennessee", 4L,
    "Illinois", 5L,
    "Indiana", 5L,
    "Michigan", 5L,
    "Minnesota", 5L,
    "Ohio", 5L,
    "Wisconsin", 5L,
    "Arkansas", 6L,
    "Louisiana", 6L,
    "New Mexico", 6L,
    "Oklahoma", 6L,
    "Texas", 6L,
    "Iowa", 7L,
    "Kansas", 7L,
    "Missouri", 7L,
    "Nebraska", 7L,
    "Colorado", 8L,
    "Montana", 8L,
    "North Dakota", 8L,
    "South Dakota", 8L,
    "Utah", 8L,
    "Wyoming", 8L,
    "Arizona", 9L,
    "California", 9L,
    "Hawaii", 9L,
    "Nevada", 9L,
    "Alaska", 10L,
    "Idaho", 10L,
    "Oregon", 10L,
    "Washington", 10L
  )

  expected_states <- sort(c(state.name, "District of Columbia"))
  if (nrow(crosswalk) != 51L ||
      anyDuplicated(crosswalk$state_or_territory) ||
      !identical(sort(crosswalk$state_or_territory), expected_states) ||
      !setequal(crosswalk$region_number, 1:10)) {
    stop("The repository-local HHS region crosswalk is invalid.", call. = FALSE)
  }

  dplyr::mutate(
    crosswalk,
    region = paste("Region", .data$region_number)
  )
}
