spencer_relative_ve_retention <- function(weeks_since_response, curve) {
  curve <- match.arg(curve, c("fast", "slow"))
  # The cubic coefficients use two-week time units. Immune-response lag is
  # applied separately by the decision engine.
  biweeks <- pmax(weeks_since_response, 0) / 2

  decline <- if (curve == "fast") {
    55 - 1.37 * biweeks + 0.18 * biweeks^2 - 0.03 * biweeks^3
  } else {
    55 - 0.50 * biweeks + 0.05 * biweeks^2 - 0.01 * biweeks^3
  }

  pmin(pmax(decline / 55, 0), 1)
}

protection_model_from_prior <- function(prior) {
  switch(
    prior,
    spencer_ferdinands_fast = "spencer_fast_relative_ve",
    spencer_slow = "spencer_slow_relative_ve",
    "exponential_effect_ratio"
  )
}
