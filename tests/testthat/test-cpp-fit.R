# Seed states, likelihood and single-specification fits against the frozen
# implementation (DESIGN.md, section 8).

# the forecast package's starting values, as a parameter list for the glue
reference_start <- function(spec, init_lambda) {
  dummy <- spec$seasonal_type == "dummy"
  long_period <- dummy && sum(spec$periods) > 16
  m <- length(spec$periods)
  list(
    lambda = if (spec$box_cox) init_lambda else 1,
    alpha = if (long_period) 1e-6 else 0.09,
    beta = if (!spec$trend) 0 else if (long_period) 5e-7 else 0.05,
    phi = if (spec$damping) 0.999 else 1,
    gamma_one = rep(if (dummy) 0.001 else 0, m),
    gamma_two = if (dummy) numeric(0) else rep(0, m),
    ar = rep(0, spec$p), ma = rep(0, spec$q)
  )
}

# -2 log L of a frozen model from its own errors. For one frozen model (the
# two-period BATS with ARMA errors on the multiseasonal series) the reference
# stored an optimiser value that disagrees with its final errors: its in-place
# update of the transition matrix skips the second seasonal block's coupling
# to the ARMA states, so the objective it optimised was not the likelihood of
# the model it returned. The errors are the authority.
reference_neg2loglik <- function(m, y) {
  value <- length(y) * log(sum(m$errors^2))
  # as.numeric strips the biasadj attribute the reference hangs on lambda
  if (!is.null(m$lambda)) value <- value - 2 * (as.numeric(m$lambda) - 1) * sum(log(y))
  value
}

guerrero_value <- function(series, period, lower = 0, upper = 1) {
  for (g in load_fixture("guerrero")) {
    if (g$series == series && g$period == period && g$lower == lower && g$upper == upper) {
      return(g$lambda)
    }
  }
  stop("no frozen Guerrero value for ", series, " period ", period)
}

for (series in fixture_series) {
  fixture <- load_fixture(series)
  y <- fixture$y

  # fixed specifications were frozen from fitSpecific* with the Guerrero
  # lambda of the plain numeric series (period 2, bounds 0 and 1)
  for (name in names(fixture$fixed)) {
    m <- fixture$fixed[[name]]
    spec <- fixture_spec(m)
    init_lambda <- guerrero_value(series, 2)

    test_that(paste(series, name, "seed states match at the starting parameters"), {
      start <- reference_start(spec, init_lambda)
      y_model <- if (spec$box_cox) tbats_call("tbats_box_cox", y, init_lambda) else y
      x0 <- tbats_call("tbats_seed_states", spec, start, y_model)
      if (spec$box_cox) {
        # the reference carries the seed to the final lambda through the
        # original scale
        x0 <- tbats_call("tbats_box_cox",
                         tbats_call("tbats_inv_box_cox", x0, init_lambda, NULL),
                         m$lambda)
      }
      expect_equal(x0, m$seed_states, tolerance = 1e-6)
    })
  }

  for (name in names(fixture_models(fixture))) {
    m <- fixture_models(fixture)[[name]]
    spec <- fixture_spec(m)
    par <- fixture_parameters(m)

    test_that(paste(series, name, "likelihood at the frozen parameters"), {
      ll <- tbats_call("tbats_neg2loglik", spec, par, y, m$seed_states)
      expect_true(ll$admissible)
      reference <- reference_neg2loglik(m, y)
      expect_equal(ll$neg2loglik, reference, tolerance = 1e-8)
      expect_equal(ll$aic, reference + 2 * (length(m$parameters$vect) + length(m$seed_states)),
                   tolerance = 1e-8)
    })
  }
}

# Fits from scratch must be no worse than the frozen optimum beyond a
# tolerance. The frozen search models used the Guerrero lambda of the series
# at its largest seasonal period (2 when non-seasonal); the fixed
# specifications used period 2.
fit_cases <- list()
for (series in fixture_series) {
  fixture <- load_fixture(series)
  search_period <- if (is.null(fixture$periods)) 2 else floor(max(fixture$periods))
  for (name in c("tbats", "bats")) {
    fit_cases[[paste(series, name)]] <- list(series = series, m = fixture[[name]], period = search_period)
  }
  for (name in names(fixture$fixed)) {
    fit_cases[[paste(series, name)]] <- list(series = series, m = fixture$fixed[[name]], period = 2)
  }
}

# a fit is accepted when its -2 log L is within this fraction of the
# reference value above it; with two Nelder-Mead restarts the port lands
# within 0.01 percent of every frozen optimum, and below most of them
fit_tolerance <- 2e-4

for (case_name in names(fit_cases)) {
  case <- fit_cases[[case_name]]
  m <- case$m
  if (isTRUE(m$degenerate)) next
  test_that(paste(case_name, "fit from scratch is no worse than the reference"), {
    spec <- fixture_spec(m)
    y <- load_fixture(case$series)$y
    init_lambda <- if (spec$box_cox) guerrero_value(case$series, case$period) else 1
    elapsed <- system.time(
      fit <- tbats_call("tbats_fit_specific", spec, y, init_lambda, list(method = "nelder_mead"), FALSE)
    )[["elapsed"]]
    reference <- reference_neg2loglik(m, y)
    expect_true(is.finite(fit$neg2loglik))
    expect_lte(fit$neg2loglik, reference + fit_tolerance * abs(reference))
    cat(sprintf("    %-45s reference %10.3f  port %10.3f  (%d evaluations, %.1fs)\n",
                case_name, reference, fit$neg2loglik, fit$evaluations, elapsed))
  })
}
