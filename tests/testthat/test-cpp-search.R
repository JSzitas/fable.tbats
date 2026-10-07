# The BATS and TBATS searches against the frozen searches: the chosen model
# must be no worse by AIC, with AIC values comparable because the likelihood
# and the parameter count are the same (DESIGN.md, section 8).

reference_aic <- function(m, y) {
  value <- length(y) * log(sum(m$errors^2))
  if (!is.null(m$lambda)) value <- value - 2 * (as.numeric(m$lambda) - 1) * sum(log(y))
  value + 2 * (length(m$parameters$vect) + length(m$seed_states))
}

describe_spec <- function(spec) {
  sprintf("%s%s%s%s periods (%s) k (%s) ARMA(%d,%d)",
          if (spec$seasonal_type == "dummy") "BATS" else "TBATS",
          if (spec$box_cox) " box-cox" else "", if (spec$trend) " trend" else "",
          if (spec$damping) " damped" else "",
          paste(spec$periods, collapse = ","), paste(spec$harmonics, collapse = ","),
          spec$p, spec$q)
}

describe_reference <- function(m) {
  sprintf("%s%s%s%s periods (%s) k (%s) ARMA(%d,%d)",
          toupper(m$class[1]), if (!is.null(m$lambda)) " box-cox" else "",
          if (!is.null(m$beta)) " trend" else "",
          if (isTRUE(m$parameters$control$use.damping)) " damped" else "",
          paste(m$seasonal_periods, collapse = ","), paste(m$k_vector, collapse = ","),
          length(m$ar), length(m$ma))
}

search_tolerance <- 2e-4

for (series in fixture_series) {
  fixture <- load_fixture(series)
  y <- fixture$y
  for (type in c("tbats", "bats")) {
    m <- fixture[[type]]
    test_that(paste(series, type, "search is no worse than the frozen search"), {
      elapsed <- system.time(
        fit <- tbats_call("tbats_search", y, fixture$periods, type, list())
      )[["elapsed"]]
      reference <- reference_aic(m, y)
      expect_true(is.finite(fit$aic))
      expect_lte(fit$aic, reference + search_tolerance * abs(reference))
      expect_equal(length(fit$fitted), length(y))
      cat(sprintf("    %-13s %-5s reference AIC %10.3f [%s]\n%30s port AIC %10.3f [%s]  %.1fs\n",
                  series, type, reference, describe_reference(m), "", fit$aic,
                  describe_spec(fit$spec), elapsed))
    })
  }
}

test_that("a constant series yields the degenerate model", {
  fixture <- load_fixture("constant")
  for (type in c("tbats", "bats")) {
    fit <- tbats_call("tbats_search", fixture$y, NULL, type, list())
    expect_equal(fit$fitted, fixture$y)
    expect_equal(fit$errors, rep(0, length(fixture$y)))
    expect_equal(fit$variance, 0)
    expect_equal(fit$aic, -Inf)
  }
})

test_that("pinned options are respected", {
  y <- load_fixture("usaccdeaths")$y
  fit <- tbats_call("tbats_search", y, 12, "tbats",
                    list(box_cox = FALSE, trend = TRUE, damping = FALSE, arma_errors = FALSE))
  expect_false(fit$spec$box_cox)
  expect_true(fit$spec$trend)
  expect_false(fit$spec$damping)
  expect_equal(fit$spec$p + fit$spec$q, 0)
})
