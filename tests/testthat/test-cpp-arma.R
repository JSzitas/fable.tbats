# Exact-likelihood ARMA against R's arima() on the frozen first-fit errors
# (DESIGN.md, sections 4 and 8).

arma_coef <- function(entry) {
  coef <- entry$coef
  list(ar = unname(coef[grep("^ar", names(coef))]),
       ma = unname(coef[grep("^ma", names(coef))]),
       mean = unname(coef["intercept"]))
}

for (series in fixture_series) {
  fixture <- load_fixture(series)
  e <- fixture$arma$errors
  grid <- Filter(function(g) isTRUE(g$ok), fixture$arma$grid)

  test_that(paste(series, "exact likelihood at R's coefficients"), {
    for (g in grid) {
      coef <- arma_coef(g)
      ll <- tbats_call("tbats_arma_loglik", e, coef$ar, coef$ma, coef$mean)
      expect_equal(ll$loglik, g$loglik, tolerance = 1e-8,
                   label = sprintf("%s ARMA(%d,%d) loglik", series, g$p, g$q))
      expect_equal(ll$sigma2, g$sigma2, tolerance = 1e-6,
                   label = sprintf("%s ARMA(%d,%d) sigma2", series, g$p, g$q))
    }
  })

  test_that(paste(series, "fits are no worse than R's by AIC"), {
    worst <- 0
    for (g in grid) {
      fit <- tbats_call("tbats_arma_fit", e, g$p, g$q, list(method = "lbfgsb"))
      expect_lte(fit$aic, g$aic + 1e-4 * abs(g$aic),
                 label = sprintf("%s ARMA(%d,%d) AIC %.4f vs R %.4f", series, g$p, g$q, fit$aic, g$aic))
      worst <- max(worst, fit$aic - g$aic)
    }
    cat(sprintf("    %-14s largest AIC excess over R across %d orders: %.5f\n", series, length(grid), worst))
  })

  test_that(paste(series, "order selection agrees with R on the frozen grid"), {
    aics <- sapply(grid, function(g) g$aic)
    r_best <- grid[[which.min(aics)]]
    sel <- tbats_call("tbats_arma_select", e, 3, 3, list(method = "lbfgsb"))
    # the likelihood is the same function to 1e-8, so AIC values compare
    # across implementations: the port's best must be no worse than R's best.
    # The orders may differ where R's optimiser stopped short on an order.
    expect_lte(sel$aic, r_best$aic + 1e-4 * abs(r_best$aic))
    cat(sprintf("    %-14s R picks (%d,%d) AIC %.3f; port picks (%d,%d) AIC %.3f\n",
                series, r_best$p, r_best$q, r_best$aic, sel$p, sel$q, sel$aic))
  })
}

test_that("the stepwise order search finds interior optima of the grid", {
  # stepwise cannot promise the grid's best when that lies on the boundary
  # of the orders behind a ridge (on two of the frozen residual series it
  # does); where the grid's best is interior it must find it, and it must
  # always fit fewer orders than the grid. The searches of
  # test-cpp-search.R show the final models are unchanged either way.
  for (series in fixture_series) {
    e <- load_fixture(series)$arma$errors
    grid <- tbats_call("tbats_arma_select", e, 5, 5, list(), "grid")
    step <- tbats_call("tbats_arma_select", e, 5, 5, list(), "stepwise")
    expect_lt(step$fits, grid$fits)
    if (grid$p < 5 && grid$q < 5) {
      expect_equal(c(step$p, step$q), c(grid$p, grid$q), label = paste(series, "stepwise order"))
    }
    cat(sprintf("    %-14s grid (%d,%d) AIC %.3f in %d fits; stepwise (%d,%d) AIC %.3f in %d fits\n",
                series, grid$p, grid$q, grid$aic, grid$fits, step$p, step$q, step$aic, step$fits))
  }
})
