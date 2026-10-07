library(dplyr)

pelt <- tsibbledata::pelt
train <- pelt %>% dplyr::filter(Year < 1930)
model <- fabletools::model(pelt, tbats = TBATS(Lynx))
fit <- model$tbats[[1]][["fit"]]

test_that("Utilities for TBATS work", {
  # innovation residuals by default, response residuals on request
  expect_equal(residuals(fit), as.numeric(fit[["fit"]][["errors"]]))
  expect_equal(residuals(fit, type = "response"), pelt$Lynx - fitted(fit))
  expect_null(residuals(fit, type = "deviance"))
  expect_length(fitted(fit), nrow(pelt))
  expect_equal(sum(residuals(model, type = "response")$.resid), sum(residuals(fit, type = "response")))
  g <- fabletools::glance(model)
  expect_true(all(is.finite(unlist(g[, c("sigma2", "log_lik", "AIC", "AICc", "BIC")]))))
  expect_equal(g$dof, nrow(pelt) - fit[["fit"]][["n_parameters"]])
  expect_output(print(fit), "^TBATS\\(")
})

test_that("Forecasts for TBATS work", {
  fcst <- fabletools::forecast(model, h = 3)
  expect_equal(nrow(fcst), 3)
  expect_true(all(is.finite(fcst$.mean)))
  # the search used to land near these values; the port lands within a few
  # percent of them
  expect_equal(fcst$.mean, c(36004, 31556, 24992), tolerance = 0.1)
  expect_true(all(distributional::variance(fcst$Lynx) > 0))
})

test_that("Refitting a TBATS works", {
  partial <- fabletools::model(train, tbats = TBATS(Lynx))
  kept <- generics::refit(partial, pelt)
  expect_equal(kept$tbats[[1]][["fit"]][["fit"]][["parameters"]],
               partial$tbats[[1]][["fit"]][["fit"]][["parameters"]])
  expect_length(kept$tbats[[1]][["fit"]][["fitted"]], nrow(pelt))
  expect_equal(kept$tbats[[1]][["fit"]][["model_summary"]], partial$tbats[[1]][["fit"]][["model_summary"]])
  searched <- generics::refit(partial, pelt, reestimate = TRUE)
  expect_equal(searched$tbats[[1]][["fit"]][["model_summary"]], fit[["model_summary"]])
})

test_that("Components of a TBATS model form a dable", {
  usa <- tsibble::as_tsibble(USAccDeaths)
  seasonal <- fabletools::model(usa, tbats = TBATS(value ~ parameters(seasonal_periods = 12, trend = TRUE)))
  cmp <- fabletools::components(seasonal)
  expect_s3_class(cmp, "dcmp_ts")
  expect_true(all(c("value", "level", "slope", "season_12", "remainder") %in% names(cmp)))
  expect_equal(nrow(cmp), nrow(usa))
  expect_equal(cmp$remainder, as.numeric(seasonal$tbats[[1]][["fit"]][["fit"]][["errors"]]))
})
