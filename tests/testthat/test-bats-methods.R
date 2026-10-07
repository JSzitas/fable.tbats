library(dplyr)

pelt <- tsibbledata::pelt
train <- pelt %>% dplyr::filter(Year < 1930)
model <- fabletools::model(pelt, bats = BATS(Lynx))
fit <- model$bats[[1]][["fit"]]

test_that("Utilities for BATS work", {
  expect_equal(residuals(fit), as.numeric(fit[["fit"]][["errors"]]))
  expect_equal(residuals(fit, type = "response"), pelt$Lynx - fitted(fit))
  expect_length(fitted(fit), nrow(pelt))
  g <- fabletools::glance(model)
  expect_true(all(is.finite(unlist(g[, c("sigma2", "log_lik", "AIC", "AICc", "BIC")]))))
  expect_output(print(fit), "^BATS\\(")
})

test_that("Forecasts for BATS work", {
  fcst <- fabletools::forecast(model, h = 3)
  expect_equal(nrow(fcst), 3)
  expect_true(all(is.finite(fcst$.mean)))
  expect_equal(fcst$.mean, c(36004, 31556, 24992), tolerance = 0.1)
})

test_that("Refitting a BATS works", {
  partial <- fabletools::model(train, bats = BATS(Lynx))
  kept <- generics::refit(partial, pelt)
  expect_s3_class(kept$bats[[1]][["fit"]], "BATS")
  expect_equal(kept$bats[[1]][["fit"]][["fit"]][["parameters"]],
               partial$bats[[1]][["fit"]][["fit"]][["parameters"]])
  expect_length(kept$bats[[1]][["fit"]][["fitted"]], nrow(pelt))
  searched <- generics::refit(partial, pelt, reestimate = TRUE)
  expect_s3_class(searched$bats[[1]][["fit"]], "BATS")
  expect_equal(searched$bats[[1]][["fit"]][["model_summary"]], fit[["model_summary"]])
})

test_that("Components of a BATS model form a dable", {
  usa <- tsibble::as_tsibble(USAccDeaths)
  seasonal <- fabletools::model(usa, bats = BATS(value ~ parameters(trend = TRUE)))
  cmp <- fabletools::components(seasonal)
  expect_s3_class(cmp, "dcmp_ts")
  expect_true(all(c("level", "slope", "season_12", "remainder") %in% names(cmp)))
  expect_equal(nrow(cmp), nrow(usa))
})
