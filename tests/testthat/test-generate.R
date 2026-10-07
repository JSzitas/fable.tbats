library(dplyr)

pelt <- tsibbledata::pelt
models <- fabletools::model(pelt, bats = BATS(Lynx), tbats = TBATS(Lynx))

# a seasonal series exercises the trigonometric (tbats) state space
usa <- tsibble::as_tsibble(USAccDeaths)
seasonal <- fabletools::model(
  usa,
  tbats = TBATS(value ~ parameters(seasonal_periods = 12, trend = TRUE))
)

test_that("generate returns one simulated path per replicate", {
  sims <- fabletools::generate(models, h = 5, times = 3)

  expect_s3_class(sims, "tbl_ts")
  expect_equal(nrow(sims), 2 * 5 * 3)
  expect_true(all(c(".model", ".rep", "Year", ".sim") %in% names(sims)))
  expect_true(all(is.finite(sims$.sim)))
  expect_equal(sort(unique(sims$.rep)), c("1", "2", "3"))
})

test_that("generate is reproducible with a seed", {
  first <- fabletools::generate(models, h = 3, times = 2, seed = 42)
  second <- fabletools::generate(models, h = 3, times = 2, seed = 42)
  expect_equal(first$.sim, second$.sim)
})

test_that("zero innovations reproduce the point forecasts", {
  h <- 4
  fcst <- fabletools::forecast(models, h = h)
  future <- tsibble::new_data(pelt, h)
  future$.innov <- 0
  sims <- fabletools::generate(models, new_data = future)

  for (name in c("bats", "tbats")) {
    expect_equal(sims$.sim[sims$.model == name],
                 fcst$.mean[fcst$.model == name])
  }
})

test_that("feeding the fitted innovations through the recursion reproduces the data", {
  # y_t = w x_{t-1} + e_t with the fitted errors e_t recovers the observed
  # series on the model scale exactly, which checks w, F, g and the state
  # alignment at once
  for (m in list(models$bats[[1]]$fit, models$tbats[[1]]$fit, seasonal$tbats[[1]]$fit)) {
    fit <- m$fit
    y <- fit$fitted + m$resid
    y_model <- if (fit$spec$box_cox) tbats_call("tbats_box_cox", y, fit$parameters$lambda) else y
    reconstructed <- tbats_call("tbats_simulate_model_scale", fit$spec, fit$parameters,
                                fit$seed_states, matrix(fit$errors, ncol = 1))
    expect_equal(as.numeric(reconstructed), y_model)
  }
})

test_that("simulated paths are centred on the forecast", {
  fcst <- fabletools::forecast(models, h = 1)
  sims <- fabletools::generate(models, h = 1, times = 2000, seed = 1)
  # innovations are symmetric on the Box-Cox scale, so the median of the
  # back-transformed paths is the point forecast
  medians <- tapply(sims$.sim, sims$.model, median)
  expect_equal(as.numeric(medians[fcst$.model]), fcst$.mean, tolerance = 0.05)
})

test_that("bootstrap resamples the innovation residuals", {
  sims <- fabletools::generate(models, h = 3, times = 4, bootstrap = TRUE, seed = 7)
  expect_equal(nrow(sims), 2 * 3 * 4)
  expect_true(all(is.finite(sims$.sim)))

  # a supplied .innov column is used as is
  fit <- models$tbats[[1]]$fit$fit
  future <- tsibble::new_data(pelt, 2)
  future$.innov <- c(10, -10)
  one_path <- fabletools::generate(models, new_data = future)
  expected <- tbats_call("tbats_simulate", fit, matrix(c(10, -10), ncol = 1))
  expect_equal(one_path$.sim[one_path$.model == "tbats"], as.numeric(expected))
})

test_that("generate works for a seasonal TBATS model", {
  sims <- fabletools::generate(seasonal, h = 24, times = 2, seed = 3)
  expect_equal(nrow(sims), 24 * 2)
  expect_true(all(is.finite(sims$.sim)))

  fcst <- fabletools::forecast(seasonal, h = 24)
  future <- tsibble::new_data(usa, 24)
  future$.innov <- 0
  sims <- fabletools::generate(seasonal, new_data = future)
  expect_equal(sims$.sim, fcst$.mean)
})
