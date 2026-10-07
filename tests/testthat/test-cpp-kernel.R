# The structured kernel against the frozen dense matrices, errors, states,
# fitted values and forecasts (DESIGN.md, section 8).

for (series in fixture_series) {
  fixture <- load_fixture(series)
  for (name in names(fixture_models(fixture))) {
    m <- fixture_models(fixture)[[name]]
    spec <- fixture_spec(m)
    par <- fixture_parameters(m)

    test_that(paste(series, name, "kernel reproduces the dense system"), {
      sys <- tbats_call("tbats_dense_system", spec, par)
      expect_equal(sys$w, m$w, tolerance = 1e-12)
      expect_equal(sys$g, m$g, tolerance = 1e-12)
      expect_equal(sys$F, m$F, tolerance = 1e-12)
      # eigenvalues of D: near-zero ones are ill conditioned (multiple zero
      # eigenvalue of the ARMA blocks), hence an absolute tolerance
      expect_lt(max(abs(sys$D_eigen_modulus - m$D_eigen_modulus)), 1e-5)
    })

    test_that(paste(series, name, "parameters pack in the frozen order"), {
      packed <- tbats_call("tbats_pack_parameters", spec, par)
      expect_equal(packed, as.numeric(m$parameters$vect), tolerance = 1e-12)
      unpacked <- tbats_call("tbats_unpack_parameters", spec, packed)
      expect_equal(tbats_call("tbats_pack_parameters", spec, unpacked), packed)
      expect_length(tbats_call("tbats_parameter_scales", spec), length(packed))
    })

    test_that(paste(series, name, "filter reproduces errors, states and fitted values"), {
      y <- fixture$y
      if (spec$box_cox) y <- tbats_call("tbats_box_cox", y, par$lambda)
      flt <- tbats_call("tbats_filter", spec, par, y, m$seed_states)
      expect_equal(flt$errors, m$errors, tolerance = 1e-8)
      expect_equal(flt$states, m$x, tolerance = 1e-8)
      fitted <- flt$predictions
      if (spec$box_cox) fitted <- tbats_call("tbats_inv_box_cox", fitted, par$lambda, NULL)
      expect_equal(fitted, m$fitted, tolerance = 1e-8)
    })

    test_that(paste(series, name, "forecast recursion matches"), {
      x_last <- m$x[, ncol(m$x)]
      fc <- tbats_call("tbats_forecast_model_scale", spec, par, x_last, m$variance, m$forecast_h)
      expect_equal(fc$mean, m$forecast_mean_model_scale, tolerance = 1e-8)
      expect_equal(fc$variance, m$forecast_variance_model_scale, tolerance = 1e-8)
      mean <- fc$mean
      if (spec$box_cox) mean <- tbats_call("tbats_inv_box_cox", mean, par$lambda, NULL)
      expect_equal(mean, m$forecast_mean, tolerance = 1e-8)
      # zero innovations simulate the forecast mean; a unit innovation at the
      # first step moves every later step by w' F^(j-1) g
      zero <- matrix(0, m$forecast_h, 2)
      sims <- tbats_call("tbats_simulate_model_scale", spec, par, x_last, zero)
      expect_equal(sims[, 1], fc$mean, tolerance = 1e-10)
      expect_equal(sims[, 2], fc$mean, tolerance = 1e-10)
    })
  }
}

test_that("Box-Cox transform and inverse round trip, including the log case", {
  y <- c(0.5, 1, 10, 250)
  for (lambda in c(-0.5, 0, 0.3, 1)) {
    z <- tbats_call("tbats_box_cox", y, lambda)
    expect_equal(tbats_call("tbats_inv_box_cox", z, lambda, NULL), y)
  }
  expect_equal(tbats_call("tbats_box_cox", y, 0), log(y))
  # bias adjustment reduces to exp(z)(1 + v/2) for the log transform
  expect_equal(tbats_call("tbats_inv_box_cox", 1, 0, 0.1), exp(1) * 1.05)
})

test_that("Guerrero's lambda matches the frozen values", {
  for (g in load_fixture("guerrero")) {
    y <- load_fixture(g$series)$y
    expect_equal(tbats_call("tbats_guerrero_lambda", y, g$period, g$lower, g$upper),
                 g$lambda, tolerance = 1e-6,
                 label = paste(g$series, g$period, g$lower, g$upper))
  }
})

test_that("specification errors are reported", {
  spec <- list(box_cox = FALSE, trend = FALSE, damping = TRUE, seasonal_type = "dummy",
               periods = numeric(0), harmonics = numeric(0), p = 0, q = 0)
  par <- list(alpha = 0.1, gamma_one = numeric(0), gamma_two = numeric(0), ar = numeric(0), ma = numeric(0))
  expect_error(tbats_call("tbats_dense_system", spec, par), "damping requires a trend")
  spec$damping <- FALSE
  spec$periods <- 12.5
  expect_error(tbats_call("tbats_dense_system", spec, par), "integers")
})
