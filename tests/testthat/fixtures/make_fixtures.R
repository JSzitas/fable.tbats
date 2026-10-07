# Freezes the behaviour of the R and Armadillo implementation (fable.tbats
# 0.5.0) into fixtures for the C++ port. The port is compared against these
# with tolerances; see src/tbats/DESIGN.md, section 8.
#
# Run from the package root with the 0.5.0 sources present:
#   Rscript tests/testthat/fixtures/make_fixtures.R
# It cannot be rerun once the vendored R implementation is deleted, which is
# the point: these files are the frozen oracle.

pkgload::load_all(quiet = TRUE)
out_dir <- "tests/testthat/fixtures"

# ---- series ---------------------------------------------------------------

synthetic <- function(n, periods, amplitudes, level, slope, ar, sd, seed) {
  set.seed(seed)
  t <- seq_len(n)
  seasonal <- Reduce(`+`, Map(function(m, a) a * sin(2 * pi * t / m), periods, amplitudes), 0)
  noise <- as.numeric(stats::arima.sim(list(ar = ar), n = n, sd = sd))
  level + slope * t + seasonal + noise
}

series <- list(
  lynx = list(y = tsibbledata::pelt$Lynx, periods = NULL),
  usaccdeaths = list(y = as.numeric(USAccDeaths), periods = 12),
  multiseasonal = list(
    y = synthetic(400, c(7, 30), c(10, 15), level = 100, slope = 0.05, ar = 0.5, sd = 3, seed = 1),
    periods = c(7, 30)
  ),
  # crosses zero, so Box-Cox must be disabled by the search
  nonpositive = list(
    y = synthetic(120, 12, 10, level = 0, slope = 0, ar = 0.3, sd = 4, seed = 2),
    periods = 12
  ),
  constant = list(y = rep(5, 50), periods = NULL)
)

# ---- snapshot of one fitted bats / tbats object ---------------------------

# transition on the model scale: y_{t} = w' x_{t-1}, x_t = F x_{t-1}
forecast_model_scale <- function(mats, x_last, variance, h) {
  x <- x_last
  mean <- numeric(h)
  for (i in seq_len(h)) {
    mean[i] <- drop(mats$w.transpose %*% x)
    x <- mats$F %*% x
  }
  # c_j = w' F^(j-1) g; var_h = sigma^2 (1 + sum_{j<h} c_j^2)
  v <- mats$g
  cj <- numeric(h - 1)
  for (j in seq_len(h - 1)) {
    cj[j] <- drop(mats$w.transpose %*% v)
    v <- mats$F %*% v
  }
  list(mean = mean, variance = variance * (1 + c(0, cumsum(cj^2))))
}

snapshot <- function(m, h = 12) {
  if (is.null(m$damping.parameter) && is.null(m$parameters)) {
    # degenerate constant-series model
    return(list(class = class(m), degenerate = TRUE, x = m$x,
                errors = as.numeric(m$errors), fitted = as.numeric(m$fitted.values),
                variance = m$variance, AIC = m$AIC, alpha = m$alpha))
  }
  mats <- state_space_matrices(m)
  D <- mats$F - mats$g %*% mats$w.transpose
  fc <- forecast_model_scale(mats, m$x[, ncol(m$x)], m$variance, h)
  list(
    class = class(m), degenerate = FALSE,
    lambda = m$lambda, alpha = m$alpha, beta = m$beta, damping = m$damping.parameter,
    gamma = m$gamma.values, gamma_one = m$gamma.one.values, gamma_two = m$gamma.two.values,
    ar = as.numeric(m$ar.coefficients), ma = as.numeric(m$ma.coefficients),
    seasonal_periods = m$seasonal.periods, k_vector = m$k.vector,
    parameters = m$parameters,
    seed_states = as.numeric(m$seed.states),
    x = m$x,
    errors = as.numeric(m$errors), fitted = as.numeric(m$fitted.values),
    variance = m$variance, likelihood = m$likelihood, AIC = m$AIC,
    w = as.numeric(mats$w.transpose), F = mats$F, g = as.numeric(mats$g),
    D_eigen_modulus = sort(Mod(eigen(D, only.values = TRUE)$values)),
    forecast_h = h,
    forecast_mean_model_scale = fc$mean,
    forecast_variance_model_scale = fc$variance,
    forecast_mean = as.numeric(forecast_tbats(m, h = h, level = 80)$mean)
  )
}

# ---- ARMA order selection oracle ------------------------------------------

arma_grid <- function(e, max_order = 3) {
  rows <- list()
  for (p in 0:max_order) for (q in 0:max_order) {
    fit <- try(suppressWarnings(stats::arima(e, order = c(p, 0, q), method = "ML")), silent = TRUE)
    rows[[length(rows) + 1]] <- if (inherits(fit, "try-error")) {
      list(p = p, q = q, ok = FALSE)
    } else {
      list(p = p, q = q, ok = TRUE, loglik = fit$loglik, aic = fit$aic,
           sigma2 = fit$sigma2, coef = fit$coef)
    }
  }
  chosen <- auto_arma(e)
  list(errors = e, grid = rows, chosen = c(p = chosen$arma[1], q = chosen$arma[2]))
}

# ---- fixed specifications (fitted without the search) ----------------------

fixed_specs <- function(y, periods) {
  positive <- all(y > 0)
  init_lambda <- if (positive) BoxCox.lambda(y, lower = 0, upper = 1) else NULL
  specs <- list()
  if (is.null(periods)) {
    specs$bats_trend <- fitSpecificBATS(y, use.box.cox = FALSE, use.beta = TRUE, use.damping = FALSE)
    if (positive) {
      specs$bats_boxcox_damped_arma11 <- fitSpecificBATS(
        y, use.box.cox = TRUE, use.beta = TRUE, use.damping = TRUE,
        ar.coefs = 0, ma.coefs = 0, init.box.cox = init_lambda)
    }
  } else {
    k1 <- rep(1, length(periods))
    specs$tbats_trend_k1 <- fitSpecificTBATS(
      y, use.box.cox = FALSE, use.beta = TRUE, use.damping = FALSE,
      seasonal.periods = periods, k.vector = k1)
    specs$bats_trend <- fitSpecificBATS(
      y, use.box.cox = FALSE, use.beta = TRUE, use.damping = FALSE,
      seasonal.periods = periods)
    if (positive) {
      specs$tbats_boxcox_damped_k2_arma11 <- fitSpecificTBATS(
        y, use.box.cox = TRUE, use.beta = TRUE, use.damping = TRUE,
        seasonal.periods = periods, k.vector = pmin(2, floor((periods - 1) / 2)),
        ar.coefs = 0, ma.coefs = 0, init.box.cox = init_lambda)
    }
  }
  lapply(specs, snapshot)
}

# ---- main -------------------------------------------------------------------

for (name in names(series)) {
  y <- series[[name]]$y
  periods <- series[[name]]$periods
  cat("==", name, "\n")
  t0 <- Sys.time()
  fixture <- list(
    package_version = as.character(packageVersion("fable.tbats")),
    y = y, periods = periods,
    tbats = snapshot(tbats(y, seasonal.periods = periods, use.parallel = FALSE)),
    bats = snapshot(bats(y, seasonal.periods = periods, use.parallel = FALSE))
  )
  if (name != "constant") {
    fixture$fixed <- fixed_specs(y, periods)
    first <- fixture$fixed[[1]]
    fixture$arma <- arma_grid(first$errors)
  }
  saveRDS(fixture, file.path(out_dir, paste0(name, ".rds")), version = 2)
  cat("   tbats:", fixture$tbats$class[1], " AIC", fixture$tbats$AIC,
      "  bats AIC", fixture$bats$AIC,
      "  elapsed", format(round(Sys.time() - t0, 1)), "\n")
}
