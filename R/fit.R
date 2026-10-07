# The fitted model object shared by BATS and TBATS, and the pieces of the
# fable interface both model classes delegate to. The fitting itself happens
# in C++ (src/tbats/tbats.h); the R side converts, stores and dispatches.

# Seasonal periods for a series: the user's choice, or the periods detected
# by find_seasonalities(), or the frequency the tsibble's index implies.
seasonal_periods_for <- function(.data, y, requested) {
  periods <- if (is.null(requested)) {
    stats::frequency(stats::as.ts(.data))
  } else if (identical(requested, "auto")) {
    find_seasonalities(y)
  } else {
    requested
  }
  periods <- as.numeric(periods)
  periods[periods > 1]
}

search_options_for <- function(parameters) {
  list(
    box_cox = parameters$use.box.cox,
    trend = parameters$use.trend,
    damping = parameters$use.damped.trend,
    arma_errors = isTRUE(parameters$use.arma.errors),
    box_cox_lower = parameters$bc.lower,
    box_cox_upper = parameters$bc.upper,
    bias_adjust = isTRUE(parameters$biasadj)
  )
}

# Trains either model class on one series; `type` is "bats" or "tbats".
train_bats_tbats <- function(.data, specials, type) {
  parameters <- specials$parameters[[1]]
  target <- tsibble::measured_vars(.data)
  y <- unclass(.data)[[target]]
  if (anyNA(y)) {
    stop(toupper(type), " does not support missing values.", call. = FALSE)
  }
  periods <- seasonal_periods_for(.data, y, parameters$seasonal.periods)
  if (type == "bats" && any(periods != round(periods))) {
    stop("BATS needs integer seasonal periods; TBATS accepts any period.", call. = FALSE)
  }
  fit <- tbats_search(y, periods, type, search_options_for(parameters))
  new_bats_tbats_fit(fit, y, .data, target, parameters, type)
}

# One constructor for the model object, used by training and both refit paths.
new_bats_tbats_fit <- function(fit, y, .data, target, parameters, type) {
  structure(
    list(
      fit = fit,
      resid = y - fit$fitted,
      fitted = fit$fitted,
      target = target,
      index = .data[[tsibble::index_var(.data)]],
      index_var = tsibble::index_var(.data),
      model_summary = model_summary(fit, type),
      model_pars = parameters,
      type = type
    ),
    class = toupper(type)
  )
}

# The summary string of the forecast package: lambda (1 without Box-Cox),
# the ARMA orders, the damping (- without a trend) and the seasonal periods,
# with their harmonics for TBATS.
model_summary <- function(fit, type) {
  spec <- fit$spec
  par <- fit$parameters
  lambda <- if (spec$box_cox) round(par$lambda, 3) else "1"
  phi <- if (spec$trend) round(par$phi, 3) else "-"
  seasonal <- if (length(spec$periods) == 0) {
    if (type == "tbats") "{-}" else "-"
  } else if (type == "tbats") {
    paste0("{", paste0("<", round(spec$periods, 2), ",", spec$harmonics, ">", collapse = ", "), "}")
  } else {
    paste0("{", paste(spec$periods, collapse = ","), "}")
  }
  sprintf("%s(%s, {%d,%d}, %s, %s)", toupper(type), lambda, spec$p, spec$q, phi, seasonal)
}

bats_tbats_residuals <- function(object, type) {
  switch(type,
         innovation = as.numeric(object[["fit"]][["errors"]]),
         response = object[["resid"]],
         NULL)
}

# The forecast distribution: normal on the model scale, mapped back through
# the inverse Box-Cox transformation when the model uses one, so quantiles
# are exact and the mean is the bias-adjusted mean the transformation implies.
forecast_bats_tbats <- function(object, new_data) {
  fit <- object[["fit"]]
  fc <- tbats_forecast(fit, nrow(new_data), 80, isTRUE(object$model_pars$biasadj))
  model_scale <- distributional::dist_normal(fc$mean_model_scale, sqrt(fc$variance_model_scale))
  if (!fit$spec$box_cox) return(model_scale)
  lambda <- fit$parameters$lambda
  distributional::dist_transformed(model_scale,
                                   transform = function(x) fabletools::inv_box_cox(x, lambda),
                                   inverse = function(x) fabletools::box_cox(x, lambda))
}

generate_bats_tbats <- function(object, new_data) {
  if (!tsibble::is_regular(new_data)) {
    stop("Simulation new_data must be regularly spaced")
  }
  fit <- object[["fit"]]
  if (!(".innov" %in% names(new_data))) {
    new_data$.innov <- stats::rnorm(NROW(new_data), sd = sqrt(fit$variance))
  }
  # one path per key (fabletools keys each path by .rep), rows in time order
  index <- new_data[[tsibble::index_var(new_data)]]
  paths <- lapply(tsibble::key_rows(new_data), function(rows) rows[order(index[rows])])
  h <- max(lengths(paths))
  innov <- matrix(0, nrow = h, ncol = length(paths))
  for (k in seq_along(paths)) {
    innov[seq_along(paths[[k]]), k] <- new_data$.innov[paths[[k]]]
  }
  sims <- tbats_simulate(fit, innov)
  sim <- numeric(NROW(new_data))
  for (k in seq_along(paths)) {
    sim[paths[[k]]] <- sims[seq_along(paths[[k]]), k]
  }
  new_data$.sim <- sim
  new_data
}

# Refit to new data: with reestimate the whole search is repeated with the
# stored options; otherwise the parameters are kept and the seed states are
# re-estimated for the new series.
refit_bats_tbats <- function(object, new_data, reestimate) {
  type <- object[["type"]]
  if (reestimate) {
    return(train_bats_tbats(new_data, list(parameters = list(object[["model_pars"]])), type))
  }
  target <- tsibble::measured_vars(new_data)
  y <- unclass(new_data)[[target]]
  if (anyNA(y)) {
    stop(toupper(type), " does not support missing values.", call. = FALSE)
  }
  fit <- tbats_refit(object[["fit"]], y, isTRUE(object$model_pars$biasadj))
  new_bats_tbats_fit(fit, y, new_data, target, object[["model_pars"]], type)
}

glance_bats_tbats <- function(x) {
  fit <- x[["fit"]]
  resid <- x[["resid"]]
  n <- length(resid)
  n_pars <- fit$n_parameters
  dof <- n - n_pars
  loglik <- -fit$neg2loglik / 2
  tsibble::tibble(
    sigma2 = sum(resid^2) / dof,
    log_lik = loglik,
    AIC = fit$aic,
    AICc = -2 * loglik + 2 * (n_pars + 1) * (n / (n - 1 - n_pars)),
    BIC = -2 * loglik + n_pars * log(n),
    dof = dof
  )
}

# Level, slope and one seasonal effect per period on the model scale, with
# the innovation as remainder, as a dable. For a Box-Cox model the
# components are on the transformed scale.
components_bats_tbats <- function(object) {
  fit <- object[["fit"]]
  cmp <- tbats_components(fit)
  columns <- c(
    stats::setNames(list(object[["index"]]), object[["index_var"]]),
    stats::setNames(list(object[["fit"]][["fitted"]] + object[["resid"]]), object[["target"]]),
    lapply(seq_len(ncol(cmp)), function(j) cmp[, j]),
    list(remainder = as.numeric(fit$errors))
  )
  names(columns)[2 + seq_len(ncol(cmp))] <- colnames(cmp)
  out <- tsibble::tsibble(!!!columns, index = !!rlang::sym(object[["index_var"]]))
  seasons <- lapply(fit$spec$periods, function(m) list(period = m, base = NA_real_))
  names(seasons) <- grep("^season_", colnames(cmp), value = TRUE)
  fabletools::as_dable(out, response = !!rlang::sym(object[["target"]]),
                       method = object[["model_summary"]], seasons = seasons)
}
