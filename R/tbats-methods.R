#' @export
print.TBATS <- function(x, ...) {
  cat(x[["model_summary"]])
}

#' Residuals of a BATS or TBATS model
#'
#' @param object A fitted BATS or TBATS model.
#' @param type \code{"innovation"} (the default) for the one step ahead errors
#' of the state space model, which are on the Box-Cox transformed scale when the
#' model uses a Box-Cox transformation, or \code{"response"} for the difference
#' between the observed and fitted values on the original scale. Other types
#' are not supported and return \code{NULL}, which fabletools reports.
#' @param ... Unused.
#' @return A numeric vector of residuals.
#' @importFrom stats residuals
#' @export
residuals.TBATS <- function(object, type = "innovation", ...) {
  bats_tbats_residuals(object, type)
}

#' @importFrom stats fitted
#' @export
fitted.TBATS <- function(object, ...) {
  object[["fitted"]]
}

#' @importFrom fabletools forecast
#' @export
forecast.TBATS <- function(object, new_data = NULL, specials = NULL, bootstrap = FALSE,
                           times = 5000, ...) {
  forecast_bats_tbats(object, new_data)
}

#' Refit a BATS or TBATS model to new data
#'
#' @param object A fitted BATS or TBATS model.
#' @param new_data A tsibble holding the new series.
#' @param specials Unused.
#' @param reestimate If \code{FALSE} (the default) the estimated parameters are
#' kept and only the seed states, the state vector before the first
#' observation, are re-estimated for the new series. If \code{TRUE} the full
#' model search is repeated with the stored options. This departs from the
#' forecast package, which reuses the old seed states unchanged.
#' @param ... Unused.
#' @return A fitted model of the same class.
#' @importFrom generics refit
#' @export
refit.TBATS <- function(object, new_data, specials = NULL, reestimate = FALSE, ...) {
  refit_bats_tbats(object, new_data, reestimate)
}

#' @importFrom generics glance
#' @export
glance.TBATS <- function(x, ...) {
  glance_bats_tbats(x)
}

#' Components of a BATS or TBATS model
#'
#' @description The level, the slope (with a trend) and one seasonal effect
#' per seasonal period at every observation, with the innovation as the
#' remainder. For a Box-Cox model the components are on the transformed
#' scale, so the response is not their sum.
#' @param object A fitted BATS or TBATS model.
#' @param ... Unused.
#' @return A dable.
#' @importFrom fabletools components
#' @export
components.TBATS <- function(object, ...) {
  components_bats_tbats(object)
}
