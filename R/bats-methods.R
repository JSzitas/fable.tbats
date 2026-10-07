#' @export
print.BATS <- function(x, ...) {
  cat(x[["model_summary"]])
}

#' @rdname residuals.TBATS
#' @importFrom stats residuals
#' @export
residuals.BATS <- function(object, type = "innovation", ...) {
  bats_tbats_residuals(object, type)
}

#' @importFrom stats fitted
#' @export
fitted.BATS <- function(object, ...) {
  object[["fitted"]]
}

#' @importFrom fabletools forecast
#' @export
forecast.BATS <- function(object, new_data = NULL, specials = NULL, bootstrap = FALSE,
                          times = 5000, ...) {
  forecast_bats_tbats(object, new_data)
}

#' @rdname refit.TBATS
#' @importFrom generics refit
#' @export
refit.BATS <- function(object, new_data, specials = NULL, reestimate = FALSE, ...) {
  refit_bats_tbats(object, new_data, reestimate)
}

#' @importFrom generics glance
#' @export
glance.BATS <- function(x, ...) {
  glance_bats_tbats(x)
}

#' @rdname components.TBATS
#' @importFrom fabletools components
#' @export
components.BATS <- function(object, ...) {
  components_bats_tbats(object)
}
