train_tbats <- function(.data, specials, ...) {
  train_bats_tbats(.data, specials, "tbats")
}

specials_tbats <- fabletools::new_specials(
  parameters = function(trend = NULL,
                        damped = NULL,
                        box_cox = NULL,
                        seasonal_periods = "auto",
                        arma_errors = TRUE,
                        bias_adj = FALSE,
                        bc_lower = 0,
                        bc_higher = 1) {
    list(
      use.box.cox = box_cox,
      use.trend = trend,
      use.damped.trend = damped,
      seasonal.periods = seasonal_periods,
      use.arma.errors = arma_errors,
      bc.lower = bc_lower,
      bc.upper = bc_higher,
      biasadj = bias_adj
    )
  },
  xreg = function(...) {
    stop("Exogenous regressors aren't supported by `TBATS()`")
  },
  .required_specials = c("parameters")
)

#' TBATS model
#'
#' @description Exponential smoothing state space model with Box-Cox
#' transformation, ARMA errors, trend and trigonometric seasonal components,
#' with automatic selection of every part, as a fable model. The fitting is
#' implemented in C++ in this package and follows De Livera, Hyndman and
#' Snyder (2011).
#' @param formula A model formula; the response may be followed by a
#' \code{parameters()} special (see details).
#' @param ... Further arguments passed to fabletools.
#' @return A model definition, used like any other fable model.
#' @details The \code{parameters()} special accepts \code{trend},
#' \code{damped} and \code{box_cox} (each \code{NULL} to let the search
#' decide, or \code{TRUE}/\code{FALSE} to pin), \code{seasonal_periods}
#' (\code{"auto"} to detect them with \code{\link{find_seasonalities}},
#' \code{NULL} for the period the index implies, or a numeric vector),
#' \code{arma_errors}, \code{bias_adj}, \code{bc_lower} and \code{bc_higher}.
#' With a Box-Cox transformation the forecast distribution is the normal
#' forecast on the transformed scale mapped back through the inverse
#' transformation, so its quantiles are exact and its mean is the
#' bias-adjusted mean; \code{bias_adj} controls the fitted values only.
#' @references De Livera, A. M., Hyndman, R. J. and Snyder, R. D. (2011).
#' Forecasting time series with complex seasonal patterns using exponential
#' smoothing. Journal of the American Statistical Association 106(496),
#' 1513-1527.
#' @export
TBATS <- function(formula, ...) {
  model_tbats <- fabletools::new_model_class("TBATS",
                                             train = train_tbats,
                                             specials = specials_tbats,
                                             check = function(.data) {
                                               if (!tsibble::is_regular(.data)) stop("Data must be regular")
                                             })
  fabletools::new_model_definition(model_tbats, !!rlang::enquo(formula), ...)
}
