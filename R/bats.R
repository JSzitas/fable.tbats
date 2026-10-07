train_bats <- function(.data, specials, ...) {
  train_bats_tbats(.data, specials, "bats")
}

specials_bats <- fabletools::new_specials(
  parameters = function(trend = NULL,
                        damped = NULL,
                        box_cox = NULL,
                        seasonal_periods = NULL,
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
    stop("Exogenous regressors aren't supported by `BATS()`")
  },
  .required_specials = c("parameters")
)

#' BATS model
#'
#' @description Exponential smoothing state space model with Box-Cox
#' transformation, ARMA errors, trend and seasonal components represented by
#' one state per position in each (integer) period, with automatic selection
#' of every part, as a fable model. See \code{\link{TBATS}} for the
#' trigonometric representation, which also accepts non-integer periods.
#' @inheritParams TBATS
#' @return A model definition, used like any other fable model.
#' @details The \code{parameters()} special is that of \code{\link{TBATS}};
#' \code{seasonal_periods} defaults to the period the index implies.
#' @export
BATS <- function(formula, ...) {
  model_bats <- fabletools::new_model_class("BATS",
                                            train = train_bats,
                                            specials = specials_bats,
                                            check = function(.data) {
                                              if (!tsibble::is_regular(.data)) stop("Data must be regular")
                                            })
  fabletools::new_model_definition(model_bats, !!rlang::enquo(formula), ...)
}
