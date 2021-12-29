train_bats <- function(.data, specials, ...) {

  # parse arguments to bats
  parameters <- specials$parameters[[1]]
  y <- unclass(.data)[[tsibble::measured_vars(.data)]]

  if( is.null( parameters$seasonal.periods )) {
    y <- stats::as.ts(.data)
  }
  else if ( parameters$seasonal.periods == "auto" ) {
    parameters$seasonal.periods <- find_seasonalities( y )
  }

  # always set use.parallel to FALSE - since nested parallelism would cause problems
  # and the ONLY situatiion where we avoid that is when someone is running a single
  # TBATS model on a single time series, OR running all models sequentially
  # by setting future::plan("sequential") - so they are probably not too worried
  # about this being slow.
  model <- do.call( forecast::bats,
                    c( list(y, use.parallel = FALSE),
                       parameters)
  )

  structure(
    list(
      fit = model,
      resid = c(y) - stats::fitted(model),
      fitted = stats::fitted(model),
      target = tsibble::measured_vars(.data),
      model_summary = as.character(model),
      model_pars = parameters
    ),
    class = "BATS"
  )
}

specials_bats <- fabletools::new_specials(
  parameters = function( trend = NULL,
                         damped = NULL,
                         box_cox = NULL,
                         seasonal_periods = NULL,
                         arma_errors = TRUE,
                         bias_adj = FALSE,
                         bc_lower = 0,
                         bc_higher = 1 ) {
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
    # This model doesn't support exogenous regressors, time to error.
    stop("Exogenous regressors aren't supported by `soothsayer()`")
  },
  .required_specials = c("parameters")
)
#' BATS model
#'
#' @description A \link{fable} wrapper for \link{forecast}[bats]
#' @param formula A TBATS model formula (see details).
#' @param ... Additional arguments (see details).
#' @return A TBATS model, analogous to other model objects within fable/fabletools.
#' @details Accepts and parses several model specials.
#' @note Maybe some other day.
#' @export
BATS <- function(formula, ...) {
  # Create a model class which combines the training method, specials, and data checks
  model_bats <- fabletools::new_model_class("BATS",
                                             # The training method (more on this later)
                                             train = train_bats,
                                             # The formula specials (the next section)
                                             specials = specials_bats,
                                             # Any checks of the unprocessed data, like gaps, ordered, regular, etc.
                                             check = function(.data) {
                                               if (!tsibble::is_regular(.data)) stop("Data must be regular")
                                             }
  )
  # Return a model definition which stores the user's model specification
  fabletools::new_model_definition(model_bats, !!rlang::enquo(formula), ...)
}
