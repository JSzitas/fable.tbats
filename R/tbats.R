train_tbats <- function(.data, specials, ...) {

  # parse arguments to tbats
  parameters <- specials$parameters[[1]]

  if( !is.null( parameters$seasonal_periods )) {
    y <- unclass(.data)[[tsibble::measured_vars(.data)]]
  }
  else {
    y <- stats::as.ts(.data)
  }

  # always set use.parallel to FALSE - since nested parallelism would cause problems
  # and the ONLY situatiion where we avoid that is when someone is running a single
  # TBATS model on a single time series, OR running all models sequentially
  # by setting future::plan("sequential") - so they are probably not too worried
  # about this being slow.
  model <- do.call( forecast::tbats,
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
    class = "TBATS"
  )
}

specials_tbats <- fabletools::new_specials(
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
#' TBATS model
#'
#' @description A \link{fable} wrapper for \link{forecast}[tbats]
#' @param formula A TBATS model formula (see details).
#' @param ... Additional arguments (see details).
#' @return A TBATS model, analogous to other model objects within fable/fabletools.
#' @details Accepts and parses several model specials.
#' @note Maybe some other day.
#' @export
TBATS <- function(formula, ...) {
  # Create a model class which combines the training method, specials, and data checks
  model_tbats <- fabletools::new_model_class("TBATS",
                                             # The training method (more on this later)
                                             train = train_tbats,
                                             # The formula specials (the next section)
                                             specials = specials_tbats,
                                             # Any checks of the unprocessed data, like gaps, ordered, regular, etc.
                                             check = function(.data) {
                                               if (!tsibble::is_regular(.data)) stop("Data must be regular")
                                             }
  )
  # Return a model definition which stores the user's model specification
  fabletools::new_model_definition(model_tbats, !!rlang::enquo(formula), ...)
}
