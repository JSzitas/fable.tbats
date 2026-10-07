#' Simulate future sample paths from a BATS or TBATS model
#'
#' @description Generates sample paths from the fitted innovations state space
#' model, starting from its final state. This is the method behind
#' \code{fabletools::generate()} for models fitted with \code{\link{BATS}} and
#' \code{\link{TBATS}}.
#' @param x A fitted BATS or TBATS model.
#' @param new_data A tsibble of future time points, one path per value of the
#' \code{.rep} key, as prepared by \code{fabletools::generate()}.
#' @param specials Unused; present for compatibility with fabletools.
#' @param ... Unused.
#' @return \code{new_data} with a \code{.sim} column holding the simulated values.
#' @details Innovations live on the model scale, i.e. the Box-Cox transformed
#' scale when the model uses a Box-Cox transformation. Unless \code{new_data}
#' already carries a \code{.innov} column, they are drawn from a normal
#' distribution with the fitted innovation variance. With
#' \code{bootstrap = TRUE}, fabletools fills \code{.innov} by resampling the
#' innovation residuals of the model (see \code{\link{residuals.TBATS}}).
#' Simulated paths are mapped back to the original scale of the series without
#' bias adjustment. Simulation always starts at the end of the training data,
#' whatever the index of \code{new_data}, which mirrors \code{forecast()} for
#' these models.
#' @importFrom generics generate
#' @export
generate.TBATS <- function(x, new_data, specials = NULL, ...) {
  generate_bats_tbats(x, new_data)
}

#' @rdname generate.TBATS
#' @export
generate.BATS <- function(x, new_data, specials = NULL, ...) {
  generate_bats_tbats(x, new_data)
}
