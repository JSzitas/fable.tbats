#' @export
print.TBATS <- function (x, ...)
{
  cat(as.character(x[["model_summary"]]))
}
#' @importFrom stats residuals
#' @export
residuals.TBATS <- function( object, ... ) {
  object[["resid"]]
}
#' @importFrom stats fitted
#' @export
fitted.TBATS <- function( object, ... ) {
  object[["fitted"]]
}
#' @importFrom fabletools forecast
#' @export
forecast.TBATS <- function( object, new_data = NULL, specials = NULL, bootstrap = FALSE,
                            times = 5000, ... ) {

  h <- nrow(new_data)
  # set the level to 80 to easily reverse standard error calculation
  fcst <- forecast::forecast( object[["fit"]], h = h, level = 80 )
  mean_fcst <- c(fcst[["mean"]])
  upper_fcst <- c(fcst[["upper"]])
  # reverse the calculations done in the forecast function to get the standard
  # error
  marg_error <- upper_fcst - mean_fcst
  # here is where the 80 comes in handy - leaving it in to make this more obvious
  st_dev = marg_error / abs(stats::qnorm((100 - 80)/200))

  distributional::dist_normal(  mean_fcst, st_dev )
}
#' @importFrom generics refit
#' @export
refit.TBATS <- function( object, new_data, specials = NULL, reestimate = FALSE,  ... ) {
  y <- unclass(new_data)[[measured_vars(new_data)]]
  model_list <- list( object[["model_pars"]] )

  if( reestimate ) {
    model_list <- c( model_list, list( model = object[["fit"]]) )
  }
  model <- do.call( forecast::tbats, c( list( y = stats::as.ts(y) ),
                               model_list )
  )
  structure(
    list(
      fit = model,
      resid = stats::residuals(model),
      fitted = stats::fitted(model),
      target = tsibble::measured_vars(new_data),
      model_summary = as.character(model),
      model_pars = model_list
    ),
    class = "TBATS"
  )
}
#' @importFrom generics components
#' @export
components.TBATS <- function( object, ... ) {
  comp <- forecast::tbats.components( object[["fit"]] )

  comp <- tsibble::as_tsibble(comp)
  comp <- tidyr::pivot_wider( as.data.frame(comp), names_from = "key")

  box_cox_lambda <- object$fit$lambda
  if (!is.null(box_cox_lambda)) {

    comp_subset <- dplyr::select( comp,
                                  tidyselect::matches("level|observed|season") )
    non_trasnformed <- dplyr::select( comp,
                                      !tidyselect::matches("level|observed|season") )

    comp <- purrr::map_dfc( comp_subset,
                            fabletools::inv_box_cox,
                            lambda = box_cox_lambda )
    comp <- dplyr::bind_cols( non_trasnformed, comp )
    comp <- dplyr::relocate( comp,
                             c("index", "observed", "level"),
                             .before = tidyselect::everything() )

  }

  comp <- dplyr::rename_with( comp, ~ object[["target"]], .cols = "observed")
  comp <- tsibble::as_tsibble(comp, index = "index")

  fabletools::as_dable( comp,
                        #index = "index",
                        response = object[["target"]],
                        # seasons =
                        method = "TBATS")
}
