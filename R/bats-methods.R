#' @export
print.BATS <- function (x, ...)
{
  cat(as.character(x[["model_summary"]]))
}
#' @importFrom stats residuals
#' @export
residuals.BATS <- function( object, ... ) {
  object[["resid"]]
}
#' @importFrom stats fitted
#' @export
fitted.BATS <- function( object, ... ) {
  object[["fitted"]]
}
#' @importFrom fabletools forecast
#' @export
forecast.BATS <- function( object, new_data = NULL, specials = NULL, bootstrap = FALSE,
                            times = 5000, ... ) {

  h <- nrow(new_data)
  # set the level to 80 to easily reverse standard error calculation
  fcst <- forecast_bats( object[["fit"]], h = h, level = 80 )
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
refit.BATS <- function( object, new_data, specials = NULL,  ... ) {
  y <- unclass(new_data)[[tsibble::measured_vars(new_data)]]
  model_list <- object[["model_pars"]]

  model <- do.call( tbats, c( list( y = stats::as.ts(y) ),
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
    class = "BATS"
  )
}

#' @importFrom generics glance
#' @export
glance.BATS <- function( x, ... ) {
  
  model_fit <- x[["fit"]]
  resid <- x[["resid"]]
  n_pars <- length(model_fit[["parameters"]][["vect"]])
  dof <- length(resid) - n_pars
  sigma2 <- sum(resid^2)/dof
  loglik <- -x$fit$likelihood/2
  AICc <- -2 * loglik + (2*(n_pars + 1)*(
    length(resid)/(length(resid) - 1 - n_pars)))
  BIC <- -2 * loglik + n_pars * log(length(resid))
  
  tsibble::tibble(sigma2 = sigma2, log_lik = loglik, AIC = x$fit$AIC,
                  AICc = AICc, BIC = BIC, dof = dof)
}
