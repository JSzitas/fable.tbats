
period <- function(x)
{
  n <- length(x)
  spec <- stats::spec.ar(c(x),plot=FALSE)
  if(max(spec$spec)>10) # Arbitrary threshold chosen by trial and error.
  {
    period <- round(1/spec$freq[which.max(spec$spec)])
    if(period==Inf) # Find next local maximum
    {
      j <- which(diff(spec$spec)>0)
      if(length(j)>0)
      {
        nextmax <- j[1] + which.max(spec$spec[j[1]:500])
        period <- round(1/spec$freq[nextmax])
      }
      else
        period <- 1
    }
  }
  else
    period <- 1
  return(period)
}

#' Find seasonalities
#'
#' @description Get seasonalities in a time series by iterative spectral density estimation
#' @param y The time series to detect seasonality in - a numeric vector.
#' @param max_iter The maximal number of iterations - of the spectral density decomposition and
#' aggregation cycles
#' @param aggregator How seasonalities are aggregated - by default **sum**.
#' @param upper_limit The highest possible seasonality to be found in a time series -
#' by default 1500.
#' @details This algorithm computes the spectral density of a time series y, using an AR process
#' @seealso [stats::spec.ar()]. If this returns a period longer than 1, the time series is
#' aggregated, using the previous period as an aggregation window. This is done by applying a
#' function to slices of the time series (by default, the **sum**). Then the first step is repeated
#' with the new, shorter time series. This is repeated until either no seasonality is found,
#' **max_ter** iterations of the algorithm have been carried out, or the **upper_limit** for
#' period length is reached. None that the upper limit is relatively liberal, and should
#' not be reached until you have minute (or smaller) samples of the data, multiple seasonalities,
#' and years of data.
#' @return A vector of seasonalities.
#' @export
find_seasonalities <- function( y, max_iter = 5, aggregator = sum, upper_limit = 1500 ) {

  periods <- list()
  for( iter in seq_len(max_iter) ) {
    last_period <- period(y)
    if( last_period <= 1 ){
      break;
    }
    periods[[iter]] <- last_period
    y <- stats::aggregate(
      stats::ts(y, freq = last_period), # where last_period is last infered periodicity
      nfrequency = 1, # nfrequency always set to 1
      FUN = aggregator # ie mean
    )
  }
  x <- cumprod( unlist(periods))
  x[ x < upper_limit ]
}
