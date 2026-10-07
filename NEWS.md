# fable.tbats 1.0.0

* The fitting is now implemented in C++ inside this package
  (`src/tbats/tbats.h`, a single header that depends only on the standard
  library and the vendored header-only libraries nlsolver and tinyqr). The
  vendored R and Armadillo code from the forecast package is gone, and so is
  the RcppArmadillo dependency; a C++17 compiler is required.
* Fits are verified against the previous implementation: for fixed
  parameters the state space kernel, seed states, likelihood and forecasts
  reproduce it to tight tolerances, and the model searches reach an equal or
  lower AIC on every test series. Sequential searches run 1.5 to 3 times
  faster than before.
* Forecasts of a Box-Cox model are now a transformed distribution: the
  normal forecast on the transformed scale mapped back through the inverse
  Box-Cox transformation. Quantiles and prediction intervals are exact (and
  asymmetric), `.mean` is the second-order bias-adjusted mean that
  distributional computes, and the median is the back-transformed point
  forecast. Previously a normal distribution was fitted to the 80 percent
  bound on the original scale. `bias_adj` now affects fitted values only.
* `components()` is available for both models.
* `refit()` with the default `reestimate = FALSE` keeps the estimated
  parameters and re-estimates the seed states for the new series;
  `reestimate = TRUE` repeats the model search. Previously the search was
  always repeated.
* Missing values are rejected with a clear message rather than silently
  trimmed to the longest complete stretch.
* Found and fixed on the way: the vendored forecast C++ refreshed the
  seasonal smoothing value of the first seasonal period only during
  optimisation, so a BATS model with two or more seasonal periods and ARMA
  errors was optimised against the wrong objective. The port is not affected.
  Also noted: the forecast package builds the gain vector for TBATS
  prediction intervals without the trend coefficient; the port's intervals
  include it.

# fable.tbats 0.5.0

* `generate()` is now supported for BATS and TBATS models, producing simulated
  future sample paths from the fitted state space model. Both parametric
  (normal innovations with the fitted variance) and bootstrap (resampled
  innovation residuals) simulation work through `fabletools::generate()`.
* `residuals()` gained a `type` argument. The default `"innovation"` returns the
  one step ahead errors of the state space model, which are on the Box-Cox
  transformed scale when the model uses a Box-Cox transformation. The previous
  behaviour, the difference between observed and fitted values, is available
  with `type = "response"`. `augment()` now reports both correctly.
