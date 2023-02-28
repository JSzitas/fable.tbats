# nocov start
msts <- function (data, seasonal.periods, ts.frequency = floor(max(seasonal.periods)), 
                  ...) 
{
  if (inherits(data, "ts") && frequency(data) == ts.frequency && 
      length(list(...)) == 0) {
    object <- data
  }
  else {
    object <- ts(data = data, frequency = ts.frequency, ...)
  }
  if (length(seasonal.periods) > 1L) {
    class(object) <- c("msts", "ts")
    attr(object, "msts") <- sort(seasonal.periods)
  }
  return(object)
}



forecast_tbats <- function (object, h, level = c(80, 95), fan = FALSE, biasadj = NULL, 
                            ...) 
{
  if (identical(class(object), "bats")) {
    return(forecast_bats(object, h, level, fan, biasadj, 
                         ...))
  }
  if (any(class(object$y) == "ts")) {
    ts.frequency <- frequency(object$y)
  }
  else {
    ts.frequency <- ifelse(!is.null(object$seasonal.periods), 
                           max(object$seasonal.periods), 1)
  }
  if (missing(h)) {
    if (is.null(object$seasonal.periods)) {
      h <- ifelse(ts.frequency == 1, 10, 2 * ts.frequency)
    }
    else {
      h <- 2 * max(object$seasonal.periods)
    }
  }
  else if (h <= 0) {
    stop("Forecast horizon out of bounds")
  }
  if (fan) {
    level <- seq(51, 99, by = 3)
  }
  else {
    if (min(level) > 0 && max(level) < 1) {
      level <- 100 * level
    }
    else if (min(level) < 0 || max(level) > 99.99) {
      stop("Confidence limit out of range")
    }
  }
  if (!is.null(object$k.vector)) {
    tau <- 2 * sum(object$k.vector)
  }
  else {
    tau <- 0
  }
  x <- matrix(0, nrow = nrow(object$x), ncol = h)
  y.forecast <- numeric(h)
  if (!is.null(object$beta)) {
    adj.beta <- 1
  }
  else {
    adj.beta <- 0
  }
  w <- .Call("makeTBATSWMatrix", smallPhi_s = object$damping.parameter, 
             kVector_s = as.integer(object$k.vector), arCoefs_s = object$ar.coefficients, 
             maCoefs_s = object$ma.coefficients, tau_s = as.integer(tau), 
             PACKAGE = "fable.tbats")
  if (!is.null(object$seasonal.periods)) {
    gamma.bold <- matrix(0, nrow = 1, ncol = tau)
    .Call("updateTBATSGammaBold", gammaBold_s = gamma.bold, 
          kVector_s = as.integer(object$k.vector), gammaOne_s = object$gamma.one.values, 
          gammaTwo_s = object$gamma.two.values, PACKAGE = "fable.tbats")
  }
  else {
    gamma.bold <- NULL
  }
  g <- matrix(0, nrow = (tau + 1 + adj.beta + object$p + object$q), 
              ncol = 1)
  if (object$p != 0) {
    g[(1 + adj.beta + tau + 1), 1] <- 1
  }
  if (object$q != 0) {
    g[(1 + adj.beta + tau + object$p + 1), 1] <- 1
  }
  .Call("updateTBATSGMatrix", g_s = g, gammaBold_s = gamma.bold, 
        alpha_s = object$alpha, beta_s = object$beta.v, PACKAGE = "fable.tbats")
  F <- makeTBATSFMatrix(alpha = object$alpha, beta = object$beta, 
                        small.phi = object$damping.parameter, seasonal.periods = object$seasonal.periods, 
                        k.vector = as.integer(object$k.vector), gamma.bold.matrix = gamma.bold, 
                        ar.coefs = object$ar.coefficients, ma.coefs = object$ma.coefficients)
  y.forecast[1] <- w$w.transpose %*% object$x[, ncol(object$x)]
  x[, 1] <- F %*% object$x[, ncol(object$x)]
  if (h > 1) {
    for (t in 2:h) {
      x[, t] <- F %*% x[, (t - 1)]
      y.forecast[t] <- w$w.transpose %*% x[, (t - 1)]
    }
  }
  lower.bounds <- upper.bounds <- matrix(NA, ncol = length(level), 
                                         nrow = h)
  variance.multiplier <- numeric(h)
  variance.multiplier[1] <- 1
  if (h > 1) {
    for (j in 1:(h - 1)) {
      if (j == 1) {
        f.running <- diag(ncol(F))
      }
      else {
        f.running <- f.running %*% F
      }
      c.j <- w$w.transpose %*% f.running %*% g
      variance.multiplier[(j + 1)] <- variance.multiplier[j] + 
        c.j^2
    }
  }
  variance <- object$variance * variance.multiplier
  st.dev <- sqrt(variance)
  for (i in 1:length(level)) {
    marg.error <- st.dev * abs(qnorm((100 - level[i])/200))
    lower.bounds[, i] <- y.forecast - marg.error
    upper.bounds[, i] <- y.forecast + marg.error
  }
  if (!is.null(object$lambda)) {
    y.forecast <- InvBoxCox(y.forecast, object$lambda, biasadj, 
                            list(level = level, upper = upper.bounds, lower = lower.bounds))
    lower.bounds <- InvBoxCox(lower.bounds, object$lambda)
    if (object$lambda < 1) {
      lower.bounds <- pmax(lower.bounds, 0)
    }
    upper.bounds <- InvBoxCox(upper.bounds, object$lambda)
  }
  start.time <- start(object$y)
  y <- ts(c(object$y, 0), start = start.time, frequency = ts.frequency)
  fcast.start.time <- end(y)
  x <- msts(object$y, seasonal.periods = (if (!is.null(object$seasonal.periods)) {
    object$seasonal.periods
  }
  else {
    ts.frequency
  }), ts.frequency = ts.frequency, start = start.time)
  fitted.values <- msts(object$fitted.values, seasonal.periods = (if (!is.null(object$seasonal.periods)) {
    object$seasonal.periods
  }
  else {
    ts.frequency
  }), ts.frequency = ts.frequency, start = start.time)
  y.forecast <- msts(y.forecast, seasonal.periods = (if (!is.null(object$seasonal.periods)) {
    object$seasonal.periods
  }
  else {
    ts.frequency
  }), ts.frequency = ts.frequency, start = fcast.start.time)
  upper.bounds <- msts(upper.bounds, seasonal.periods = (if (!is.null(object$seasonal.periods)) {
    object$seasonal.periods
  }
  else {
    ts.frequency
  }), ts.frequency = ts.frequency, start = fcast.start.time)
  lower.bounds <- msts(lower.bounds, seasonal.periods = (if (!is.null(object$seasonal.periods)) {
    object$seasonal.periods
  }
  else {
    ts.frequency
  }), ts.frequency = ts.frequency, start = fcast.start.time)
  colnames(upper.bounds) <- colnames(lower.bounds) <- paste0(level, 
                                                             "%")
  forecast.object <- list(model = object, mean = y.forecast, 
                          level = level, x = x, series = object$series, upper = upper.bounds, 
                          lower = lower.bounds, fitted = fitted.values, method = as.character(object), 
                          residuals = object$errors)
  if (is.null(object$series)) {
    forecast.object$series <- deparse(object$call$y)
  }
  class(forecast.object) <- "forecast"
  return(forecast.object)
}

forecast_bats <- function (object, h, level = c(80, 95), fan = FALSE, biasadj = NULL, 
                           ...) 
{
  if (any(class(object$y) == "ts")) {
    ts.frequency <- frequency(object$y)
  }
  else {
    ts.frequency <- ifelse(!is.null(object$seasonal.periods), 
                           max(object$seasonal.periods), 1)
  }
  if (missing(h)) {
    if (is.null(object$seasonal.periods)) {
      h <- ifelse(ts.frequency == 1, 10, 2 * ts.frequency)
    }
    else {
      h <- 2 * max(object$seasonal.periods)
    }
  }
  else if (h <= 0) {
    stop("Forecast horizon out of bounds")
  }
  if (fan) {
    level <- seq(51, 99, by = 3)
  }
  else {
    if (min(level) > 0 && max(level) < 1) {
      level <- 100 * level
    }
    else if (min(level) < 0 || max(level) > 99.99) {
      stop("Confidence limit out of range")
    }
  }
  x <- matrix(0, nrow = nrow(object$x), ncol = h)
  y.forecast <- numeric(h)
  w <- .Call("makeBATSWMatrix", smallPhi_s = object$damping.parameter, 
             sPeriods_s = object$seasonal.periods, arCoefs_s = object$ar.coefficients, 
             maCoefs_s = object$ma.coefficients, PACKAGE = "fable.tbats")
  g <- .Call("makeBATSGMatrix", object$alpha, object$beta, 
             object$gamma.values, object$seasonal.periods, length(object$ar.coefficients), 
             length(object$ma.coefficients), PACKAGE = "fable.tbats")
  F <- makeFMatrix(alpha = object$alpha, beta = object$beta, 
                   small.phi = object$damping.parameter, seasonal.periods = object$seasonal.periods, 
                   gamma.bold.matrix = g$gamma.bold.matrix, ar.coefs = object$ar.coefficients, 
                   ma.coefs = object$ma.coefficients)
  y.forecast[1] <- w$w.transpose %*% object$x[, ncol(object$x)]
  x[, 1] <- F %*% object$x[, ncol(object$x)]
  if (h > 1) {
    for (t in 2:h) {
      x[, t] <- F %*% x[, (t - 1)]
      y.forecast[t] <- w$w.transpose %*% x[, (t - 1)]
    }
  }
  lower.bounds <- upper.bounds <- matrix(NA, ncol = length(level), 
                                         nrow = h)
  variance.multiplier <- numeric(h)
  variance.multiplier[1] <- 1
  if (h > 1) {
    for (j in 1:(h - 1)) {
      if (j == 1) {
        f.running <- diag(ncol(F))
      }
      else {
        f.running <- f.running %*% F
      }
      c.j <- w$w.transpose %*% f.running %*% g$g
      variance.multiplier[(j + 1)] <- variance.multiplier[j] + 
        c.j^2
    }
  }
  variance <- object$variance * variance.multiplier
  st.dev <- sqrt(variance)
  for (i in 1:length(level)) {
    marg.error <- st.dev * abs(qnorm((100 - level[i])/200))
    lower.bounds[, i] <- y.forecast - marg.error
    upper.bounds[, i] <- y.forecast + marg.error
  }
  if (!is.null(object$lambda)) {
    y.forecast <- InvBoxCox(y.forecast, object$lambda, biasadj, 
                            list(level = level, upper = upper.bounds, lower = lower.bounds))
    lower.bounds <- InvBoxCox(lower.bounds, object$lambda)
    if (object$lambda < 1) {
      lower.bounds <- pmax(lower.bounds, 0)
    }
    upper.bounds <- InvBoxCox(upper.bounds, object$lambda)
  }
  start.time <- start(object$y)
  y <- ts(c(object$y, 0), start = start.time, frequency = ts.frequency)
  fcast.start.time <- end(y)
  x <- msts(object$y, seasonal.periods = (if (!is.null(object$seasonal.periods)) {
    object$seasonal.periods
  }
  else {
    ts.frequency
  }), ts.frequency = ts.frequency, start = start.time)
  fitted.values <- msts(object$fitted.values, seasonal.periods = (if (!is.null(object$seasonal.periods)) {
    object$seasonal.periods
  }
  else {
    ts.frequency
  }), ts.frequency = ts.frequency, start = start.time)
  y.forecast <- msts(y.forecast, seasonal.periods = (if (!is.null(object$seasonal.periods)) {
    object$seasonal.periods
  }
  else {
    ts.frequency
  }), ts.frequency = ts.frequency, start = fcast.start.time)
  upper.bounds <- msts(upper.bounds, seasonal.periods = (if (!is.null(object$seasonal.periods)) {
    object$seasonal.periods
  }
  else {
    ts.frequency
  }), ts.frequency = ts.frequency, start = fcast.start.time)
  lower.bounds <- msts(lower.bounds, seasonal.periods = (if (!is.null(object$seasonal.periods)) {
    object$seasonal.periods
  }
  else {
    ts.frequency
  }), ts.frequency = ts.frequency, start = fcast.start.time)
  colnames(upper.bounds) <- colnames(lower.bounds) <- paste0(level, 
                                                             "%")
  forecast.object <- list(model = object, mean = y.forecast, 
                          level = level, x = x, series = object$series, upper = upper.bounds, 
                          lower = lower.bounds, fitted = fitted.values, method = as.character(object), 
                          residuals = object$errors)
  if (is.null(object$series)) {
    forecast.object$series <- deparse(object$call$y)
  }
  class(forecast.object) <- "forecast"
  return(forecast.object)
}

fitSpecificBATS <- function (y, use.box.cox, use.beta, use.damping, seasonal.periods = NULL, 
                             starting.params = NULL, x.nought = NULL, ar.coefs = NULL, 
                             ma.coefs = NULL, init.box.cox = NULL, bc.lower = 0, bc.upper = 1, 
                             biasadj = FALSE) 
{
  if (!is.null(seasonal.periods)) {
    seasonal.periods <- as.integer(sort(seasonal.periods))
  }
  if (is.null(starting.params)) {
    if (!is.null(ar.coefs)) {
      p <- length(ar.coefs)
    }
    else {
      p <- 0
    }
    if (!is.null(ma.coefs)) {
      q <- length(ma.coefs)
    }
    else {
      q <- 0
    }
    if (sum(seasonal.periods) > 16) {
      alpha <- (1e-06)
    }
    else {
      alpha <- 0.09
    }
    if (use.beta) {
      if (sum(seasonal.periods) > 16) {
        beta.v <- (5e-07)
      }
      else {
        beta.v <- 0.05
      }
      b <- 0
      if (use.damping) {
        small.phi <- 0.999
      }
      else {
        small.phi <- 1
      }
    }
    else {
      beta.v <- NULL
      b <- NULL
      small.phi <- NULL
      use.damping <- FALSE
    }
    if (!is.null(seasonal.periods)) {
      gamma.v <- rep(0.001, length(seasonal.periods))
      s.vector <- numeric(sum(seasonal.periods))
    }
    else {
      gamma.v <- NULL
      s.vector <- NULL
    }
    if (use.box.cox) {
      if (!is.null(init.box.cox)) {
        lambda <- init.box.cox
      }
      else {
        lambda <- BoxCox.lambda(y, lower = 0, upper = 1.5)
      }
      y.transformed <- BoxCox(y, lambda = lambda)
      lambda <- attr(y.transformed, "lambda")
    }
    else {
      lambda <- NULL
    }
  }
  else {
    paramz <- unParameterise(starting.params$vect, starting.params$control)
    lambda <- paramz$lambda
    alpha <- paramz$alpha
    beta.v <- paramz$beta
    b <- 0
    small.phi <- paramz$small.phi
    gamma.v <- paramz$gamma.v
    if (!is.null(seasonal.periods)) {
      s.vector <- numeric(sum(seasonal.periods))
    }
    else {
      s.vector <- NULL
    }
    if (!is.null(ar.coefs)) {
      p <- length(ar.coefs)
    }
    else {
      p <- 0
    }
    if (!is.null(ma.coefs)) {
      q <- length(ma.coefs)
    }
    else {
      q <- 0
    }
  }
  if (is.null(x.nought)) {
    if (!is.null(ar.coefs)) {
      d.vector <- numeric(length(ar.coefs))
    }
    else {
      d.vector <- NULL
    }
    if (!is.null(ma.coefs)) {
      epsilon.vector <- numeric(length(ma.coefs))
    }
    else {
      epsilon.vector <- NULL
    }
    x.nought <- makeXMatrix(l = 0, b = b, s.vector = s.vector, 
                            d.vector = d.vector, epsilon.vector = epsilon.vector)$x
  }
  param.vector <- parameterise(alpha = alpha, beta.v = beta.v, 
                               small.phi = small.phi, gamma.v = gamma.v, lambda = lambda, 
                               ar.coefs = ar.coefs, ma.coefs = ma.coefs)
  par.scale <- makeParscaleBATS(param.vector$control)
  w <- .Call("makeBATSWMatrix", smallPhi_s = small.phi, sPeriods_s = seasonal.periods, 
             arCoefs_s = ar.coefs, maCoefs_s = ma.coefs, PACKAGE = "fable.tbats")
  g <- .Call("makeBATSGMatrix", as.numeric(alpha), beta.v, 
             gamma.v, seasonal.periods, as.integer(p), as.integer(q), 
             PACKAGE = "fable.tbats")
  F <- makeFMatrix(alpha = alpha, beta = beta.v, small.phi = small.phi, 
                   seasonal.periods = seasonal.periods, gamma.bold.matrix = g$gamma.bold.matrix, 
                   ar.coefs = ar.coefs, ma.coefs = ma.coefs)
  D <- F - g$g %*% w$w.transpose
  if (use.box.cox) {
    y.transformed <- BoxCox(y, lambda = lambda)
    lambda <- attr(y.transformed, "lambda")
    y.tilda <- calcModel(y.transformed, x.nought, F, g$g, 
                         w)$e
  }
  else {
    y.tilda <- calcModel(y, x.nought, F, g$g, w)$e
  }
  w.tilda.transpose <- matrix(0, nrow = length(y), ncol = ncol(w$w.transpose))
  w.tilda.transpose[1, ] <- w$w.transpose
  w.tilda.transpose <- .Call("calcWTilda", wTildaTransposes = w.tilda.transpose, 
                             Ds = D, PACKAGE = "fable.tbats")
  if (!is.null(seasonal.periods)) {
    list.cut.w <- cutW(use.beta = use.beta, w.tilda.transpose = w.tilda.transpose, 
                       seasonal.periods = seasonal.periods, p = p, q = q)
    w.tilda.transpose <- list.cut.w$matrix
    mask.vector <- list.cut.w$mask.vector
    coefs <- lm(t(y.tilda) ~ w.tilda.transpose - 1)$coefficients
    x.nought <- calcSeasonalSeeds(use.beta = use.beta, coefs = coefs, 
                                  seasonal.periods = seasonal.periods, mask.vector = mask.vector, 
                                  p = p, q = q)
  }
  else {
    if ((p != 0) | (q != 0)) {
      end.cut <- ncol(w.tilda.transpose)
      start.cut <- end.cut - (p + q) + 1
      w.tilda.transpose <- w.tilda.transpose[, -c(start.cut:end.cut)]
    }
    x.nought <- lm(t(y.tilda) ~ w.tilda.transpose - 1)$coefficients
    x.nought <- matrix(x.nought, nrow = length(x.nought), 
                       ncol = 1)
    if ((p != 0) | (q != 0)) {
      arma.seed.states <- numeric((p + q))
      arma.seed.states <- matrix(arma.seed.states, nrow = length(arma.seed.states), 
                                 ncol = 1)
      x.nought <- rbind(x.nought, arma.seed.states)
    }
  }
  opt.env <- new.env()
  assign("F", F, envir = opt.env)
  assign("w.transpose", w$w.transpose, envir = opt.env)
  assign("g", g$g, envir = opt.env)
  assign("gamma.bold.matrix", g$gamma.bold.matrix, envir = opt.env)
  assign("y", matrix(y, nrow = 1, ncol = length(y)), envir = opt.env)
  assign("y.hat", matrix(0, nrow = 1, ncol = length(y)), envir = opt.env)
  assign("e", matrix(0, nrow = 1, ncol = length(y)), envir = opt.env)
  assign("x", matrix(0, nrow = length(x.nought), ncol = length(y)), 
         envir = opt.env)
  if (!is.null(seasonal.periods)) {
    tau <- sum(seasonal.periods)
  }
  else {
    tau <- 0
  }
  if (use.box.cox) {
    assign("x.nought.untransformed", InvBoxCox(x.nought, 
                                               lambda = lambda), envir = opt.env)
    optim.like <- optim(par = param.vector$vect, fn = calcLikelihood, 
                        method = "Nelder-Mead", opt.env = opt.env, use.beta = use.beta, 
                        use.small.phi = use.damping, seasonal.periods = seasonal.periods, 
                        p = p, q = q, tau = tau, bc.lower = bc.lower, bc.upper = bc.upper, 
                        control = list(maxit = (100 * length(param.vector$vect)^2), 
                                       parscale = par.scale))
    paramz <- unParameterise(optim.like$par, param.vector$control)
    lambda <- paramz$lambda
    alpha <- paramz$alpha
    beta.v <- paramz$beta
    small.phi <- paramz$small.phi
    gamma.v <- paramz$gamma.v
    ar.coefs <- paramz$ar.coefs
    ma.coefs <- paramz$ma.coefs
    x.nought <- BoxCox(opt.env$x.nought.untransformed, lambda = lambda)
    lambda <- attr(x.nought, "lambda")
    w <- .Call("makeBATSWMatrix", smallPhi_s = small.phi, 
               sPeriods_s = seasonal.periods, arCoefs_s = ar.coefs, 
               maCoefs_s = ma.coefs, PACKAGE = "fable.tbats")
    g <- .Call("makeBATSGMatrix", as.numeric(alpha), beta.v, 
               gamma.v, seasonal.periods, as.integer(p), as.integer(q), 
               PACKAGE = "fable.tbats")
    F <- makeFMatrix(alpha = alpha, beta = beta.v, small.phi = small.phi, 
                     seasonal.periods = seasonal.periods, gamma.bold.matrix = g$gamma.bold.matrix, 
                     ar.coefs = ar.coefs, ma.coefs = ma.coefs)
    y.transformed <- BoxCox(y, lambda = lambda)
    lambda <- attr(y.transformed, "lambda")
    fitted.values.and.errors <- calcModel(y.transformed, 
                                          x.nought, F, g$g, w)
    e <- fitted.values.and.errors$e
    variance <- sum((e * e))/length(y)
    fitted.values <- InvBoxCox(fitted.values.and.errors$y.hat, 
                               lambda = lambda, biasadj, variance)
    attr(lambda, "biasadj") <- biasadj
  }
  else {
    if (length(param.vector$vect) > 1) {
      optim.like <- optim(par = param.vector$vect, fn = calcLikelihoodNOTransformed, 
                          method = "Nelder-Mead", opt.env = opt.env, x.nought = x.nought, 
                          use.beta = use.beta, use.small.phi = use.damping, 
                          seasonal.periods = seasonal.periods, p = p, q = q, 
                          tau = tau, control = list(maxit = (100 * length(param.vector$vect)^2), 
                                                    parscale = par.scale))
    }
    else {
      optim.like <- optim(par = param.vector$vect, fn = calcLikelihoodNOTransformed, 
                          method = "BFGS", opt.env = opt.env, x.nought = x.nought, 
                          use.beta = use.beta, use.small.phi = use.damping, 
                          seasonal.periods = seasonal.periods, p = p, q = q, 
                          tau = tau, control = list(parscale = par.scale))
    }
    paramz <- unParameterise(optim.like$par, param.vector$control)
    lambda <- paramz$lambda
    alpha <- paramz$alpha
    beta.v <- paramz$beta
    small.phi <- paramz$small.phi
    gamma.v <- paramz$gamma.v
    ar.coefs <- paramz$ar.coefs
    ma.coefs <- paramz$ma.coefs
    w <- .Call("makeBATSWMatrix", smallPhi_s = small.phi, 
               sPeriods_s = seasonal.periods, arCoefs_s = ar.coefs, 
               maCoefs_s = ma.coefs, PACKAGE = "fable.tbats")
    g <- .Call("makeBATSGMatrix", as.numeric(alpha), beta.v, 
               gamma.v, seasonal.periods, as.integer(p), as.integer(q), 
               PACKAGE = "fable.tbats")
    F <- makeFMatrix(alpha = alpha, beta = beta.v, small.phi <- small.phi, 
                     seasonal.periods = seasonal.periods, gamma.bold.matrix = g$gamma.bold.matrix, 
                     ar.coefs = ar.coefs, ma.coefs = ma.coefs)
    fitted.values.and.errors <- calcModel(y, x.nought, F, 
                                          g$g, w)
    e <- fitted.values.and.errors$e
    fitted.values <- fitted.values.and.errors$y.hat
    variance <- sum((e * e))/length(y)
  }
  likelihood <- optim.like$value
  aic <- likelihood + 2 * (length(param.vector$vect) + nrow(x.nought))
  model.for.output <- list(lambda = lambda, alpha = alpha, 
                           beta = beta.v, damping.parameter = small.phi, gamma.values = gamma.v, 
                           ar.coefficients = ar.coefs, ma.coefficients = ma.coefs, 
                           likelihood = likelihood, optim.return.code = optim.like$convergence, 
                           variance = variance, AIC = aic, parameters = list(vect = optim.like$par, 
                                                                             control = param.vector$control), seed.states = x.nought, 
                           fitted.values = c(fitted.values), errors = c(e), x = fitted.values.and.errors$x, 
                           seasonal.periods = seasonal.periods, y = y)
  class(model.for.output) <- "bats"
  return(model.for.output)
}




filterSpecifics <- function (y, box.cox, trend, damping, seasonal.periods, use.arma.errors, 
                             force.seasonality = FALSE, init.box.cox = NULL, bc.lower = 0, 
                             bc.upper = 1, biasadj = FALSE, ...) 
{
  if (!trend && damping) {
    return(list(AIC = Inf))
  }
  first.model <- fitSpecificBATS(y, use.box.cox = box.cox, 
                                 use.beta = trend, use.damping = damping, seasonal.periods = seasonal.periods, 
                                 init.box.cox = init.box.cox, bc.lower = bc.lower, bc.upper = bc.upper, 
                                 biasadj = biasadj)
  if (!is.null(seasonal.periods) && !force.seasonality) {
    non.seasonal.model <- fitSpecificBATS(y, use.box.cox = box.cox, 
                                          use.beta = trend, use.damping = damping, seasonal.periods = NULL, 
                                          init.box.cox = init.box.cox, bc.lower = bc.lower, 
                                          bc.upper = bc.upper, biasadj = biasadj)
    if (first.model$AIC > non.seasonal.model$AIC) {
      seasonal.periods <- NULL
      first.model <- non.seasonal.model
    }
  }
  if (use.arma.errors) {
    suppressWarnings(arma <- auto_arma(as.numeric(first.model$errors), 
                                       d = 0, ...))
    p <- arma$arma[1]
    q <- arma$arma[2]
    if (p != 0 || q != 0) {
      if (p != 0) {
        ar.coefs <- numeric(p)
      }
      else {
        ar.coefs <- NULL
      }
      if (q != 0) {
        ma.coefs <- numeric(q)
      }
      else {
        ma.coefs <- NULL
      }
      starting.params <- first.model$parameters
      second.model <- fitSpecificBATS(y, use.box.cox = box.cox, 
                                      use.beta = trend, use.damping = damping, seasonal.periods = seasonal.periods, 
                                      ar.coefs = ar.coefs, ma.coefs = ma.coefs, init.box.cox = init.box.cox, 
                                      bc.lower = bc.lower, bc.upper = bc.upper, biasadj = biasadj)
      if (second.model$AIC < first.model$AIC) {
        return(second.model)
      }
      else {
        return(first.model)
      }
    }
    else {
      return(first.model)
    }
  }
  else {
    return(first.model)
  }
}

fitSpecificTBATS <- function (y, use.box.cox, use.beta, use.damping, seasonal.periods = NULL, 
                              k.vector = NULL, starting.params = NULL, x.nought = NULL, 
                              ar.coefs = NULL, ma.coefs = NULL, init.box.cox = NULL, bc.lower = 0, 
                              bc.upper = 1, biasadj = FALSE) 
{
  if (!is.null(seasonal.periods)) {
    seasonal.periods <- sort(seasonal.periods)
  }
  if (is.null(starting.params)) {
    if (!is.null(ar.coefs)) {
      p <- length(ar.coefs)
    }
    else {
      p <- 0
    }
    if (!is.null(ma.coefs)) {
      q <- length(ma.coefs)
    }
    else {
      q <- 0
    }
    alpha <- 0.09
    if (use.beta) {
      adj.beta <- 1
      beta.v <- 0.05
      b <- 0
      if (use.damping) {
        small.phi <- 0.999
      }
      else {
        small.phi <- 1
      }
    }
    else {
      adj.beta <- 0
      beta.v <- NULL
      b <- NULL
      small.phi <- NULL
      use.damping <- FALSE
    }
    if (!is.null(seasonal.periods)) {
      gamma.one.v <- rep(0, length(k.vector))
      gamma.two.v <- rep(0, length(k.vector))
      s.vector <- numeric(2 * sum(k.vector))
      k.vector <- as.integer(k.vector)
    }
    else {
      gamma.one.v <- NULL
      gamma.two.v <- NULL
      s.vector <- NULL
    }
    if (use.box.cox) {
      if (!is.null(init.box.cox)) {
        lambda <- init.box.cox
      }
      else {
        lambda <- BoxCox.lambda(y, lower = 0, upper = 1.5)
      }
      y.transformed <- BoxCox(y, lambda = lambda)
      lambda <- attr(y.transformed, "lambda")
    }
    else {
      lambda <- NULL
    }
  }
  else {
    paramz <- unParameteriseTBATS(starting.params$vect, starting.params$control)
    lambda <- paramz$lambda
    alpha <- paramz$alpha
    beta.v <- paramz$beta
    if (!is.null(beta.v)) {
      adj.beta <- 1
    }
    else {
      adj.beta <- 0
    }
    b <- 0
    small.phi <- paramz$small.phi
    gamma.one.v <- paramz$gamma.one.v
    gamma.two.v <- paramz$gamma.two.v
    if (!is.null(seasonal.periods)) {
      s.vector <- numeric(2 * sum(k.vector))
    }
    else {
      s.vector <- NULL
    }
    if (!is.null(ar.coefs)) {
      p <- length(ar.coefs)
    }
    else {
      p <- 0
    }
    if (!is.null(ma.coefs)) {
      q <- length(ma.coefs)
    }
    else {
      q <- 0
    }
  }
  if (is.null(x.nought)) {
    if (!is.null(ar.coefs)) {
      d.vector <- numeric(length(ar.coefs))
    }
    else {
      d.vector <- NULL
    }
    if (!is.null(ma.coefs)) {
      epsilon.vector <- numeric(length(ma.coefs))
    }
    else {
      epsilon.vector <- NULL
    }
    x.nought <- makeXMatrix(l = 0, b = b, s.vector = s.vector, 
                            d.vector = d.vector, epsilon.vector = epsilon.vector)$x
  }
  param.vector <- parameterise(alpha = alpha, beta.v = beta.v, 
                               small.phi = small.phi, gamma.v = cbind(gamma.one.v, gamma.two.v), 
                               lambda = lambda, ar.coefs = ar.coefs, ma.coefs = ma.coefs)
  par.scale <- makeParscale(param.vector$control)
  if (!is.null(seasonal.periods)) {
    tau <- as.integer(2 * sum(k.vector))
  }
  else {
    tau <- as.integer(0)
  }
  w <- .Call("makeTBATSWMatrix", smallPhi_s = small.phi, kVector_s = k.vector, 
             arCoefs_s = ar.coefs, maCoefs_s = ma.coefs, tau_s = tau, 
             PACKAGE = "fable.tbats")
  if (!is.null(seasonal.periods)) {
    gamma.bold <- matrix(0, nrow = 1, ncol = (2 * sum(k.vector)))
    .Call("updateTBATSGammaBold", gammaBold_s = gamma.bold, 
          kVector_s = k.vector, gammaOne_s = gamma.one.v, gammaTwo_s = gamma.two.v, 
          PACKAGE = "fable.tbats")
  }
  else {
    gamma.bold <- NULL
  }
  g <- matrix(0, nrow = ((2 * sum(k.vector)) + 1 + adj.beta + 
                           p + q), ncol = 1)
  if (p != 0) {
    g[(1 + adj.beta + tau + 1), 1] <- 1
  }
  if (q != 0) {
    g[(1 + adj.beta + tau + p + 1), 1] <- 1
  }
  .Call("updateTBATSGMatrix", g_s = g, gammaBold_s = gamma.bold, 
        alpha_s = alpha, beta_s = beta.v, PACKAGE = "fable.tbats")
  F <- makeTBATSFMatrix(alpha = alpha, beta = beta.v, small.phi = small.phi, 
                        seasonal.periods = seasonal.periods, k.vector = k.vector, 
                        gamma.bold.matrix = gamma.bold, ar.coefs = ar.coefs, 
                        ma.coefs = ma.coefs)
  D <- F - g %*% w$w.transpose
  opt.env <- new.env()
  assign("F", F, envir = opt.env)
  assign("w.transpose", w$w.transpose, envir = opt.env)
  assign("g", g, envir = opt.env)
  assign("gamma.bold", gamma.bold, envir = opt.env)
  assign("k.vector", k.vector, envir = opt.env)
  assign("y", matrix(y, nrow = 1, ncol = length(y)), envir = opt.env)
  assign("y.hat", matrix(0, nrow = 1, ncol = length(y)), envir = opt.env)
  assign("e", matrix(0, nrow = 1, ncol = length(y)), envir = opt.env)
  assign("x", matrix(0, nrow = length(x.nought), ncol = length(y)), 
         envir = opt.env)
  if (use.box.cox) {
    y.transformed <- BoxCox(y, lambda = lambda)
    lambda <- attr(y.transformed, "lambda")
    .Call("calcTBATSFaster", ys = matrix(y.transformed, nrow = 1, 
                                         ncol = length(y.transformed)), yHats = opt.env$y.hat, 
          wTransposes = opt.env$w.transpose, Fs = opt.env$F, 
          xs = opt.env$x, gs = opt.env$g, es = opt.env$e, xNought_s = x.nought, 
          PACKAGE = "fable.tbats")
    y.tilda <- opt.env$e
  }
  else {
    .Call("calcTBATSFaster", ys = opt.env$y, yHats = opt.env$y.hat, 
          wTransposes = opt.env$w.transpose, Fs = opt.env$F, 
          xs = opt.env$x, gs = opt.env$g, es = opt.env$e, xNought_s = x.nought, 
          PACKAGE = "fable.tbats")
    y.tilda <- opt.env$e
  }
  w.tilda.transpose <- matrix(0, nrow = length(y), ncol = ncol(w$w.transpose))
  w.tilda.transpose[1, ] <- w$w.transpose
  w.tilda.transpose <- .Call("calcWTilda", wTildaTransposes = w.tilda.transpose, 
                             Ds = D, PACKAGE = "fable.tbats")
  if ((p != 0) | (q != 0)) {
    end.cut <- ncol(w.tilda.transpose)
    start.cut <- end.cut - (p + q) + 1
    w.tilda.transpose <- w.tilda.transpose[, -c(start.cut:end.cut)]
  }
  x.nought <- lm(t(y.tilda) ~ w.tilda.transpose - 1)$coefficients
  x.nought <- matrix(x.nought, nrow = length(x.nought), ncol = 1)
  if ((p != 0) | (q != 0)) {
    arma.seed.states <- numeric((p + q))
    arma.seed.states <- matrix(arma.seed.states, nrow = length(arma.seed.states), 
                               ncol = 1)
    x.nought <- rbind(x.nought, arma.seed.states)
  }
  if (use.box.cox) {
    assign("x.nought.untransformed", InvBoxCox(x.nought, 
                                               lambda = lambda), envir = opt.env)
    optim.like <- optim(par = param.vector$vect, fn = calcLikelihoodTBATS, 
                        method = "Nelder-Mead", opt.env = opt.env, use.beta = use.beta, 
                        use.small.phi = use.damping, seasonal.periods = seasonal.periods, 
                        param.control = param.vector$control, p = p, q = q, 
                        tau = tau, bc.lower = bc.lower, bc.upper = bc.upper, 
                        control = list(maxit = (100 * length(param.vector$vect)^2), 
                                       parscale = par.scale))
    paramz <- unParameteriseTBATS(optim.like$par, param.vector$control)
    lambda <- paramz$lambda
    alpha <- paramz$alpha
    beta.v <- paramz$beta
    small.phi <- paramz$small.phi
    gamma.one.v <- paramz$gamma.one.v
    gamma.two.v <- paramz$gamma.two.v
    if (!is.null(paramz$ar.coefs)) {
      p <- length(paramz$ar.coefs)
      ar.coefs <- matrix(paramz$ar.coefs, nrow = 1, ncol = p)
    }
    else {
      ar.coefs <- NULL
      p <- 0
    }
    if (!is.null(paramz$ma.coefs)) {
      q <- length(paramz$ma.coefs)
      ma.coefs <- matrix(paramz$ma.coefs, nrow = 1, ncol = q)
    }
    else {
      ma.coefs <- NULL
      q <- 0
    }
    x.nought <- BoxCox(opt.env$x.nought.untransformed, lambda = lambda)
    lambda <- attr(x.nought, "lambda")
    w <- .Call("makeTBATSWMatrix", smallPhi_s = small.phi, 
               kVector_s = k.vector, arCoefs_s = ar.coefs, maCoefs_s = ma.coefs, 
               tau_s = tau, PACKAGE = "fable.tbats")
    if (!is.null(gamma.bold)) {
      .Call("updateTBATSGammaBold", gammaBold_s = gamma.bold, 
            kVector_s = k.vector, gammaOne_s = gamma.one.v, 
            gammaTwo_s = gamma.two.v, PACKAGE = "fable.tbats")
    }
    .Call("updateTBATSGMatrix", g_s = g, gammaBold_s = gamma.bold, 
          alpha_s = alpha, beta_s = beta.v, PACKAGE = "fable.tbats")
    .Call("updateFMatrix", F, small.phi, alpha, beta.v, gamma.bold, 
          ar.coefs, ma.coefs, tau, PACKAGE = "fable.tbats")
    y.transformed <- BoxCox(y, lambda = lambda)
    lambda <- attr(y.transformed, "lambda")
    fitted.values.and.errors <- calcModel(y.transformed, 
                                          x.nought, F, g, w)
    e <- fitted.values.and.errors$e
    variance <- sum((e * e))/length(y)
    fitted.values <- InvBoxCox(fitted.values.and.errors$y.hat, 
                               lambda = lambda, biasadj, variance)
    attr(lambda, "biasadj") <- biasadj
    ee <- y - fitted.values
  }
  else {
    if (length(param.vector$vect) > 1) {
      optim.like <- optim(par = param.vector$vect, fn = calcLikelihoodNOTransformedTBATS, 
                          method = "Nelder-Mead", opt.env = opt.env, x.nought = x.nought, 
                          use.beta = use.beta, use.small.phi = use.damping, 
                          seasonal.periods = seasonal.periods, param.control = param.vector$control, 
                          p = p, q = q, tau = tau, control = list(maxit = (100 * 
                                                                             length(param.vector$vect)^2), parscale = par.scale))
    }
    else {
      optim.like <- optim(par = param.vector$vect, fn = calcLikelihoodNOTransformedTBATS, 
                          method = "BFGS", opt.env = opt.env, x.nought = x.nought, 
                          use.beta = use.beta, use.small.phi = use.damping, 
                          seasonal.periods = seasonal.periods, param.control = param.vector$control, 
                          p = p, q = q, tau = tau, control = list(parscale = par.scale))
    }
    paramz <- unParameteriseTBATS(optim.like$par, param.vector$control)
    lambda <- paramz$lambda
    alpha <- paramz$alpha
    beta.v <- paramz$beta
    small.phi <- paramz$small.phi
    gamma.one.v <- paramz$gamma.one.v
    gamma.two.v <- paramz$gamma.two.v
    if (!is.null(paramz$ar.coefs)) {
      p <- length(paramz$ar.coefs)
      ar.coefs <- matrix(paramz$ar.coefs, nrow = 1, ncol = p)
    }
    else {
      ar.coefs <- NULL
      p <- 0
    }
    if (!is.null(paramz$ma.coefs)) {
      q <- length(paramz$ma.coefs)
      ma.coefs <- matrix(paramz$ma.coefs, nrow = 1, ncol = q)
    }
    else {
      ma.coefs <- NULL
      q <- 0
    }
    w <- .Call("makeTBATSWMatrix", smallPhi_s = small.phi, 
               kVector_s = k.vector, arCoefs_s = ar.coefs, maCoefs_s = ma.coefs, 
               tau_s = tau, PACKAGE = "fable.tbats")
    if (!is.null(gamma.bold)) {
      .Call("updateTBATSGammaBold", gammaBold_s = gamma.bold, 
            kVector_s = k.vector, gammaOne_s = gamma.one.v, 
            gammaTwo_s = gamma.two.v, PACKAGE = "fable.tbats")
    }
    .Call("updateTBATSGMatrix", g_s = g, gammaBold_s = gamma.bold, 
          alpha_s = alpha, beta_s = beta.v, PACKAGE = "fable.tbats")
    .Call("updateFMatrix", F, small.phi, alpha, beta.v, gamma.bold, 
          ar.coefs, ma.coefs, tau, PACKAGE = "fable.tbats")
    fitted.values.and.errors <- calcModel(y, x.nought, F, 
                                          g, w)
    e <- fitted.values.and.errors$e
    fitted.values <- fitted.values.and.errors$y.hat
    variance <- sum((e * e))/length(y)
  }
  likelihood <- optim.like$value
  aic <- likelihood + 2 * (length(param.vector$vect) + nrow(x.nought))
  fits <- ts(c(fitted.values))
  e <- ts(c(e))
  tsp(fits) <- tsp(e) <- tsp(y)
  model.for.output <- list(lambda = lambda, alpha = alpha, 
                           beta = beta.v, damping.parameter = small.phi, gamma.one.values = gamma.one.v, 
                           gamma.two.values = gamma.two.v, ar.coefficients = ar.coefs, 
                           ma.coefficients = ma.coefs, likelihood = likelihood, 
                           optim.return.code = optim.like$convergence, variance = variance, 
                           AIC = aic, parameters = list(vect = optim.like$par, control = param.vector$control), 
                           seed.states = x.nought, fitted.values = fits, errors = e, 
                           x = fitted.values.and.errors$x, seasonal.periods = seasonal.periods, 
                           k.vector = k.vector, y = y, p = p, q = q)
  class(model.for.output) <- c("tbats", "bats")
  return(model.for.output)
}

makeParscale <- function (control) 
{
  if (control$use.box.cox) {
    parscale <- c(0.001, 0.01)
  }
  else {
    parscale <- 0.01
  }
  if (control$use.beta) {
    if (control$use.damping) {
      parscale <- c(parscale, 0.01, 0.01)
    }
    else {
      parscale <- c(parscale, 0.01)
    }
  }
  if (control$length.gamma > 0) {
    parscale <- c(parscale, rep(1e-05, control$length.gamma))
  }
  if ((control$p != 0) | (control$q != 0)) {
    parscale <- c(parscale, rep(0.1, (control$p + control$q)))
  }
  return(parscale)
}

makeParscaleBATS <- function (control) 
{
  if (control$use.box.cox) {
    parscale <- c(0.001, 0.1)
  }
  else {
    parscale <- 0.1
  }
  if (control$use.beta) {
    if (control$use.damping) {
      parscale <- c(parscale, 0.01, 0.01)
    }
    else {
      parscale <- c(parscale, 0.01)
    }
  }
  if (control$length.gamma > 0) {
    parscale <- c(parscale, rep(0.01, control$length.gamma))
  }
  if ((control$p != 0) | (control$q != 0)) {
    parscale <- c(parscale, rep(0.1, (control$p + control$q)))
  }
  return(parscale)
}

filterTBATSSpecifics <- function (y, box.cox, trend, damping, seasonal.periods, k.vector, 
                                  use.arma.errors, aux.model = NULL, init.box.cox = NULL, bc.lower = 0, 
                                  bc.upper = 1, biasadj = FALSE, ...) 
{
  if (is.null(aux.model)) {
    first.model <- try(fitSpecificTBATS(y, use.box.cox = box.cox, 
                                        use.beta = trend, use.damping = damping, seasonal.periods = seasonal.periods, 
                                        k.vector = k.vector, init.box.cox = init.box.cox, 
                                        bc.lower = bc.lower, bc.upper = bc.upper, biasadj = biasadj), 
                       silent = TRUE)
  }
  else {
    first.model <- aux.model
  }
  if (is.element("try-error", class(first.model))) {
    first.model <- list(AIC = Inf)
  }
  if (use.arma.errors) {
    suppressWarnings(arma <- try(auto_arma(as.numeric(first.model$errors), 
                                           d = 0, ...), silent = TRUE))
    if (!is.element("try-error", class(arma))) {
      p <- arma$arma[1]
      q <- arma$arma[2]
      if ((p != 0) || (q != 0)) {
        if (p != 0) {
          ar.coefs <- numeric(p)
        }
        else {
          ar.coefs <- NULL
        }
        if (q != 0) {
          ma.coefs <- numeric(q)
        }
        else {
          ma.coefs <- NULL
        }
        starting.params <- first.model$parameters
        second.model <- try(fitSpecificTBATS(y, use.box.cox = box.cox, 
                                             use.beta = trend, use.damping = damping, seasonal.periods = seasonal.periods, 
                                             k.vector = k.vector, ar.coefs = ar.coefs, ma.coefs = ma.coefs, 
                                             init.box.cox = init.box.cox, bc.lower = bc.lower, 
                                             bc.upper = bc.upper, biasadj = biasadj), silent = TRUE)
        if (is.element("try-error", class(second.model))) {
          second.model <- list(AIC = Inf)
        }
        if (second.model$AIC < first.model$AIC) {
          return(second.model)
        }
        else {
          return(first.model)
        }
      }
      else {
        return(first.model)
      }
    }
    else {
      return(first.model)
    }
  }
  else {
    return(first.model)
  }
}

is.constant <- function (x) 
{
  x <- as.numeric(x)
  y <- rep(x[1], length(x))
  return(isTRUE(all.equal(x, y)))
}

BoxCox <- function (x, lambda) 
{
  if (lambda == "auto") {
    lambda <- BoxCox.lambda(x, lower = -0.9)
  }
  if (lambda < 0) {
    x[x < 0] <- NA
  }
  if (lambda == 0) {
    out <- log(x)
  }
  else {
    out <- (sign(x) * abs(x)^lambda - 1)/lambda
  }
  if (!is.null(colnames(x))) {
    colnames(out) <- colnames(x)
  }
  attr(out, "lambda") <- lambda
  return(out)
}

InvBoxCox <- function (x, lambda, biasadj = FALSE, fvar = NULL) 
{
  if (lambda < 0) {
    x[x > -1/lambda] <- NA
  }
  if (lambda == 0) {
    out <- exp(x)
  }
  else {
    xx <- x * lambda + 1
    out <- sign(xx) * abs(xx)^(1/lambda)
  }
  if (!is.null(colnames(x))) {
    colnames(out) <- colnames(x)
  }
  if (is.null(biasadj)) {
    biasadj <- attr(lambda, "biasadj")
  }
  if (!is.logical(biasadj)) {
    warning("biasadj information not found, defaulting to FALSE.")
    biasadj <- FALSE
  }
  if (biasadj) {
    if (is.null(fvar)) {
      stop("fvar must be provided when biasadj=TRUE")
    }
    if (is.list(fvar)) {
      level <- max(fvar$level)
      if (NCOL(fvar$upper) > 1 && NCOL(fvar$lower)) {
        i <- match(level, fvar$level)
        fvar$upper <- fvar$upper[, i]
        fvar$lower <- fvar$lower[, i]
      }
      if (level > 1) {
        level <- level/100
      }
      level <- mean(c(level, 1))
      fvar <- as.numeric((fvar$upper - fvar$lower)/stats::qnorm(level)/2)^2
    }
    if (NCOL(fvar) > 1) {
      fvar <- diag(fvar)
    }
    out <- out * (1 + 0.5 * as.numeric(fvar) * (1 - lambda)/(out)^(2 * 
                                                                     lambda))
  }
  return(out)
}

parameterise <- function (alpha, beta.v = NULL, small.phi = 1, gamma.v = NULL, 
                          lambda = NULL, ar.coefs = NULL, ma.coefs = NULL) 
{
  if (!is.null(lambda)) {
    param.vector <- cbind(lambda, alpha)
    use.box.cox <- TRUE
  }
  else {
    param.vector <- alpha
    use.box.cox <- FALSE
  }
  if (!is.null(beta.v)) {
    use.beta <- TRUE
    if (is.null(small.phi)) {
      use.damping <- FALSE
    }
    else if (small.phi != 1) {
      param.vector <- cbind(param.vector, small.phi)
      use.damping <- TRUE
    }
    else {
      use.damping <- FALSE
    }
    param.vector <- cbind(param.vector, beta.v)
  }
  else {
    use.beta <- FALSE
    use.damping <- FALSE
  }
  if (!is.null(gamma.v)) {
    gamma.v <- matrix(gamma.v, nrow = 1, ncol = length(gamma.v))
    param.vector <- cbind(param.vector, gamma.v)
    length.gamma <- length(gamma.v)
  }
  else {
    length.gamma <- 0
  }
  if (!is.null(ar.coefs)) {
    ar.coefs <- matrix(ar.coefs, nrow = 1, ncol = length(ar.coefs))
    param.vector <- cbind(param.vector, ar.coefs)
    p <- length(ar.coefs)
  }
  else {
    p <- 0
  }
  if (!is.null(ma.coefs)) {
    ma.coefs <- matrix(ma.coefs, nrow = 1, ncol = length(ma.coefs))
    param.vector <- cbind(param.vector, ma.coefs)
    q <- length(ma.coefs)
  }
  else {
    q <- 0
  }
  control <- list(use.beta = use.beta, use.box.cox = use.box.cox, 
                  use.damping = use.damping, length.gamma = length.gamma, 
                  p = p, q = q)
  return(list(vect = as.numeric(param.vector), control = control))
}

makeParscaleBATS <- function (control) 
{
  if (control$use.box.cox) {
    parscale <- c(0.001, 0.1)
  }
  else {
    parscale <- 0.1
  }
  if (control$use.beta) {
    if (control$use.damping) {
      parscale <- c(parscale, 0.01, 0.01)
    }
    else {
      parscale <- c(parscale, 0.01)
    }
  }
  if (control$length.gamma > 0) {
    parscale <- c(parscale, rep(0.01, control$length.gamma))
  }
  if ((control$p != 0) | (control$q != 0)) {
    parscale <- c(parscale, rep(0.1, (control$p + control$q)))
  }
  return(parscale)
}


guerrero <- function (x, lower = -1, upper = 2, nonseasonal.length = 2) 
{
  if (any(x <= 0, na.rm = TRUE)) 
    warning("Guerrero's method for selecting a Box-Cox parameter (lambda) is given for strictly positive data.")
  return(optimize(guer.cv, c(lower, upper), x = x, nonseasonal.length = nonseasonal.length)$minimum)
}

guer.cv <- function (lam, x, nonseasonal.length = 2) 
{
  period <- round(max(nonseasonal.length, frequency(x)))
  nobsf <- length(x)
  nyr <- floor(nobsf/period)
  nobst <- floor(nyr * period)
  x.mat <- matrix(x[(nobsf - nobst + 1):nobsf], period, nyr)
  x.mean <- apply(x.mat, 2, mean, na.rm = TRUE)
  x.sd <- apply(x.mat, 2, sd, na.rm = TRUE)
  x.rat <- x.sd/x.mean^(1 - lam)
  return(sd(x.rat, na.rm = TRUE)/mean(x.rat, na.rm = TRUE))
}

BoxCox.lambda <- function (x, lower = -1, upper = 2) 
{
  if (any(x <= 0, na.rm = TRUE)) {
    lower <- max(lower, 0)
  }
  if (length(x) <= 2 * frequency(x)) {
    return(1)
  }
  return(guerrero(x, lower, upper))
}

tbats <- function (y, use.box.cox = NULL, use.trend = NULL, use.damped.trend = NULL, 
                   seasonal.periods = NULL, use.arma.errors = TRUE, use.parallel = length(y) > 
                     1000, num.cores = 2, bc.lower = 0, bc.upper = 1, biasadj = FALSE, 
                   model = NULL, ...) 
{
  if (!is.numeric(y) || NCOL(y) > 1) {
    stop("y should be a univariate time series")
  }
  seriesname <- deparse(substitute(y))
  origy <- y
  attr_y <- attributes(origy)
  if (is.null(seasonal.periods)) {
    if (any(class(y) == "msts")) {
      seasonal.periods <- sort(attr(y, "msts"))
    }
    else if (class(y) == "ts") {
      seasonal.periods <- frequency(y)
    }
    else {
      y <- as.ts(y)
      seasonal.periods <- 1
    }
  }
  else {
    if (!any(class(y) == "ts")) {
      y <- msts(y, seasonal.periods)
    }
  }
  seasonal.periods <- unique(pmax(seasonal.periods, 1))
  if (all(seasonal.periods == 1)) {
    seasonal.periods <- NULL
  }
  ny <- length(y)
  y <- na.contiguous(y)
  if (ny != length(y)) {
    warning("Missing values encountered. Using longest contiguous portion of time series")
    if (!is.null(attr_y$tsp)) {
      attr_y$tsp[1:2] <- range(time(y))
    }
  }
  if (!is.null(model)) {
    if (is.element("tbats", class(model))) {
      refitModel <- try(fitPreviousTBATSModel(y, model = model), 
                        silent = TRUE)
    }
    else if (is.element("bats", class(model))) {
      refitModel <- bats(origy, model = model)
    }
    return(refitModel)
  }
  if (is.constant(y)) {
    fit <- list(y = y, x = matrix(y, nrow = 1, ncol = ny), 
                errors = y * 0, fitted.values = y, seed.states = matrix(y[1]), 
                AIC = -Inf, likelihood = -Inf, variance = 0, alpha = 0.9999, 
                method = "TBATS", call = match.call())
    return(structure(fit, class = "bats"))
  }
  if (any((y <= 0))) {
    use.box.cox <- FALSE
  }
  non.seasonal.model <- bats(as.numeric(y), use.box.cox = use.box.cox, 
                             use.trend = use.trend, use.damped.trend = use.damped.trend, 
                             use.arma.errors = use.arma.errors, use.parallel = use.parallel, 
                             num.cores = num.cores, bc.lower = bc.lower, bc.upper = bc.upper, 
                             biasadj = biasadj, ...)
  if (is.null(seasonal.periods)) {
    non.seasonal.model$call <- match.call()
    attributes(non.seasonal.model$fitted.values) <- attributes(non.seasonal.model$errors) <- attributes(origy)
    non.seasonal.model$y <- origy
    return(non.seasonal.model)
  }
  else {
    seasonal.mask <- (seasonal.periods == 1)
    seasonal.periods <- seasonal.periods[!seasonal.mask]
  }
  if (is.null(use.box.cox)) {
    use.box.cox <- c(FALSE, TRUE)
  }
  if (any(use.box.cox)) {
    init.box.cox <- BoxCox.lambda(y, lower = bc.lower, upper = bc.upper)
  }
  else {
    init.box.cox <- NULL
  }
  if (is.null(use.trend)) {
    use.trend <- c(FALSE, TRUE)
  }
  else if (use.trend == FALSE) {
    use.damped.trend <- FALSE
  }
  if (is.null(use.damped.trend)) {
    use.damped.trend <- c(FALSE, TRUE)
  }
  model.params <- logical(length = 3)
  model.params[1] <- any(use.box.cox)
  model.params[2] <- any(use.trend)
  model.params[3] <- any(use.damped.trend)
  y <- as.numeric(y)
  n <- length(y)
  k.vector <- rep(1, length(seasonal.periods))
  if (use.parallel) {
    if (is.null(num.cores)) {
      num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
    }
    clus <- makeCluster(num.cores)
  }
  best.model <- try(fitSpecificTBATS(y, use.box.cox = model.params[1], 
                                     use.beta = model.params[2], use.damping = model.params[3], 
                                     seasonal.periods = seasonal.periods, k.vector = k.vector, 
                                     init.box.cox = init.box.cox, bc.lower = bc.lower, bc.upper = bc.upper, 
                                     biasadj = biasadj), silent = TRUE)
  if (is.element("try-error", class(best.model))) {
    best.model <- list(AIC = Inf)
  }
  for (i in 1:length(seasonal.periods)) {
    if (seasonal.periods[i] == 2) {
      next
    }
    max.k <- floor(((seasonal.periods[i] - 1)/2))
    if (i != 1) {
      current.k <- 2
      while (current.k <= max.k) {
        if (seasonal.periods[i]%%current.k != 0) {
          current.k <- current.k + 1
          next
        }
        latter <- seasonal.periods[i]/current.k
        if (any(((seasonal.periods[1:(i - 1)]%%latter) == 
                 0))) {
          max.k <- current.k - 1
          break
        }
        else {
          current.k <- current.k + 1
        }
      }
    }
    if (max.k == 1) {
      next
    }
    if (max.k <= 6) {
      k.vector[i] <- max.k
      best.model$AIC <- Inf
      repeat {
        new.model <- try(fitSpecificTBATS(y, use.box.cox = model.params[1], 
                                          use.beta = model.params[2], use.damping = model.params[3], 
                                          seasonal.periods = seasonal.periods, k.vector = k.vector, 
                                          init.box.cox = init.box.cox, bc.lower = bc.lower, 
                                          bc.upper = bc.upper, biasadj = biasadj), silent = TRUE)
        if (is.element("try-error", class(new.model))) {
          new.model <- list(AIC = Inf)
        }
        if (new.model$AIC > best.model$AIC) {
          k.vector[i] <- k.vector[i] + 1
          break
        }
        else {
          if (k.vector[i] == 1) {
            break
          }
          k.vector[i] <- k.vector[i] - 1
          best.model <- new.model
        }
      }
      next
    }
    else {
      step.up.k <- k.vector
      step.down.k <- k.vector
      step.up.k[i] <- 7
      step.down.k[i] <- 5
      k.vector[i] <- 6
      if (use.parallel) {
        k.control.array <- rbind(step.up.k, step.down.k, 
                                 k.vector)
        models.list <- clusterApplyLB(clus, c(1:3), parFitSpecificTBATS, 
                                      y = y, box.cox = model.params[1], trend = model.params[2], 
                                      damping = model.params[3], seasonal.periods = seasonal.periods, 
                                      k.control.matrix = k.control.array, init.box.cox = init.box.cox, 
                                      bc.lower = bc.lower, bc.upper = bc.upper, biasadj = biasadj)
        up.model <- models.list[[1]]
        level.model <- models.list[[3]]
        down.model <- models.list[[2]]
      }
      else {
        up.model <- try(fitSpecificTBATS(y, use.box.cox = model.params[1], 
                                         use.beta = model.params[2], use.damping = model.params[3], 
                                         seasonal.periods = seasonal.periods, k.vector = step.up.k, 
                                         init.box.cox = init.box.cox, bc.lower = bc.lower, 
                                         bc.upper = bc.upper, biasadj = biasadj), silent = TRUE)
        if (is.element("try-error", class(up.model))) {
          up.model <- list(AIC = Inf)
        }
        level.model <- try(fitSpecificTBATS(y, use.box.cox = model.params[1], 
                                            use.beta = model.params[2], use.damping = model.params[3], 
                                            seasonal.periods = seasonal.periods, k.vector = k.vector, 
                                            init.box.cox = init.box.cox, bc.lower = bc.lower, 
                                            bc.upper = bc.upper, biasadj = biasadj), silent = TRUE)
        if (is.element("try-error", class(level.model))) {
          level.model <- list(AIC = Inf)
        }
        down.model <- try(fitSpecificTBATS(y, use.box.cox = model.params[1], 
                                           use.beta = model.params[2], use.damping = model.params[3], 
                                           seasonal.periods = seasonal.periods, k.vector = step.down.k, 
                                           init.box.cox = init.box.cox, bc.lower = bc.lower, 
                                           bc.upper = bc.upper, biasadj = biasadj), silent = TRUE)
        if (is.element("try-error", class(down.model))) {
          down.model <- list(AIC = Inf)
        }
      }
      aic.vector <- c(up.model$AIC, level.model$AIC, down.model$AIC)
      if (min(aic.vector) == down.model$AIC) {
        best.model <- down.model
        k.vector[i] <- 5
        repeat {
          k.vector[i] <- k.vector[i] - 1
          down.model <- try(fitSpecificTBATS(y = y, use.box.cox = model.params[1], 
                                             use.beta = model.params[2], use.damping = model.params[3], 
                                             seasonal.periods = seasonal.periods, k.vector = k.vector, 
                                             init.box.cox = init.box.cox, bc.lower = bc.lower, 
                                             bc.upper = bc.upper, biasadj = biasadj), 
                            silent = TRUE)
          if (is.element("try-error", class(down.model))) {
            down.model <- list(AIC = Inf)
          }
          if (down.model$AIC > best.model$AIC) {
            k.vector[i] <- k.vector[i] + 1
            break
          }
          else {
            best.model <- down.model
          }
          if (k.vector[i] == 1) {
            break
          }
        }
      }
      else if (min(aic.vector) == level.model$AIC) {
        best.model <- level.model
        next
      }
      else {
        best.model <- up.model
        k.vector[i] <- 7
        repeat {
          k.vector[i] <- k.vector[i] + 1
          up.model <- try(fitSpecificTBATS(y, model.params[1], 
                                           model.params[2], model.params[3], seasonal.periods, 
                                           k.vector, init.box.cox = init.box.cox, bc.lower = bc.lower, 
                                           bc.upper = bc.upper, biasadj = biasadj), 
                          silent = TRUE)
          if (is.element("try-error", class(up.model))) {
            up.model <- list(AIC = Inf)
          }
          if (up.model$AIC > best.model$AIC) {
            k.vector[i] <- k.vector[i] - 1
            break
          }
          else {
            best.model <- up.model
          }
          if (k.vector[i] == max.k) {
            break
          }
        }
      }
    }
  }
  aux.model <- best.model
  if (non.seasonal.model$AIC < best.model$AIC) {
    best.model <- non.seasonal.model
  }
  if ((length(use.box.cox) == 1) && use.trend[1] && (length(use.trend) == 
                                                     1) && (length(use.damped.trend) == 1) && (use.parallel)) {
    use.parallel <- FALSE
    stopCluster(clus)
  }
  else if ((length(use.box.cox) == 1) && !use.trend[1] && (length(use.trend) == 
                                                           1) && (use.parallel)) {
    use.parallel <- FALSE
    stopCluster(clus)
  }
  if (use.parallel) {
    control.array <- NULL
    for (box.cox in use.box.cox) {
      for (trend in use.trend) {
        for (damping in use.damped.trend) {
          if (!trend && damping) {
            next
          }
          control.line <- c(box.cox, trend, damping)
          if (!is.null(control.array)) {
            control.array <- rbind(control.array, control.line)
          }
          else {
            control.array <- control.line
          }
        }
      }
    }
    models.list <- clusterApplyLB(clus, c(1:nrow(control.array)), 
                                  parFilterTBATSSpecifics, y = y, control.array = control.array, 
                                  model.params = model.params, seasonal.periods = seasonal.periods, 
                                  k.vector = k.vector, use.arma.errors = use.arma.errors, 
                                  aux.model = aux.model, init.box.cox = init.box.cox, 
                                  bc.lower = bc.lower, bc.upper = bc.upper, biasadj = biasadj, 
                                  ...)
    stopCluster(clus)
    aics <- numeric(nrow(control.array))
    for (i in 1:nrow(control.array)) {
      aics[i] <- models.list[[i]]$AIC
    }
    best.number <- which.min(aics)
    best.seasonal.model <- models.list[[best.number]]
    if (best.seasonal.model$AIC < best.model$AIC) {
      best.model <- best.seasonal.model
    }
  }
  else {
    for (box.cox in use.box.cox) {
      for (trend in use.trend) {
        for (damping in use.damped.trend) {
          if (all((model.params == c(box.cox, trend, 
                                     damping)))) {
            new.model <- filterTBATSSpecifics(y, box.cox, 
                                              trend, damping, seasonal.periods, k.vector, 
                                              use.arma.errors, aux.model = aux.model, 
                                              init.box.cox = init.box.cox, bc.lower = bc.lower, 
                                              bc.upper = bc.upper, biasadj = biasadj, 
                                              ...)
          }
          else if (trend || !damping) {
            new.model <- filterTBATSSpecifics(y, box.cox, 
                                              trend, damping, seasonal.periods, k.vector, 
                                              use.arma.errors, init.box.cox = init.box.cox, 
                                              bc.lower = bc.lower, bc.upper = bc.upper, 
                                              biasadj = biasadj, ...)
          }
          if (new.model$AIC < best.model$AIC) {
            best.model <- new.model
          }
        }
      }
    }
  }
  best.model$call <- match.call()
  attributes(best.model$fitted.values) <- attributes(best.model$errors) <- attr_y
  best.model$y <- origy
  best.model$series <- seriesname
  best.model$method <- "TBATS"
  return(best.model)
}

bats <- function (y, use.box.cox = NULL, use.trend = NULL, use.damped.trend = NULL, 
                  seasonal.periods = NULL, use.arma.errors = TRUE, use.parallel = length(y) > 
                    1000, num.cores = 2, bc.lower = 0, bc.upper = 1, biasadj = FALSE, 
                  model = NULL, ...) 
{
  if (!is.numeric(y) || NCOL(y) > 1) {
    stop("y should be a univariate time series")
  }
  seriesname <- deparse(substitute(y))
  origy <- y
  attr_y <- attributes(origy)
  if (is.null(seasonal.periods)) {
    if (any(class(y) == "msts")) {
      seasonal.periods <- attr(y, "msts")
    }
    else if (class(y) == "ts") {
      seasonal.periods <- frequency(y)
    }
    else {
      y <- as.ts(y)
      seasonal.periods <- 1
    }
    seasonal.periods <- seasonal.periods[seasonal.periods < 
                                           length(y)]
    if (length(seasonal.periods) == 0L) 
      seasonal.periods <- 1
  }
  else {
    if (!any(class(y) == "ts")) {
      y <- msts(y, seasonal.periods)
    }
  }
  seasonal.periods <- unique(pmax(seasonal.periods, 1))
  if (all(seasonal.periods == 1)) {
    seasonal.periods <- NULL
  }
  ny <- length(y)
  y <- na.contiguous(y)
  if (ny != length(y)) {
    warning("Missing values encountered. Using longest contiguous portion of time series")
    if (!is.null(attr_y$tsp)) {
      attr_y$tsp[1:2] <- range(time(y))
    }
  }
  if (!is.null(model)) {
    refitModel <- try(fitPreviousBATSModel(y, model = model), 
                      silent = TRUE)
    return(refitModel)
  }
  if (is.constant(y)) {
    fit <- list(y = y, x = matrix(y, nrow = 1, ncol = ny), 
                errors = y * 0, fitted.values = y, seed.states = matrix(y[1]), 
                AIC = -Inf, likelihood = -Inf, variance = 0, alpha = 0.9999, 
                method = "BATS", call = match.call())
    return(structure(fit, class = "bats"))
  }
  if (any((y <= 0))) {
    use.box.cox <- FALSE
  }
  if ((!is.null(use.box.cox)) && (!is.null(use.trend)) && (use.parallel)) {
    if (use.trend && (!is.null(use.damped.trend))) {
      use.parallel <- FALSE
    }
    else if (use.trend == FALSE) {
      use.parallel <- FALSE
    }
  }
  if (!is.null(seasonal.periods)) {
    seasonal.mask <- (seasonal.periods == 1)
    seasonal.periods <- seasonal.periods[!seasonal.mask]
  }
  if (is.null(seasonal.periods) && !is.null(use.box.cox) && 
      !is.null(use.trend)) {
    use.parallel <- FALSE
  }
  if (is.null(use.box.cox)) {
    use.box.cox <- c(FALSE, TRUE)
  }
  if (any(use.box.cox)) {
    init.box.cox <- BoxCox.lambda(y, lower = bc.lower, upper = bc.upper)
  }
  else {
    init.box.cox <- NULL
  }
  if (is.null(use.trend)) {
    use.trend <- c(FALSE, TRUE)
  }
  else if (use.trend == FALSE) {
    use.damped.trend <- FALSE
  }
  if (is.null(use.damped.trend)) {
    use.damped.trend <- c(FALSE, TRUE)
  }
  y <- as.numeric(y)
  if (use.parallel) {
    control.array <- NULL
    for (box.cox in use.box.cox) {
      for (trend in use.trend) {
        for (damping in use.damped.trend) {
          if (!trend && damping) {
            next
          }
          control.line <- c(box.cox, trend, damping)
          if (!is.null(control.array)) {
            control.array <- rbind(control.array, control.line)
          }
          else {
            control.array <- control.line
          }
        }
      }
    }
    if (is.null(num.cores)) {
      num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
    }
    clus <- makeCluster(num.cores)
    models.list <- clusterApplyLB(clus, c(1:nrow(control.array)), 
                                  parFilterSpecifics, y = y, control.array = control.array, 
                                  seasonal.periods = seasonal.periods, use.arma.errors = use.arma.errors, 
                                  init.box.cox = init.box.cox, bc.lower = bc.lower, 
                                  bc.upper = bc.upper, biasadj = biasadj, ...)
    stopCluster(clus)
    aics <- numeric(nrow(control.array))
    for (i in 1:nrow(control.array)) {
      aics[i] <- models.list[[i]]$AIC
    }
    best.number <- which.min(aics)
    best.model <- models.list[[best.number]]
  }
  else {
    best.aic <- Inf
    best.model <- NULL
    for (box.cox in use.box.cox) {
      for (trend in use.trend) {
        for (damping in use.damped.trend) {
          current.model <- try(filterSpecifics(y, box.cox = box.cox, 
                                               trend = trend, damping = damping, seasonal.periods = seasonal.periods, 
                                               use.arma.errors = use.arma.errors, init.box.cox = init.box.cox, 
                                               bc.lower = bc.lower, bc.upper = bc.upper, 
                                               biasadj = biasadj, ...), silent = FALSE)#TRUE)
          if (!("try-error" %in% class(current.model))) {
            if (current.model$AIC < best.aic) {
              best.aic <- current.model$AIC
              best.model <- current.model
            }
          }
        }
      }
    }
  }
  if (is.null(best.model)) 
    stop("Unable to fit a model")
  best.model$call <- match.call()
  if (best.model$optim.return.code != 0) {
    warning("optim() did not converge.")
  }
  attributes(best.model$fitted.values) <- attributes(best.model$errors) <- attr_y
  best.model$y <- origy
  best.model$series <- seriesname
  best.model$method <- "BATS"
  return(best.model)
}

auto_arma <- function( x, max_p = 5, max_q = 5, ... ) {
  
  mdls <- vector("list", (max_p+1)*(max_q+1))
  id <- 1
  for( p in 0:max_p ) {
    for( q in 0:max_q) {
      mdls[[id]] <- try(suppressWarnings(stats::arima(x, order = c(p, 0, q), method = "ML")),
                        silent = TRUE)
      id <- id+1
    }
  }
  best_mdl <- which.min(sapply(mdls, function(i){i$aic}))
  return(mdls[[best_mdl]])
}

makeXMatrix <- function (l, b = NULL, s.vector = NULL, d.vector = NULL, epsilon.vector = NULL) 
{
  x.transpose <- matrix(l, nrow = 1, ncol = 1)
  if (!is.null(b)) {
    x.transpose <- cbind(x.transpose, matrix(b, nrow = 1, 
                                             ncol = 1))
  }
  if (!is.null(s.vector)) {
    x.transpose <- cbind(x.transpose, matrix(s.vector, nrow = 1, 
                                             ncol = length(s.vector)))
  }
  if (!is.null(d.vector)) {
    x.transpose <- cbind(x.transpose, matrix(d.vector, nrow = 1, 
                                             ncol = length(d.vector)))
  }
  if (!is.null(epsilon.vector)) {
    x.transpose <- cbind(x.transpose, matrix(epsilon.vector, 
                                             nrow = 1, ncol = length(epsilon.vector)))
  }
  x <- t(x.transpose)
  return(list(x = x, x.transpose = x.transpose))
}

makeFMatrix <- function (alpha, beta = NULL, small.phi = NULL, seasonal.periods = NULL, 
                         gamma.bold.matrix = NULL, ar.coefs = NULL, ma.coefs = NULL) 
{
  F <- matrix(1, nrow = 1, ncol = 1)
  if (!is.null(beta)) {
    F <- cbind(F, matrix(small.phi, nrow = 1, ncol = 1))
  }
  if (!is.null(seasonal.periods)) {
    tau <- sum(seasonal.periods)
    zero.tau <- matrix(0, nrow = 1, ncol = tau)
    F <- cbind(F, zero.tau)
  }
  if (!is.null(ar.coefs)) {
    p <- length(ar.coefs)
    ar.coefs <- matrix(ar.coefs, nrow = 1, ncol = p)
    alpha.phi <- alpha * ar.coefs
    F <- cbind(F, alpha.phi)
  }
  if (!is.null(ma.coefs)) {
    q <- length(ma.coefs)
    ma.coefs <- matrix(ma.coefs, nrow = 1, ncol = q)
    alpha.theta <- alpha * ma.coefs
    F <- cbind(F, alpha.theta)
  }
  if (!is.null(beta)) {
    beta.row <- matrix(c(0, small.phi), nrow = 1, ncol = 2)
    if (!is.null(seasonal.periods)) {
      beta.row <- cbind(beta.row, zero.tau)
    }
    if (!is.null(ar.coefs)) {
      beta.phi <- beta * ar.coefs
      beta.row <- cbind(beta.row, beta.phi)
    }
    if (!is.null(ma.coefs)) {
      beta.theta <- beta * ma.coefs
      beta.row <- cbind(beta.row, beta.theta)
    }
    F <- rbind(F, beta.row)
  }
  if (!is.null(seasonal.periods)) {
    seasonal.row <- t(zero.tau)
    if (!is.null(beta)) {
      seasonal.row <- cbind(seasonal.row, seasonal.row)
    }
    for (i in seasonal.periods) {
      if (i == seasonal.periods[1]) {
        a.row.one <- matrix(0, nrow = 1, ncol = i)
        a.row.one[i] <- 1
        a.row.two <- cbind(diag((i - 1)), matrix(0, nrow = (i - 
                                                              1), ncol = 1))
        A <- rbind(a.row.one, a.row.two)
      }
      else {
        old.A.rows <- dim(A)[1]
        old.A.columns <- dim(A)[2]
        a.row.one <- matrix(0, nrow = 1, ncol = i)
        a.row.one[i] <- 1
        a.row.two <- cbind(diag((i - 1)), matrix(0, nrow = (i - 
                                                              1), ncol = 1))
        Ai <- rbind(a.row.one, a.row.two)
        A <- rbind(A, matrix(0, nrow = dim(Ai)[1], ncol = old.A.columns))
        A <- cbind(A, matrix(0, nrow = dim(A)[1], ncol = dim(Ai)[2]))
        A[((old.A.rows + 1):(old.A.rows + dim(Ai)[1])), 
          ((old.A.columns + 1):(old.A.columns + dim(Ai)[2]))] <- Ai
      }
    }
    seasonal.row <- cbind(seasonal.row, A)
    if (!is.null(ar.coefs)) {
      B <- t(gamma.bold.matrix) %*% ar.coefs
      seasonal.row <- cbind(seasonal.row, B)
    }
    if (!is.null(ma.coefs)) {
      C <- t(gamma.bold.matrix) %*% ma.coefs
      seasonal.row <- cbind(seasonal.row, C)
    }
    F <- rbind(F, seasonal.row)
  }
  if (!is.null(ar.coefs)) {
    ar.rows <- matrix(0, nrow = p, ncol = 1)
    if (!is.null(beta)) {
      ar.rows <- cbind(ar.rows, ar.rows)
    }
    if (!is.null(seasonal.periods)) {
      ar.seasonal.zeros <- matrix(0, nrow = p, ncol = tau)
      ar.rows <- cbind(ar.rows, ar.seasonal.zeros)
    }
    ident <- diag((p - 1))
    ident <- cbind(ident, matrix(0, nrow = (p - 1), ncol = 1))
    ar.part <- rbind(ar.coefs, ident)
    ar.rows <- cbind(ar.rows, ar.part)
    if (!is.null(ma.coefs)) {
      ma.in.ar <- matrix(0, nrow = p, ncol = q)
      ma.in.ar[1, ] <- ma.coefs
      ar.rows <- cbind(ar.rows, ma.in.ar)
    }
    F <- rbind(F, ar.rows)
  }
  if (!is.null(ma.coefs)) {
    ma.rows <- matrix(0, nrow = q, ncol = 1)
    if (!is.null(beta)) {
      ma.rows <- cbind(ma.rows, ma.rows)
    }
    if (!is.null(seasonal.periods)) {
      ma.seasonal <- matrix(0, nrow = q, ncol = tau)
      ma.rows <- cbind(ma.rows, ma.seasonal)
    }
    if (!is.null(ar.coefs)) {
      ar.in.ma <- matrix(0, nrow = q, ncol = p)
      ma.rows <- cbind(ma.rows, ar.in.ma)
    }
    ident <- diag((q - 1))
    ident <- cbind(ident, matrix(0, nrow = (q - 1), ncol = 1))
    ma.part <- rbind(matrix(0, nrow = 1, ncol = q), ident)
    ma.rows <- cbind(ma.rows, ma.part)
    F <- rbind(F, ma.rows)
  }
  return(F)
}

calcModel <- function (y, x.nought, F, g, w) 
{
  length.ts <- length(y)
  x <- matrix(0, nrow = length(x.nought), ncol = length.ts)
  y.hat <- matrix(0, nrow = 1, ncol = length.ts)
  e <- matrix(0, nrow = 1, ncol = length.ts)
  y.hat[, 1] <- w$w.transpose %*% x.nought
  e[, 1] <- y[1] - y.hat[, 1]
  x[, 1] <- F %*% x.nought + g %*% e[, 1]
  y <- matrix(y, nrow = 1, ncol = length.ts)
  loop <- .Call("calcBATS", ys = y, yHats = y.hat, wTransposes = w$w.transpose, 
                Fs = F, xs = x, gs = g, es = e, PACKAGE = "fable.tbats")
  return(list(y.hat = loop$y.hat, e = loop$e, x = loop$x))
}

calcLikelihood <- function (param.vector, opt.env, use.beta, use.small.phi, seasonal.periods, 
                            p = 0, q = 0, tau = 0, bc.lower = 0, bc.upper = 1) 
{
  box.cox.parameter <- param.vector[1]
  alpha <- param.vector[2]
  if (use.beta) {
    if (use.small.phi) {
      small.phi <- param.vector[3]
      beta.v <- param.vector[4]
      gamma.start <- 5
    }
    else {
      small.phi <- 1
      beta.v <- param.vector[3]
      gamma.start <- 4
    }
  }
  else {
    small.phi <- NULL
    beta.v <- NULL
    gamma.start <- 3
  }
  if (!is.null(seasonal.periods)) {
    gamma.vector <- param.vector[gamma.start:(gamma.start + 
                                                length(seasonal.periods) - 1)]
    final.gamma.pos <- gamma.start + length(gamma.vector) - 
      1
  }
  else {
    gamma.vector <- NULL
    final.gamma.pos <- gamma.start - 1
  }
  if (p != 0) {
    ar.coefs <- matrix(param.vector[(final.gamma.pos + 1):(final.gamma.pos + 
                                                             p)], nrow = 1, ncol = p)
  }
  else {
    ar.coefs <- NULL
  }
  if (q != 0) {
    ma.coefs <- matrix(param.vector[(final.gamma.pos + p + 
                                       1):length(param.vector)], nrow = 1, ncol = q)
  }
  else {
    ma.coefs <- NULL
  }
  x.nought <- BoxCox(opt.env$x.nought.untransformed, lambda = box.cox.parameter)
  lambda <- attr(x.nought, "lambda")
  .Call("updateWtransposeMatrix", wTranspose_s = opt.env$w.transpose, 
        smallPhi_s = small.phi, tau_s = as.integer(tau), arCoefs_s = ar.coefs, 
        maCoefs_s = ma.coefs, p_s = as.integer(p), q_s = as.integer(q), 
        PACKAGE = "fable.tbats")
  .Call("updateGMatrix", g_s = opt.env$g, gammaBold_s = opt.env$gamma.bold.matrix, 
        alpha_s = alpha, beta_s = beta.v, gammaVector_s = gamma.vector, 
        seasonalPeriods_s = seasonal.periods, PACKAGE = "fable.tbats")
  .Call("updateFMatrix", opt.env$F, small.phi, alpha, beta.v, 
        opt.env$gamma.bold.matrix, ar.coefs, ma.coefs, tau, PACKAGE = "fable.tbats")
  mat.transformed.y <- BoxCox(opt.env$y, box.cox.parameter)
  lambda <- attr(mat.transformed.y, "lambda")
  n <- ncol(opt.env$y)
  .Call("calcBATSFaster", ys = mat.transformed.y, yHats = opt.env$y.hat, 
        wTransposes = opt.env$w.transpose, Fs = opt.env$F, xs = opt.env$x, 
        gs = opt.env$g, es = opt.env$e, xNought_s = x.nought, 
        sPeriods_s = seasonal.periods, betaV = beta.v, tau_s = as.integer(tau), 
        p_s = as.integer(p), q_s = as.integer(q), PACKAGE = "fable.tbats")
  log.likelihood <- n * log(sum(opt.env$e^2)) - 2 * (box.cox.parameter - 
                                                       1) * sum(log(opt.env$y))
  assign("D", (opt.env$F - opt.env$g %*% opt.env$w.transpose), 
         envir = opt.env)
  if (checkAdmissibility(opt.env, box.cox = box.cox.parameter, 
                         small.phi = small.phi, ar.coefs = ar.coefs, ma.coefs = ma.coefs, 
                         tau = tau, bc.lower = bc.lower, bc.upper = bc.upper)) {
    return(log.likelihood)
  }
  else {
    return(10^20)
  }
}

checkAdmissibility <- function (opt.env, box.cox = NULL, small.phi = NULL, ar.coefs = NULL, 
                                ma.coefs = NULL, tau = 0, bc.lower = 0, bc.upper = 1) 
{
  if (!is.null(box.cox)) {
    if ((box.cox <= bc.lower) | (box.cox >= bc.upper)) {
      return(FALSE)
    }
  }
  if (!is.null(small.phi)) {
    if (((small.phi < 0.8) | (small.phi > 1))) {
      return(FALSE)
    }
  }
  if (!is.null(ar.coefs)) {
    arlags <- which(abs(ar.coefs) > 1e-08)
    if (length(arlags) > 0L) {
      p <- max(arlags)
      if (min(Mod(polyroot(c(1, -ar.coefs[1L:p])))) < 1 + 
          0.01) {
        return(FALSE)
      }
    }
  }
  if (!is.null(ma.coefs)) {
    malags <- which(abs(ma.coefs) > 1e-08)
    if (length(malags) > 0L) {
      q <- max(malags)
      if (min(Mod(polyroot(c(1, ma.coefs[1L:q])))) < 1 + 
          0.01) {
        return(FALSE)
      }
    }
  }
  D.eigen.values <- eigen(opt.env$D, symmetric = FALSE, only.values = TRUE)$values
  return(all(abs(D.eigen.values) < 1 + 0.01))
}

unParameterise <- function (param.vector, control) 
{
  if (control$use.box.cox) {
    lambda <- param.vector[1]
    alpha <- param.vector[2]
    if (control$use.beta) {
      if (control$use.damping) {
        small.phi <- param.vector[3]
        beta <- param.vector[4]
        gamma.start <- 5
      }
      else {
        small.phi <- 1
        beta <- param.vector[3]
        gamma.start <- 4
      }
    }
    else {
      small.phi <- NULL
      beta <- NULL
      gamma.start <- 3
    }
    if (control$length.gamma > 0) {
      gamma.vector <- param.vector[gamma.start:(gamma.start + 
                                                  control$length.gamma - 1)]
      final.gamma.pos <- gamma.start + control$length.gamma - 
        1
    }
    else {
      gamma.vector <- NULL
      final.gamma.pos <- gamma.start - 1
    }
    if (control$p != 0) {
      ar.coefs <- param.vector[(final.gamma.pos + 1):(final.gamma.pos + 
                                                        control$p)]
    }
    else {
      ar.coefs <- NULL
    }
    if (control$q != 0) {
      ma.coefs <- param.vector[(final.gamma.pos + control$p + 
                                  1):length(param.vector)]
    }
    else {
      ma.coefs <- NULL
    }
  }
  else {
    lambda <- NULL
    alpha <- param.vector[1]
    if (control$use.beta) {
      if (control$use.damping) {
        small.phi <- param.vector[2]
        beta <- param.vector[3]
        gamma.start <- 4
      }
      else {
        small.phi <- 1
        beta <- param.vector[2]
        gamma.start <- 3
      }
    }
    else {
      small.phi <- NULL
      beta <- NULL
      gamma.start <- 2
    }
    if (control$length.gamma > 0) {
      gamma.vector <- param.vector[gamma.start:(gamma.start + 
                                                  control$length.gamma - 1)]
      final.gamma.pos <- gamma.start + control$length.gamma - 
        1
    }
    else {
      gamma.vector <- NULL
      final.gamma.pos <- gamma.start - 1
    }
    if (control$p != 0) {
      ar.coefs <- param.vector[(final.gamma.pos + 1):(final.gamma.pos + 
                                                        control$p)]
    }
    else {
      ar.coefs <- NULL
    }
    if (control$q != 0) {
      ma.coefs <- param.vector[(final.gamma.pos + control$p + 
                                  1):length(param.vector)]
    }
    else {
      ma.coefs <- NULL
    }
  }
  return(list(lambda = lambda, alpha = alpha, beta = beta, 
              small.phi = small.phi, gamma.v = gamma.vector, ar.coefs = ar.coefs, 
              ma.coefs = ma.coefs))
}

calcLikelihoodNOTransformed <- function (param.vector, opt.env, x.nought, use.beta, use.small.phi, 
                                         seasonal.periods, p = 0, q = 0, tau = 0) 
{
  alpha <- param.vector[1]
  if (use.beta) {
    if (use.small.phi) {
      small.phi <- param.vector[2]
      beta.v <- param.vector[3]
      gamma.start <- 4
    }
    else {
      small.phi <- 1
      beta.v <- param.vector[2]
      gamma.start <- 3
    }
  }
  else {
    small.phi <- NULL
    beta.v <- NULL
    gamma.start <- 2
  }
  if (!is.null(seasonal.periods)) {
    gamma.vector <- param.vector[gamma.start:(gamma.start + 
                                                length(seasonal.periods) - 1)]
    final.gamma.pos <- gamma.start + length(gamma.vector) - 
      1
  }
  else {
    gamma.vector <- NULL
    final.gamma.pos <- gamma.start - 1
  }
  if (p != 0) {
    ar.coefs <- matrix(param.vector[(final.gamma.pos + 1):(final.gamma.pos + 
                                                             p)], nrow = 1, ncol = p)
  }
  else {
    ar.coefs <- NULL
  }
  if (q != 0) {
    ma.coefs <- matrix(param.vector[(final.gamma.pos + p + 
                                       1):length(param.vector)], nrow = 1, ncol = q)
  }
  else {
    ma.coefs <- NULL
  }
  .Call("updateWtransposeMatrix", wTranspose_s = opt.env$w.transpose, 
        smallPhi_s = small.phi, tau_s = as.integer(tau), arCoefs_s = ar.coefs, 
        maCoefs_s = ma.coefs, p_s = as.integer(p), q_s = as.integer(q), 
        PACKAGE = "fable.tbats")
  .Call("updateGMatrix", g_s = opt.env$g, gammaBold_s = opt.env$gamma.bold.matrix, 
        alpha_s = alpha, beta_s = beta.v, gammaVector_s = gamma.vector, 
        seasonalPeriods_s = seasonal.periods, PACKAGE = "fable.tbats")
  .Call("updateFMatrix", opt.env$F, small.phi, alpha, beta.v, 
        opt.env$gamma.bold.matrix, ar.coefs, ma.coefs, tau, PACKAGE = "fable.tbats")
  n <- ncol(opt.env$y)
  .Call("calcBATSFaster", ys = opt.env$y, yHats = opt.env$y.hat, 
        wTransposes = opt.env$w.transpose, Fs = opt.env$F, xs = opt.env$x, 
        gs = opt.env$g, es = opt.env$e, xNought_s = x.nought, 
        sPeriods_s = seasonal.periods, betaV = beta.v, tau_s = as.integer(tau), 
        p_s = as.integer(p), q_s = as.integer(q), PACKAGE = "fable.tbats")
  log.likelihood <- n * log(sum(opt.env$e * opt.env$e))
  assign("D", (opt.env$F - opt.env$g %*% opt.env$w.transpose), 
         envir = opt.env)
  if (checkAdmissibility(opt.env = opt.env, box.cox = NULL, 
                         small.phi = small.phi, ar.coefs = ar.coefs, ma.coefs = ma.coefs, 
                         tau = tau)) {
    return(log.likelihood)
  }
  else {
    return(10^20)
  }
}

findGCD <- function (larger, smaller) 
{
  remainder <- larger%%smaller
  if (remainder != 0) {
    return(findGCD(smaller, remainder))
  }
  else {
    return(smaller)
  }
}

cutW <- function (use.beta, w.tilda.transpose, seasonal.periods, p = 0, 
                  q = 0) 
{
  mask.vector <- numeric(length(seasonal.periods))
  i <- length(seasonal.periods)
  while (i > 1) {
    for (j in 1:(i - 1)) {
      if ((seasonal.periods[i]%%seasonal.periods[j]) == 
          0) {
        mask.vector[j] <- 1
      }
    }
    i <- i - 1
  }
  if (length(seasonal.periods) > 1) {
    for (s in length(seasonal.periods):2) {
      for (j in (s - 1):1) {
        hcf <- findGCD(seasonal.periods[s], seasonal.periods[j])
        if (hcf != 1) {
          if ((mask.vector[s] != 1) && (mask.vector[j] != 
                                        1)) {
            mask.vector[s] <- hcf * -1
          }
        }
      }
    }
  }
  w.pos.counter <- 1
  w.pos <- 1
  if (use.beta) {
    w.pos <- w.pos + 1
  }
  for (s in seasonal.periods) {
    if (mask.vector[w.pos.counter] == 1) {
      w.tilda.transpose <- w.tilda.transpose[, -((w.pos + 
                                                    1):(w.pos + s))]
    }
    else if (mask.vector[w.pos.counter] < 0) {
      w.pos <- w.pos + s
      w.tilda.transpose <- w.tilda.transpose[, -c((w.pos + 
                                                     mask.vector[w.pos.counter] + 1):w.pos)]
      w.pos <- w.pos + mask.vector[w.pos.counter]
    }
    else {
      w.pos <- w.pos + s
      w.tilda.transpose <- w.tilda.transpose[, -w.pos]
      w.pos <- w.pos - 1
    }
    w.pos.counter <- w.pos.counter + 1
  }
  if ((p != 0) | (q != 0)) {
    end.cut <- ncol(w.tilda.transpose)
    start.cut <- end.cut - (p + q) + 1
    w.tilda.transpose <- w.tilda.transpose[, -c(start.cut:end.cut)]
  }
  return(list(matrix = w.tilda.transpose, mask.vector = mask.vector))
}

calcSeasonalSeeds <- function (use.beta, coefs, seasonal.periods, mask.vector, p = 0, 
                               q = 0) 
{
  x.pos.counter <- 1
  sum.k <- 0
  if (use.beta) {
    x.pos <- 2
    new.x.nought <- matrix(coefs[1:2], nrow = 2, ncol = 1)
  }
  else {
    x.pos <- 1
    new.x.nought <- matrix(coefs[1], nrow = 1, ncol = 1)
  }
  x.pos.counter <- 1
  for (s in seasonal.periods) {
    if (mask.vector[x.pos.counter] == 1) {
      season <- matrix(0, nrow = s, ncol = 1)
      new.x.nought <- rbind(new.x.nought, season)
    }
    else if (mask.vector[x.pos.counter] < 0) {
      extract <- coefs[(x.pos + 1):(x.pos + s + mask.vector[x.pos.counter])]
      k <- sum(extract)
      sum.k <- sum.k + k/s
      current.periodicity <- extract - k/s
      current.periodicity <- matrix(current.periodicity, 
                                    nrow = length(current.periodicity), ncol = 1)
      additional <- matrix(-k/s, nrow = (-1 * mask.vector[x.pos.counter]), 
                           ncol = 1)
      current.periodicity <- rbind(current.periodicity, 
                                   additional)
      new.x.nought <- rbind(new.x.nought, current.periodicity)
      x.pos <- x.pos + s + mask.vector[x.pos.counter]
    }
    else {
      k <- sum(coefs[(x.pos + 1):(x.pos + s - 1)])
      sum.k <- sum.k + k/s
      current.periodicity <- coefs[(x.pos + 1):(x.pos + 
                                                  s - 1)] - k/s
      current.periodicity <- c(current.periodicity, -k/s)
      current.periodicity <- matrix(current.periodicity, 
                                    nrow = length(current.periodicity), ncol = 1)
      new.x.nought <- rbind(new.x.nought, current.periodicity)
      x.pos <- x.pos + s - 1
    }
    x.pos.counter <- x.pos.counter + 1
  }
  if ((p != 0) | (q != 0)) {
    arma.seed.states <- numeric((p + q))
    arma.seed.states <- matrix(arma.seed.states, nrow = length(arma.seed.states), 
                               ncol = 1)
    x.nought <- rbind(new.x.nought, arma.seed.states)
  }
  else {
    x.nought <- new.x.nought
  }
  return(x.nought)
}

unParameteriseTBATS <- function (param.vector, control) 
{
  if (control$use.box.cox) {
    lambda <- param.vector[1]
    alpha <- param.vector[2]
    if (control$use.beta) {
      if (control$use.damping) {
        small.phi <- param.vector[3]
        beta <- param.vector[4]
        gamma.start <- 5
      }
      else {
        small.phi <- 1
        beta <- param.vector[3]
        gamma.start <- 4
      }
    }
    else {
      small.phi <- NULL
      beta <- NULL
      gamma.start <- 3
    }
    if (control$length.gamma > 0) {
      gamma.one.vector <- param.vector[gamma.start:(gamma.start + 
                                                      (control$length.gamma/2) - 1)]
      gamma.two.vector <- param.vector[(gamma.start + (control$length.gamma/2)):(gamma.start + 
                                                                                   (control$length.gamma) - 1)]
      final.gamma.pos <- gamma.start + control$length.gamma - 
        1
    }
    else {
      gamma.one.vector <- NULL
      gamma.two.vector <- NULL
      final.gamma.pos <- gamma.start - 1
    }
    if (control$p != 0) {
      ar.coefs <- param.vector[(final.gamma.pos + 1):(final.gamma.pos + 
                                                        control$p)]
    }
    else {
      ar.coefs <- NULL
    }
    if (control$q != 0) {
      ma.coefs <- param.vector[(final.gamma.pos + control$p + 
                                  1):length(param.vector)]
    }
    else {
      ma.coefs <- NULL
    }
  }
  else {
    lambda <- NULL
    alpha <- param.vector[1]
    if (control$use.beta) {
      if (control$use.damping) {
        small.phi <- param.vector[2]
        beta <- param.vector[3]
        gamma.start <- 4
      }
      else {
        small.phi <- 1
        beta <- param.vector[2]
        gamma.start <- 3
      }
    }
    else {
      small.phi <- NULL
      beta <- NULL
      gamma.start <- 2
    }
    if (control$length.gamma > 0) {
      gamma.one.vector <- param.vector[gamma.start:(gamma.start + 
                                                      (control$length.gamma/2) - 1)]
      gamma.two.vector <- param.vector[(gamma.start + (control$length.gamma/2)):(gamma.start + 
                                                                                   (control$length.gamma) - 1)]
      final.gamma.pos <- gamma.start + control$length.gamma - 
        1
    }
    else {
      gamma.one.vector <- NULL
      gamma.two.vector <- NULL
      final.gamma.pos <- gamma.start - 1
    }
    if (control$p != 0) {
      ar.coefs <- param.vector[(final.gamma.pos + 1):(final.gamma.pos + 
                                                        control$p)]
    }
    else {
      ar.coefs <- NULL
    }
    if (control$q != 0) {
      ma.coefs <- param.vector[(final.gamma.pos + control$p + 
                                  1):length(param.vector)]
    }
    else {
      ma.coefs <- NULL
    }
  }
  return(list(lambda = lambda, alpha = alpha, beta = beta, 
              small.phi = small.phi, gamma.one.v = gamma.one.vector, 
              gamma.two.v = gamma.two.vector, ar.coefs = ar.coefs, 
              ma.coefs = ma.coefs))
}

makeTBATSFMatrix <- function (alpha, beta = NULL, small.phi = NULL, seasonal.periods = NULL, 
                              k.vector = NULL, gamma.bold.matrix = NULL, ar.coefs = NULL, 
                              ma.coefs = NULL) 
{
  F <- matrix(1, nrow = 1, ncol = 1)
  if (!is.null(beta)) {
    F <- cbind(F, matrix(small.phi, nrow = 1, ncol = 1))
  }
  if (!is.null(seasonal.periods)) {
    tau <- sum(k.vector) * 2
    zero.tau <- matrix(0, nrow = 1, ncol = tau)
    F <- cbind(F, zero.tau)
  }
  if (!is.null(ar.coefs)) {
    p <- length(ar.coefs)
    ar.coefs <- matrix(ar.coefs, nrow = 1, ncol = p)
    alpha.phi <- alpha * ar.coefs
    F <- cbind(F, alpha.phi)
  }
  if (!is.null(ma.coefs)) {
    q <- length(ma.coefs)
    ma.coefs <- matrix(ma.coefs, nrow = 1, ncol = q)
    alpha.theta <- alpha * ma.coefs
    F <- cbind(F, alpha.theta)
  }
  if (!is.null(beta)) {
    beta.row <- matrix(c(0, small.phi), nrow = 1, ncol = 2)
    if (!is.null(seasonal.periods)) {
      beta.row <- cbind(beta.row, zero.tau)
    }
    if (!is.null(ar.coefs)) {
      beta.phi <- beta * ar.coefs
      beta.row <- cbind(beta.row, beta.phi)
    }
    if (!is.null(ma.coefs)) {
      beta.theta <- beta * ma.coefs
      beta.row <- cbind(beta.row, beta.theta)
    }
    F <- rbind(F, beta.row)
  }
  if (!is.null(seasonal.periods)) {
    seasonal.row <- t(zero.tau)
    if (!is.null(beta)) {
      seasonal.row <- cbind(seasonal.row, seasonal.row)
    }
    A <- matrix(0, tau, tau)
    last.pos <- 0
    for (i in 1:length(k.vector)) {
      if (seasonal.periods[i] != 2) {
        C <- .Call("makeCIMatrix", k_s = as.integer(k.vector[i]), 
                   m_s = as.double(seasonal.periods[i]), PACKAGE = "fable.tbats")
      }
      else {
        C <- matrix(0, 1, 1)
      }
      S <- .Call("makeSIMatrix", k_s = as.integer(k.vector[i]), 
                 m_s = as.double(seasonal.periods[i]), PACKAGE = "fable.tbats")
      Ai <- .Call("makeAIMatrix", C_s = C, S_s = S, k_s = as.integer(k.vector[i]), 
                  PACKAGE = "fable.tbats")
      A[(last.pos + 1):(last.pos + (2 * k.vector[i])), 
        (last.pos + 1):(last.pos + (2 * k.vector[i]))] <- Ai
      last.pos <- last.pos + (2 * k.vector[i])
    }
    seasonal.row <- cbind(seasonal.row, A)
    if (!is.null(ar.coefs)) {
      B <- t(gamma.bold.matrix) %*% ar.coefs
      seasonal.row <- cbind(seasonal.row, B)
    }
    if (!is.null(ma.coefs)) {
      C <- t(gamma.bold.matrix) %*% ma.coefs
      seasonal.row <- cbind(seasonal.row, C)
    }
    F <- rbind(F, seasonal.row)
  }
  if (!is.null(ar.coefs)) {
    ar.rows <- matrix(0, nrow = p, ncol = 1)
    if (!is.null(beta)) {
      ar.rows <- cbind(ar.rows, ar.rows)
    }
    if (!is.null(seasonal.periods)) {
      ar.seasonal.zeros <- matrix(0, nrow = p, ncol = tau)
      ar.rows <- cbind(ar.rows, ar.seasonal.zeros)
    }
    ident <- diag((p - 1))
    ident <- cbind(ident, matrix(0, nrow = (p - 1), ncol = 1))
    ar.part <- rbind(ar.coefs, ident)
    ar.rows <- cbind(ar.rows, ar.part)
    if (!is.null(ma.coefs)) {
      ma.in.ar <- matrix(0, nrow = p, ncol = q)
      ma.in.ar[1, ] <- ma.coefs
      ar.rows <- cbind(ar.rows, ma.in.ar)
    }
    F <- rbind(F, ar.rows)
  }
  if (!is.null(ma.coefs)) {
    ma.rows <- matrix(0, nrow = q, ncol = 1)
    if (!is.null(beta)) {
      ma.rows <- cbind(ma.rows, ma.rows)
    }
    if (!is.null(seasonal.periods)) {
      ma.seasonal <- matrix(0, nrow = q, ncol = tau)
      ma.rows <- cbind(ma.rows, ma.seasonal)
    }
    if (!is.null(ar.coefs)) {
      ar.in.ma <- matrix(0, nrow = q, ncol = p)
      ma.rows <- cbind(ma.rows, ar.in.ma)
    }
    ident <- diag((q - 1))
    ident <- cbind(ident, matrix(0, nrow = (q - 1), ncol = 1))
    ma.part <- rbind(matrix(0, nrow = 1, ncol = q), ident)
    ma.rows <- cbind(ma.rows, ma.part)
    F <- rbind(F, ma.rows)
  }
  return(F)
}

calcLikelihoodTBATS <- function (param.vector, opt.env, use.beta, use.small.phi, seasonal.periods, 
                                 param.control, p = 0, q = 0, tau = 0, bc.lower = 0, bc.upper = 1) 
{
  paramz <- unParameteriseTBATS(param.vector, param.control)
  box.cox.parameter <- paramz$lambda
  alpha <- paramz$alpha
  beta.v <- paramz$beta
  small.phi <- paramz$small.phi
  gamma.one.v <- paramz$gamma.one.v
  gamma.two.v <- paramz$gamma.two.v
  ar.coefs <- paramz$ar.coefs
  ma.coefs <- paramz$ma.coefs
  if (!is.null(paramz$ar.coefs)) {
    p <- length(paramz$ar.coefs)
    ar.coefs <- matrix(paramz$ar.coefs, nrow = 1, ncol = p)
  }
  else {
    ar.coefs <- NULL
    p <- 0
  }
  if (!is.null(paramz$ma.coefs)) {
    q <- length(paramz$ma.coefs)
    ma.coefs <- matrix(paramz$ma.coefs, nrow = 1, ncol = q)
  }
  else {
    ma.coefs <- NULL
    q <- 0
  }
  x.nought <- BoxCox(opt.env$x.nought.untransformed, lambda = box.cox.parameter)
  lambda <- attr(x.nought, "lambda")
  .Call("updateWtransposeMatrix", wTranspose_s = opt.env$w.transpose, 
        smallPhi_s = small.phi, tau_s = as.integer(tau), arCoefs_s = ar.coefs, 
        maCoefs_s = ma.coefs, p_s = as.integer(p), q_s = as.integer(q), 
        PACKAGE = "fable.tbats")
  if (!is.null(opt.env$gamma.bold)) {
    .Call("updateTBATSGammaBold", gammaBold_s = opt.env$gamma.bold, 
          kVector_s = opt.env$k.vector, gammaOne_s = gamma.one.v, 
          gammaTwo_s = gamma.two.v)
  }
  .Call("updateTBATSGMatrix", g_s = opt.env$g, gammaBold_s = opt.env$gamma.bold, 
        alpha_s = alpha, beta_s = beta.v, PACKAGE = "fable.tbats")
  .Call("updateFMatrix", opt.env$F, small.phi, alpha, beta.v, 
        opt.env$gamma.bold, ar.coefs, ma.coefs, tau, PACKAGE = "fable.tbats")
  mat.transformed.y <- BoxCox(opt.env$y, box.cox.parameter)
  lambda <- attr(mat.transformed.y, "lambda")
  n <- ncol(opt.env$y)
  .Call("calcTBATSFaster", ys = mat.transformed.y, yHats = opt.env$y.hat, 
        wTransposes = opt.env$w.transpose, Fs = opt.env$F, xs = opt.env$x, 
        gs = opt.env$g, es = opt.env$e, xNought_s = x.nought, 
        PACKAGE = "fable.tbats")
  log.likelihood <- n * log(sum(opt.env$e^2)) - 2 * (box.cox.parameter - 
                                                       1) * sum(log(opt.env$y))
  if (is.na(log.likelihood)) {
    return(Inf)
  }
  assign("D", (opt.env$F - opt.env$g %*% opt.env$w.transpose), 
         envir = opt.env)
  if (checkAdmissibility(opt.env, box.cox = box.cox.parameter, 
                         small.phi = small.phi, ar.coefs = ar.coefs, ma.coefs = ma.coefs, 
                         tau = sum(seasonal.periods), bc.lower = bc.lower, bc.upper = bc.upper)) {
    return(log.likelihood)
  }
  else {
    return(Inf)
  }
}

calcLikelihoodNOTransformedTBATS <- function (param.vector, opt.env, x.nought, use.beta, use.small.phi, 
                                              seasonal.periods, param.control, p = 0, q = 0, tau = 0) 
{
  paramz <- unParameteriseTBATS(param.vector, param.control)
  box.cox.parameter <- paramz$lambda
  alpha <- paramz$alpha
  beta.v <- paramz$beta
  small.phi <- paramz$small.phi
  gamma.one.v <- paramz$gamma.one.v
  gamma.two.v <- paramz$gamma.two.v
  if (!is.null(paramz$ar.coefs)) {
    p <- length(paramz$ar.coefs)
    ar.coefs <- matrix(paramz$ar.coefs, nrow = 1, ncol = p)
  }
  else {
    ar.coefs <- NULL
    p <- 0
  }
  if (!is.null(paramz$ma.coefs)) {
    q <- length(paramz$ma.coefs)
    ma.coefs <- matrix(paramz$ma.coefs, nrow = 1, ncol = q)
  }
  else {
    ma.coefs <- NULL
    q <- 0
  }
  .Call("updateWtransposeMatrix", wTranspose_s = opt.env$w.transpose, 
        smallPhi_s = small.phi, tau_s = as.integer(tau), arCoefs_s = ar.coefs, 
        maCoefs_s = ma.coefs, p_s = as.integer(p), q_s = as.integer(q), 
        PACKAGE = "fable.tbats")
  if (!is.null(opt.env$gamma.bold)) {
    .Call("updateTBATSGammaBold", gammaBold_s = opt.env$gamma.bold, 
          kVector_s = opt.env$k.vector, gammaOne_s = gamma.one.v, 
          gammaTwo_s = gamma.two.v)
  }
  .Call("updateTBATSGMatrix", g_s = opt.env$g, gammaBold_s = opt.env$gamma.bold, 
        alpha_s = alpha, beta_s = beta.v, PACKAGE = "fable.tbats")
  .Call("updateFMatrix", opt.env$F, small.phi, alpha, beta.v, 
        opt.env$gamma.bold, ar.coefs, ma.coefs, tau, PACKAGE = "fable.tbats")
  n <- ncol(opt.env$y)
  .Call("calcTBATSFaster", ys = opt.env$y, yHats = opt.env$y.hat, 
        wTransposes = opt.env$w.transpose, Fs = opt.env$F, xs = opt.env$x, 
        gs = opt.env$g, es = opt.env$e, xNought_s = x.nought, 
        PACKAGE = "fable.tbats")
  log.likelihood <- n * log(sum(opt.env$e * opt.env$e))
  if (is.na(log.likelihood)) {
    return(Inf)
  }
  assign("D", (opt.env$F - opt.env$g %*% opt.env$w.transpose), 
         envir = opt.env)
  if (checkAdmissibility(opt.env = opt.env, box.cox = NULL, 
                         small.phi = small.phi, ar.coefs = ar.coefs, ma.coefs = ma.coefs, 
                         tau = tau)) {
    return(log.likelihood)
  }
  else {
    return(Inf)
  }
}

parFitSpecificTBATS <- function (control.number, y, box.cox, trend, damping, seasonal.periods, 
                                 k.control.matrix, init.box.cox = NULL, bc.lower = 0, bc.upper = 1, 
                                 biasadj = FALSE) 
{
  k.vector <- k.control.matrix[control.number, ]
  model <- try(fitSpecificTBATS(y, use.box.cox = box.cox, use.beta = trend, 
                                use.damping = damping, seasonal.periods = seasonal.periods, 
                                k.vector = k.vector, init.box.cox = init.box.cox, bc.lower = bc.lower, 
                                bc.upper = bc.upper, biasadj = biasadj), silent = TRUE)
  if (is.element("try-error", class(model))) {
    model <- list(AIC = Inf)
  }
  return(model)
}

parFilterSpecifics <- function (control.number, control.array, y, seasonal.periods, 
                                use.arma.errors, force.seasonality = FALSE, init.box.cox = NULL, 
                                bc.lower = 0, bc.upper = 1, biasadj = FALSE, ...) 
{
  box.cox <- control.array[control.number, 1]
  trend <- control.array[control.number, 2]
  damping <- control.array[control.number, 3]
  if (!trend && damping) {
    return(list(AIC = Inf))
  }
  first.model <- fitSpecificBATS(y, use.box.cox = box.cox, 
                                 use.beta = trend, use.damping = damping, seasonal.periods = seasonal.periods, 
                                 init.box.cox = init.box.cox, bc.lower = bc.lower, bc.upper = bc.upper, 
                                 biasadj = biasadj)
  if (!is.null(seasonal.periods) && !force.seasonality) {
    non.seasonal.model <- fitSpecificBATS(y, use.box.cox = box.cox, 
                                          use.beta = trend, use.damping = damping, seasonal.periods = NULL, 
                                          init.box.cox = init.box.cox, bc.lower = bc.lower, 
                                          bc.upper = bc.upper, biasadj = biasadj)
    if (first.model$AIC > non.seasonal.model$AIC) {
      seasonal.periods <- NULL
      first.model <- non.seasonal.model
    }
  }
  if (use.arma.errors) {
    suppressWarnings(arma <- auto_arma(as.numeric(first.model$errors), 
                                       d = 0, ...))
    p <- arma$arma[1]
    q <- arma$arma[2]
    if (p != 0 || q != 0) {
      if (p != 0) {
        ar.coefs <- numeric(p)
      }
      else {
        ar.coefs <- NULL
      }
      if (q != 0) {
        ma.coefs <- numeric(q)
      }
      else {
        ma.coefs <- NULL
      }
      starting.params <- first.model$parameters
      second.model <- fitSpecificBATS(y, use.box.cox = box.cox, 
                                      use.beta = trend, use.damping = damping, seasonal.periods = seasonal.periods, 
                                      ar.coefs = ar.coefs, ma.coefs = ma.coefs, init.box.cox = init.box.cox, 
                                      bc.lower = bc.lower, bc.upper = bc.upper, biasadj = biasadj)
      if (second.model$AIC < first.model$AIC) {
        return(second.model)
      }
      else {
        return(first.model)
      }
    }
    else {
      return(first.model)
    }
  }
  else {
    return(first.model)
  }
}

parFilterTBATSSpecifics <- function (control.number, y, control.array, model.params, seasonal.periods, 
                                     k.vector, use.arma.errors, aux.model = NULL, init.box.cox = NULL, 
                                     bc.lower = 0, bc.upper = 1, biasadj = FALSE, ...) 
{
  box.cox <- control.array[control.number, 1]
  trend <- control.array[control.number, 2]
  damping <- control.array[control.number, 3]
  if (!all((model.params == c(box.cox, trend, damping)))) {
    first.model <- try(fitSpecificTBATS(y, use.box.cox = box.cox, 
                                        use.beta = trend, use.damping = damping, seasonal.periods = seasonal.periods, 
                                        k.vector = k.vector, init.box.cox = init.box.cox, 
                                        bc.lower = bc.lower, bc.upper = bc.upper, biasadj = biasadj), 
                       silent = TRUE)
  }
  else {
    first.model <- aux.model
  }
  if (is.element("try-error", class(first.model))) {
    first.model <- list(AIC = Inf)
  }
  if (use.arma.errors) {
    suppressWarnings(arma <- try(auto_arma(as.numeric(first.model$errors), 
                                           d = 0, ...), silent = TRUE))
    if (!is.element("try-error", class(arma))) {
      p <- arma$arma[1]
      q <- arma$arma[2]
      if ((p != 0) || (q != 0)) {
        if (p != 0) {
          ar.coefs <- numeric(p)
        }
        else {
          ar.coefs <- NULL
        }
        if (q != 0) {
          ma.coefs <- numeric(q)
        }
        else {
          ma.coefs <- NULL
        }
        starting.params <- first.model$parameters
        second.model <- try(fitSpecificTBATS(y, use.box.cox = box.cox, 
                                             use.beta = trend, use.damping = damping, seasonal.periods = seasonal.periods, 
                                             k.vector = k.vector, ar.coefs = ar.coefs, ma.coefs = ma.coefs, 
                                             init.box.cox = init.box.cox, bc.lower = bc.lower, 
                                             bc.upper = bc.upper, biasadj = biasadj), silent = TRUE)
        if (is.element("try-error", class(second.model))) {
          second.model <- list(AIC = Inf)
        }
        if (second.model$AIC < first.model$AIC) {
          return(second.model)
        }
        else {
          return(first.model)
        }
      }
      else {
        return(first.model)
      }
    }
    else {
      return(first.model)
    }
  }
  else {
    return(first.model)
  }
}

fitPreviousBATSModel <- function (y, model, biasadj = FALSE) 
{
  seasonal.periods <- model$seasonal.periods
  if (is.null(seasonal.periods) == FALSE) {
    seasonal.periods <- as.integer(sort(seasonal.periods))
  }
  paramz <- unParameterise(model$parameters$vect, model$parameters$control)
  lambda <- paramz$lambda
  alpha <- paramz$alpha
  beta.v <- paramz$beta
  small.phi <- paramz$small.phi
  gamma.v <- paramz$gamma.v
  ar.coefs <- paramz$ar.coefs
  ma.coefs <- paramz$ma.coefs
  p <- length(ar.coefs)
  q <- length(ma.coefs)
  w <- .Call("makeBATSWMatrix", smallPhi_s = small.phi, sPeriods_s = seasonal.periods, 
             arCoefs_s = ar.coefs, maCoefs_s = ma.coefs, PACKAGE = "fable.tbats")
  g <- .Call("makeBATSGMatrix", as.numeric(alpha), beta.v, 
             gamma.v, seasonal.periods, as.integer(p), as.integer(q), 
             PACKAGE = "fable.tbats")
  F <- makeFMatrix(alpha = alpha, beta = beta.v, small.phi <- small.phi, 
                   seasonal.periods = seasonal.periods, gamma.bold.matrix = g$gamma.bold.matrix, 
                   ar.coefs = ar.coefs, ma.coefs = ma.coefs)
  y.touse <- y
  if (!is.null(lambda)) {
    y.touse <- BoxCox(y, lambda = lambda)
    lambda <- attr(y.touse, "lambda")
  }
  fitted.values.and.errors <- calcModel(y.touse, model$seed.states, 
                                        F, g$g, w)
  e <- fitted.values.and.errors$e
  fitted.values <- fitted.values.and.errors$y.hat
  variance <- sum((e * e))/length(y)
  if (!is.null(lambda)) {
    fitted.values <- InvBoxCox(fitted.values, lambda = lambda, 
                               biasadj, variance)
  }
  model.for.output <- model
  model.for.output$variance <- variance
  model.for.output$fitted.values <- c(fitted.values)
  model.for.output$errors <- c(e)
  model.for.output$x <- fitted.values.and.errors$x
  model.for.output$y <- y
  attributes(model.for.output$fitted.values) <- attributes(model.for.output$errors) <- attributes(y)
  return(model.for.output)
}

fitPreviousTBATSModel <- function (y, model, biasadj = FALSE) 
{
  seasonal.periods <- model$seasonal.periods
  if (is.null(seasonal.periods) == FALSE) {
    seasonal.periods <- sort(seasonal.periods)
  }
  paramz <- unParameteriseTBATS(model$parameters$vect, model$parameters$control)
  lambda <- paramz$lambda
  alpha <- paramz$alpha
  beta.v <- paramz$beta
  if (!is.null(beta.v)) {
    adj.beta <- 1
  }
  else {
    adj.beta <- 0
  }
  small.phi <- paramz$small.phi
  gamma.one.v <- paramz$gamma.one.v
  gamma.two.v <- paramz$gamma.two.v
  if (!is.null(paramz$ar.coefs)) {
    p <- length(paramz$ar.coefs)
    ar.coefs <- matrix(paramz$ar.coefs, nrow = 1, ncol = p)
  }
  else {
    ar.coefs <- NULL
    p <- 0
  }
  if (!is.null(paramz$ma.coefs)) {
    q <- length(paramz$ma.coefs)
    ma.coefs <- matrix(paramz$ma.coefs, nrow = 1, ncol = q)
  }
  else {
    ma.coefs <- NULL
    q <- 0
  }
  if (!is.null(seasonal.periods)) {
    tau <- as.integer(2 * sum(model$k.vector))
    gamma.bold <- matrix(0, nrow = 1, ncol = (2 * sum(model$k.vector)))
  }
  else {
    tau <- as.integer(0)
    gamma.bold <- NULL
  }
  g <- matrix(0, nrow = ((2 * sum(model$k.vector)) + 1 + adj.beta + 
                           p + q), ncol = 1)
  if (p != 0) {
    g[(1 + adj.beta + tau + 1), 1] <- 1
  }
  if (q != 0) {
    g[(1 + adj.beta + tau + p + 1), 1] <- 1
  }
  y.touse <- y
  if (is.null(lambda) == FALSE) {
    y.touse <- BoxCox(y, lambda = lambda)
    lambda <- attr(y.touse, "lambda")
  }
  w <- .Call("makeTBATSWMatrix", smallPhi_s = small.phi, kVector_s = model$k.vector, 
             arCoefs_s = ar.coefs, maCoefs_s = ma.coefs, tau_s = tau, 
             PACKAGE = "fable.tbats")
  if (!is.null(gamma.bold)) {
    .Call("updateTBATSGammaBold", gammaBold_s = gamma.bold, 
          kVector_s = model$k.vector, gammaOne_s = gamma.one.v, 
          gammaTwo_s = gamma.two.v, PACKAGE = "fable.tbats")
  }
  .Call("updateTBATSGMatrix", g_s = g, gammaBold_s = gamma.bold, 
        alpha_s = alpha, beta_s = beta.v, PACKAGE = "fable.tbats")
  F <- makeTBATSFMatrix(alpha = alpha, beta = beta.v, small.phi = small.phi, 
                        seasonal.periods = seasonal.periods, k.vector = model$k.vector, 
                        gamma.bold.matrix = gamma.bold, ar.coefs = ar.coefs, 
                        ma.coefs = ma.coefs)
  .Call("updateFMatrix", F, small.phi, alpha, beta.v, gamma.bold, 
        ar.coefs, ma.coefs, tau, PACKAGE = "fable.tbats")
  fitted.values.and.errors <- calcModel(y.touse, model$seed.states, 
                                        F, g, w)
  e <- fitted.values.and.errors$e
  fitted.values <- fitted.values.and.errors$y.hat
  variance <- sum((e * e))/length(y)
  if (!is.null(lambda)) {
    fitted.values <- InvBoxCox(fitted.values, lambda = lambda, 
                               biasadj, variance)
  }
  model.for.output <- model
  model.for.output$variance <- variance
  model.for.output$fitted.values <- ts(c(fitted.values))
  model.for.output$errors <- ts(c(e))
  tsp(model.for.output$fitted.values) <- tsp(model.for.output$errors) <- tsp(y)
  model.for.output$x <- fitted.values.and.errors$x
  model.for.output$y <- y
  return(model.for.output)
}
#' @exportS3Method 
as.character.tbats <- function (x, ...) 
{
  name <- "TBATS("
  if (!is.null(x$lambda)) {
    name <- paste(name, round(x$lambda, digits = 3), sep = "")
  }
  else {
    name <- paste(name, "1", sep = "")
  }
  name <- paste(name, ", {", sep = "")
  if (!is.null(x$ar.coefficients)) {
    name <- paste(name, length(x$ar.coefficients), sep = "")
  }
  else {
    name <- paste(name, "0", sep = "")
  }
  name <- paste(name, ",", sep = "")
  if (!is.null(x$ma.coefficients)) {
    name <- paste(name, length(x$ma.coefficients), sep = "")
  }
  else {
    name <- paste(name, "0", sep = "")
  }
  name <- paste(name, "}, ", sep = "")
  if (!is.null(x$damping.parameter)) {
    name <- paste(name, round(x$damping.parameter, digits = 3), 
                  ",", sep = "")
  }
  else {
    name <- paste(name, "-,", sep = "")
  }
  if (!is.null(x$seasonal.periods)) {
    name <- paste(name, " {", sep = "")
    M <- length(x$seasonal.periods)
    for (i in 1:M) {
      name <- paste(name, "<", round(x$seasonal.periods[i], 
                                     2), ",", x$k.vector[i], ">", sep = "")
      if (i < M) {
        name <- paste(name, ", ", sep = "")
      }
      else {
        name <- paste(name, "})", sep = "")
      }
    }
  }
  else {
    name <- paste(name, "{-})", sep = "")
  }
  return(name)
}
#' @exportS3Method 
as.character.bats <- function (x, ...) 
{
  name <- "BATS("
  if (!is.null(x$lambda)) {
    name <- paste(name, round(x$lambda, digits = 3), sep = "")
  }
  else {
    name <- paste(name, "1", sep = "")
  }
  name <- paste(name, ", {", sep = "")
  if (!is.null(x$ar.coefficients)) {
    name <- paste(name, length(x$ar.coefficients), sep = "")
  }
  else {
    name <- paste(name, "0", sep = "")
  }
  name <- paste(name, ",", sep = "")
  if (!is.null(x$ma.coefficients)) {
    name <- paste(name, length(x$ma.coefficients), sep = "")
  }
  else {
    name <- paste(name, "0", sep = "")
  }
  name <- paste(name, "}, ", sep = "")
  if (!is.null(x$damping.parameter)) {
    name <- paste(name, round(x$damping.parameter, digits = 3), 
                  sep = "")
  }
  else {
    name <- paste(name, "-", sep = "")
  }
  name <- paste(name, ", ", sep = "")
  if (!is.null(x$seasonal.periods)) {
    name <- paste(name, "{", sep = "")
    for (i in x$seasonal.periods) {
      name <- paste(name, i, sep = "")
      if (i != x$seasonal.periods[length(x$seasonal.periods)]) {
        name <- paste(name, ",", sep = "")
      }
      else {
        name <- paste(name, "})", sep = "")
      }
    }
  }
  else {
    name <- paste(name, "-)", sep = "")
  }
  return(name)
}
# nocov end

