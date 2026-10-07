# Frozen behaviour of the R and Armadillo implementation (fable.tbats 0.5.0);
# see fixtures/make_fixtures.R and src/tbats/DESIGN.md, section 8.

# pkgload::load_all() compiles with -O0, which makes the search tests about
# six times slower than R CMD check's -O2 build; run
# pkgbuild::compile_dll(debug = FALSE) first to test at full speed.
fixture_dir <- testthat::test_path("fixtures")

load_fixture <- function(name) {
  readRDS(file.path(fixture_dir, paste0(name, ".rds")))
}

fixture_series <- c("lynx", "usaccdeaths", "multiseasonal", "nonpositive")

# every non-degenerate fitted model held by a fixture, named
fixture_models <- function(fixture) {
  models <- c(list(tbats = fixture$tbats, bats = fixture$bats), fixture$fixed)
  Filter(function(m) !isTRUE(m$degenerate), models)
}

# the specification list the C++ glue expects, from a fixture snapshot
fixture_spec <- function(m) {
  control <- m$parameters$control
  trigonometric <- m$class[1] == "tbats"
  list(
    box_cox = isTRUE(control$use.box.cox),
    box_cox_lower = 0, box_cox_upper = 1,
    trend = isTRUE(control$use.beta),
    damping = isTRUE(control$use.damping),
    seasonal_type = if (trigonometric) "trigonometric" else "dummy",
    periods = as.numeric(m$seasonal_periods),
    harmonics = if (trigonometric) as.numeric(m$k_vector) else numeric(0),
    p = length(m$ar), q = length(m$ma)
  )
}

fixture_parameters <- function(m) {
  trigonometric <- m$class[1] == "tbats"
  list(
    lambda = if (is.null(m$lambda)) 1 else m$lambda,
    alpha = m$alpha,
    beta = if (is.null(m$beta)) 0 else m$beta,
    phi = if (isTRUE(m$parameters$control$use.damping)) m$damping else 1,
    gamma_one = if (trigonometric) as.numeric(m$gamma_one) else as.numeric(m$gamma),
    gamma_two = if (trigonometric) as.numeric(m$gamma_two) else numeric(0),
    ar = as.numeric(m$ar), ma = as.numeric(m$ma)
  )
}

tbats_call <- function(name, ...) {
  do.call(utils::getFromNamespace(name, "fable.tbats"), list(...))
}
