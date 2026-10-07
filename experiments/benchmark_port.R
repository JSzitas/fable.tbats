# Like-for-like timing of the C++ port against the vendored R and Armadillo
# implementation: the same series, the same search, sequential in both
# (use.parallel = FALSE in R; the port is sequential). Compile the package at
# -O2 first (pkgbuild::compile_dll(debug = FALSE)); pkgload's default is -O0.
#
#   Rscript experiments/benchmark_port.R

pkgload::load_all(quiet = TRUE)
fixture_dir <- "tests/testthat/fixtures"

synthetic <- function(n, periods, amplitudes, level, slope, ar, sd, seed) {
  set.seed(seed)
  t <- seq_len(n)
  seasonal <- Reduce(`+`, Map(function(m, a) a * sin(2 * pi * t / m), periods, amplitudes), 0)
  level + slope * t + seasonal + as.numeric(stats::arima.sim(list(ar = ar), n = n, sd = sd))
}

cases <- list(
  list(name = "lynx", y = readRDS(file.path(fixture_dir, "lynx.rds"))$y, periods = NULL, types = c("tbats", "bats")),
  list(name = "usaccdeaths", y = readRDS(file.path(fixture_dir, "usaccdeaths.rds"))$y, periods = 12, types = c("tbats", "bats")),
  list(name = "nonpositive", y = readRDS(file.path(fixture_dir, "nonpositive.rds"))$y, periods = 12, types = c("tbats", "bats")),
  list(name = "multiseasonal", y = readRDS(file.path(fixture_dir, "multiseasonal.rds"))$y, periods = c(7, 30), types = c("tbats", "bats")),
  list(name = "hourly_weekly", y = synthetic(1000, c(24, 168), c(10, 20), level = 200, slope = 0.01, ar = 0.6, sd = 4, seed = 5),
       periods = c(24, 168), types = "tbats"),  # BATS would carry 192 dummy states
  list(name = "long_nonseasonal", y = synthetic(5000, 2, 0, level = 50, slope = 0.001, ar = 0.7, sd = 2, seed = 6),
       periods = NULL, types = c("tbats", "bats"))
)

results <- list()
for (case in cases) {
  for (type in case$types) {
    reference <- if (type == "tbats") tbats else bats
    t_ref <- system.time(m_ref <- reference(case$y, seasonal.periods = case$periods, use.parallel = FALSE))[["elapsed"]]
    t_port <- system.time(m_port <- .Call("tbats_search", case$y, case$periods, type, list(), PACKAGE = "fable.tbats"))[["elapsed"]]
    row <- data.frame(series = case$name, n = length(case$y), periods = paste(case$periods, collapse = ","),
                      search = type, reference_s = t_ref, port_s = t_port, speedup = t_ref / t_port,
                      reference_aic = m_ref$AIC, port_aic = m_port$aic,
                      reference_model = as.character(m_ref))
    results[[length(results) + 1]] <- row
    cat(sprintf("%-17s n %5d %-6s %-5s  R %7.1fs  port %6.1fs  x%5.1f   AIC R %10.2f  port %10.2f\n",
                case$name, length(case$y), row$periods, type, t_ref, t_port, t_ref / t_port, m_ref$AIC, m_port$aic))
  }
}
results <- do.call(rbind, results)
saveRDS(results, "experiments/benchmark_port.rds")
print(results[, c("series", "n", "search", "reference_s", "port_s", "speedup", "reference_aic", "port_aic")], row.names = FALSE)
