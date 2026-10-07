pelt <- tsibbledata::pelt
usa <- tsibble::as_tsibble(USAccDeaths)

test_that("Making a TBATS model works", {
  model <- TBATS(Lynx)
  expect_equal(model[["model"]], "TBATS")
  expect_s3_class(model, "mdl_defn")
  expect_s3_class(model, "R6")
})

test_that("TBATS can be trained", {
  model <- fabletools::model(pelt, tbats = TBATS(Lynx))
  expect_equal(class(model), c("mdl_df", "tbl_df", "tbl", "data.frame"))
  fit <- model$tbats[[1]][["fit"]]
  expect_s3_class(fit, "TBATS")
  # yearly data with no seasonality: TBATS(lambda, {p,q}, phi, {-})
  expect_match(fit[["model_summary"]], "^TBATS\\([0-9.]+, \\{[0-9],[0-9]\\}, (-|[0-9.]+), \\{-\\}\\)$")
  expect_length(fit[["fitted"]], nrow(pelt))
  expect_length(fit[["fit"]][["errors"]], nrow(pelt))
  expect_true(is.finite(fit[["fit"]][["aic"]]))
})

test_that("Passing arguments to TBATS works", {
  model <- fabletools::model(pelt, tbats = TBATS(Lynx ~ parameters(box_cox = FALSE)))
  spec <- model$tbats[[1]][["fit"]][["fit"]][["spec"]]
  expect_false(spec$box_cox)
  expect_match(model$tbats[[1]][["fit"]][["model_summary"]], "^TBATS\\(1, ")

  model <- fabletools::model(pelt, tbats = TBATS(Lynx ~ parameters(box_cox = FALSE, arma_errors = FALSE, trend = TRUE)))
  spec <- model$tbats[[1]][["fit"]][["fit"]][["spec"]]
  expect_equal(spec$p + spec$q, 0)
  expect_true(spec$trend)
  expect_match(model$tbats[[1]][["fit"]][["model_summary"]], "\\{0,0\\}")
})

test_that("TBATS uses the seasonal periods it is given or detects", {
  given <- fabletools::model(usa, tbats = TBATS(value ~ parameters(seasonal_periods = 12)))
  spec <- given$tbats[[1]][["fit"]][["fit"]][["spec"]]
  expect_equal(spec$periods, 12)
  expect_gte(spec$harmonics, 1)
  expect_match(given$tbats[[1]][["fit"]][["model_summary"]], "<12,[0-9]+>")

  detected <- fabletools::model(usa, tbats = TBATS(value))
  expect_equal(detected$tbats[[1]][["fit"]][["fit"]][["spec"]][["periods"]], 12)
})

test_that("missing values are rejected", {
  with_na <- pelt
  with_na$Lynx[5] <- NA
  expect_warning(fabletools::model(with_na, tbats = TBATS(Lynx)), "missing values")
})
