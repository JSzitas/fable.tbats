pelt <- tsibbledata::pelt
usa <- tsibble::as_tsibble(USAccDeaths)

test_that("Making a BATS model works", {
  model <- BATS(Lynx)
  expect_equal(model[["model"]], "BATS")
  expect_s3_class(model, "mdl_defn")
  expect_s3_class(model, "R6")
})

test_that("BATS can be trained", {
  model <- fabletools::model(pelt, bats = BATS(Lynx))
  expect_equal(class(model), c("mdl_df", "tbl_df", "tbl", "data.frame"))
  fit <- model$bats[[1]][["fit"]]
  expect_s3_class(fit, "BATS")
  # yearly data with no seasonality: BATS(lambda, {p,q}, phi, -)
  expect_match(fit[["model_summary"]], "^BATS\\([0-9.]+, \\{[0-9],[0-9]\\}, (-|[0-9.]+), -\\)$")
  expect_length(fit[["fitted"]], nrow(pelt))
})

test_that("Passing arguments to BATS works", {
  model <- fabletools::model(pelt, bats = BATS(Lynx ~ parameters(box_cox = FALSE)))
  expect_false(model$bats[[1]][["fit"]][["fit"]][["spec"]][["box_cox"]])
  expect_match(model$bats[[1]][["fit"]][["model_summary"]], "^BATS\\(1, ")

  model <- fabletools::model(pelt, bats = BATS(Lynx ~ parameters(box_cox = FALSE, arma_errors = FALSE)))
  spec <- model$bats[[1]][["fit"]][["fit"]][["spec"]]
  expect_equal(spec$p + spec$q, 0)
})

test_that("BATS uses the period the index implies and needs integer periods", {
  monthly <- fabletools::model(usa, bats = BATS(value))
  spec <- monthly$bats[[1]][["fit"]][["fit"]][["spec"]]
  expect_equal(spec$periods, 12)
  expect_equal(spec$seasonal_type, "dummy")
  expect_match(monthly$bats[[1]][["fit"]][["model_summary"]], "\\{12\\}")
  expect_warning(fabletools::model(usa, bats = BATS(value ~ parameters(seasonal_periods = 12.5))),
                 "integer seasonal periods")
})
