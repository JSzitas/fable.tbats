test_that("Making a BATS model works", {

  model <- BATS(Lynx)

  expect_equal( model[["model"]], "BATS")
  expect_s3_class( model, "mdl_defn" )
  expect_s3_class( model, "R6" )
})

test_that("BATS can be trained", {
  pelt <- tsibbledata::pelt

  model <- fabletools::model(pelt,  tbats = BATS(Lynx) )

  expect_equal( class(model), c("mdl_df", "tbl_df", "tbl", "data.frame") )
  expect_equal( as.character(model$tbats[[1]][[1]][["fit"]]),
                "BATS(0.168, {3,2}, 1, -)" )
})

test_that("Passing arguments to BATS works", {

  pelt <- tsibbledata::pelt

  model <- fabletools::model(pelt,  tbats = BATS(Lynx ~ parameters(box_cox = FALSE)) )

  expect_equal( class(model), c("mdl_df", "tbl_df", "tbl", "data.frame") )
  expect_equal( as.character(model$tbats[[1]][[1]][["fit"]]),
                "BATS(1, {2,1}, -, -)" )

  model <- fabletools::model(pelt,  tbats = BATS(Lynx ~ parameters(box_cox = FALSE, arma_errors = FALSE)) )

  expect_equal( class(model), c("mdl_df", "tbl_df", "tbl", "data.frame") )
  expect_equal( as.character(model$tbats[[1]][[1]][["fit"]]),
                "BATS(1, {0,0}, 0.8, -)" )
})
