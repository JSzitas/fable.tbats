test_that("Making a TBATS model works", {

  model <- TBATS(Lynx)

  expect_equal( model[["model"]], "TBATS")
  expect_s3_class( model, "mdl_defn" )
  expect_s3_class( model, "R6" )
})

test_that("TBATS can be trained", {
  pelt <- tsibbledata::pelt

  model <- fabletools::model(pelt,  tbats = TBATS(Lynx) )

  expect_equal( class(model), c("mdl_df", "tbl_df", "tbl", "data.frame") )
  expect_equal( as.character(model$tbats[[1]][[1]][["fit"]]),
                "BATS(0.168, {3,2}, 1, -)" )
})

test_that("Passing arguments to TBATS works", {

  pelt <- tsibbledata::pelt

  model <- fabletools::model(pelt,  tbats = TBATS(Lynx ~ parameters(box_cox = FALSE)) )

  expect_equal( class(model), c("mdl_df", "tbl_df", "tbl", "data.frame") )
  expect_equal( as.character(model$tbats[[1]][[1]][["fit"]]),
                "BATS(1, {2,1}, -, -)" )

  model <- fabletools::model(pelt,  tbats = TBATS(Lynx ~ parameters(box_cox = FALSE, arma_errors = FALSE)) )

  expect_equal( class(model), c("mdl_df", "tbl_df", "tbl", "data.frame") )
  expect_equal( as.character(model$tbats[[1]][[1]][["fit"]]),
                "TBATS(1, {0,0}, 0.8, {<10,2>})" )
})
