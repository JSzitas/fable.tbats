library(dplyr)

pelt <- tsibbledata::pelt
train <- pelt %>%
  dplyr::filter( Year < 1930 )
test <- pelt %>%
  dplyr::filter( Year >= 1930 )
model <- fabletools::model(pelt,  tbats = TBATS(Lynx) )

test_that("Utilities for TBATS work", {
  # residuals
  expect_equal( sum(residuals(model[[1]][[1]][["fit"]])),
                81778.24,
                tolerance = 0.05
                )
  # fitted
  expect_equal( sum(fitted(model[[1]][[1]][["fit"]])),
                2464540,
                tolerance = 0.05
  )
})

test_that( "Forecasts for TBATS work", {
  fcsts <- fabletools::forecast(model, h = 3)$.mean
  expect_equal( fcsts,
                c(41090.52, 41973.19, 38789.32),
                tolerance = 0.05)
})

test_that( "Refitting a TBATS works", {

  model <- fabletools::model(train,  tbats = TBATS(Lynx) )
  expect_equal( as.character(model[[1]][[1]][["fit"]][["fit"]]),
                "BATS(0.15, {3,2}, -, -)")

  model <- refit(model, pelt)
  expect_equal( as.character(model[[1]][[1]][["fit"]][["fit"]]),
                "BATS(0.168, {3,2}, 1, -)")
})
