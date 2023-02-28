library(dplyr)

pelt <- tsibbledata::pelt
train <- pelt %>%
  dplyr::filter( Year < 1930 )
test <- pelt %>%
  dplyr::filter( Year >= 1930 )
model <- fabletools::model(pelt,  tbats = BATS(Lynx) )

test_that("Utilities for BATS work", {
  # residuals
  expect_equal( sum(residuals(model[[1]][[1]][["fit"]])),
                115493.1,
                tolerance = 0.05
  )
  # fitted
  expect_equal( sum(fitted(model[[1]][[1]][["fit"]])),
                2463137,
                tolerance = 0.05
  )
})

test_that( "Forecasts for BATS work", {
  fcsts <- fabletools::forecast(model, h = 3)$.mean
  expect_equal( fcsts,
                c(36004, 31556, 24992),
                tolerance = 0.05)
})

test_that( "Refitting a BATS works", {

  model <- fabletools::model(train,  tbats = BATS(Lynx) )
  expect_equal( as.character(model[[1]][[1]][["fit"]][["fit"]]),
                "BATS(0.251, {2,5}, 1, -)")

  model <- refit(model, pelt)
  expect_equal( as.character(model[[1]][[1]][["fit"]][["fit"]]),
                "BATS(0.43, {3,2}, 0.879, -)")
})
