# Tests unitaires — module allometry.R

test_that("fit_length_weight recovers a and b on synthetic data", {

  set.seed(1)
  L <- runif(500, 5, 40)
  a_true <- 0.01; b_true <- 3.05
  W <- a_true * L^b_true * exp(stats::rnorm(length(L), 0, 0.03))

  df <- data.frame(L = L, W = W)

  lw <- fit_length_weight(df, "L", "W")

  expect_s3_class(lw, "FishStockLW")
  expect_equal(lw$b, b_true, tolerance = 0.05)
  expect_equal(lw$a, a_true, tolerance = 0.15 * a_true)
  expect_true(lw$r2 > 0.98)

})


test_that("predict.FishStockLW returns consistent weights", {

  df <- data.frame(L = runif(200, 5, 30),
                   W = 0.02 * runif(200, 5, 30)^3)
  lw <- fit_length_weight(df, "L", "W")

  w20 <- predict(lw, length = 20)
  expect_true(is.numeric(w20) && w20 > 0)

})


test_that("fit_length_weight requires a minimum sample size", {

  df <- data.frame(L = 1:5, W = (1:5)^3)
  expect_error(fit_length_weight(df, "L", "W"))

})