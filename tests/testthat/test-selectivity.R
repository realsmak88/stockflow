test_that("fit_gillnet_selectivity fits the Millar model (TropFishR)", {

  skip_if_not_installed("TropFishR")

  data(gillnet, package = "TropFishR")
  sel <- fit_gillnet_selectivity(gillnet, x0 = c(60, 4), plot = FALSE)

  expect_s3_class(sel, "FishStockSelectivity")
  expect_equal(sel$gear, "gillnet")
  ## the estimated modal size is finite and within the size range
  mode1 <- sel$estimates[1, 1]
  expect_true(is.finite(mode1))
  expect_true(mode1 > min(gillnet$midLengths) &&
              mode1 < max(gillnet$midLengths))

})


test_that("gillnet_mesh_for_Lc inverts the mesh-size relationship", {

  skip_if_not_installed("TropFishR")

  data(gillnet, package = "TropFishR")
  sel <- fit_gillnet_selectivity(gillnet, x0 = c(60, 4), plot = FALSE)

  m <- gillnet_mesh_for_Lc(sel, Lc = c(50, 60))
  expect_equal(nrow(m), 2)
  ## mesh size increases with target size
  expect_true(m$mesh[2] > m$mesh[1])
  ## proportionality: Lc doubles the mesh at constant k
  expect_equal(m$mesh[2] / m$mesh[1], 60 / 50, tolerance = 1e-6)

})


test_that("retention_curve produces an increasing ogive (TropFishR)", {

  skip_if_not_installed("TropFishR")

  rc <- retention_curve(Lt = seq(20, 70, 5), L50 = 40, L75 = 45)
  expect_s3_class(rc, "data.frame")
  ## retention at L50 = 0.5
  expect_equal(rc$retention[rc$Lt == 40], 0.5, tolerance = 1e-6)
  ## retention at L75 = 0.75
  expect_equal(rc$retention[rc$Lt == 45], 0.75, tolerance = 1e-6)
  ## monotonically increasing
  expect_true(all(diff(rc$retention) >= -1e-9))

  ## knife_edge
  ke <- retention_curve(Lt = c(30, 45), type = "knife_edge", Lc = 40)
  expect_equal(ke$retention, c(0, 1))

  ## missing arguments -> error
  expect_error(retention_curve(Lt = 1:5, L50 = 3), "requires L50 and L75")

})


test_that("recommend_sizes ranks sizes vs biological references", {

  adv <- recommend_sizes(candidates = c(10, 13, 15), Linf = 42,
                         L50_maturity = 13, MK = 1.4)
  expect_s3_class(adv, "FishStockSizeAdvice")

  ## Lc = 10 < L50 -> does not protect reproduction
  expect_false(adv$advice$protege_reproduction[adv$advice$Lc == 10])
  ## Lc = 13 = L50 -> protects
  expect_true(adv$advice$protege_reproduction[adv$advice$Lc == 13])
  ## "insufficient" verdict for Lc < L50
  expect_match(adv$advice$verdict[adv$advice$Lc == 10], "Insufficient")

  ## Lopt computed by lbb_reference_points: 3/(3+MK)*Linf
  expect_equal(adv$reference$Lopt, 3 / (3 + 1.4) * 42, tolerance = 1e-6)

  ## invalid inputs
  expect_error(recommend_sizes(candidates = numeric(0), Linf = 42),
               "non-empty numeric vector")
  expect_error(recommend_sizes(candidates = 10, Linf = -1),
               "positive number")

})