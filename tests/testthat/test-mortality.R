# Unit tests — module mortality.R
# Runnable without TropFishR (F/E arithmetic); parts depending
# on TropFishR are skipped if the package is not installed.

test_that("estimate_FE correctly computes F, E and F/M", {

  res <- estimate_FE(Z = 1.2, M = 0.4, E_target = 0.5)

  expect_equal(res$F, 0.8)
  expect_equal(res$E, 0.8 / 1.2, tolerance = 1e-8)
  expect_equal(res$FM, 0.8 / 0.4, tolerance = 1e-8)
  expect_true(grepl("Overexploitation", res$statut))

})


test_that("calculate_F_E (alias) remains compatible", {

  out <- calculate_F_E(Z = 1.0, M = 0.5)

  expect_true(all(c("Z", "M", "F", "E") %in% names(out)))
  expect_equal(out$F, 0.5)
  expect_equal(out$E, 0.5)

})


test_that("estimate_FE detects underexploitation", {

  res <- estimate_FE(Z = 0.6, M = 0.5)
  expect_true(grepl("Underexploitation", res$statut))

})


test_that(".growth_pars extracts Linf and K from explicit Linf/K", {

  gp <- .growth_pars(Linf = 80, K = 0.2)
  expect_equal(gp$Linf, 80)
  expect_equal(gp$K, 0.2)

})


test_that(".growth_pars extracts from a $par list", {

  fake <- list(par = list(Linf = 75, K = 0.3))
  gp <- .growth_pars(fake)
  expect_equal(gp$Linf, 75)
  expect_equal(gp$K, 0.3)

})


test_that(".growth_pars extracts from a flat list(Linf, K, t0)", {

  gp <- .growth_pars(list(Linf = 42, K = 0.45, t0 = -0.3))
  expect_equal(gp$Linf, 42)
  expect_equal(gp$K, 0.45)
  ## t_anchor absent from a flat list -> NA (a clean scalar, not numeric(0))
  expect_true(is.na(gp$t_anchor))
  expect_length(gp$t_anchor, 1L)

})


test_that(".growth_pars extracts from a FishStockGrowthAnalysis ($best)", {

  fake <- list(best = list(Linf = 50, K = 0.25, t_anchor = 0.4))
  class(fake) <- "FishStockGrowthAnalysis"
  gp <- .growth_pars(fake)
  expect_equal(gp$Linf, 50)
  expect_equal(gp$K, 0.25)
  expect_equal(gp$t_anchor, 0.4)

})


test_that("estimate_M_all works (if TropFishR is available)", {

  skip_if_not_installed("TropFishR")

  m <- estimate_M_all(Linf = 86.4, K = 0.18, temp = 25)

  expect_s3_class(m, "FishStockM")
  expect_true(nrow(m$table) >= 1)
  expect_true(all(m$table$M > 0))
  expect_true("geomean" %in% names(m$consensus))

})
