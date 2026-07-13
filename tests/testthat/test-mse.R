# Unit tests — mse.R module
# The translation of measures into modifiers (FM/SL50/SL95) does not
# depend on any external package; the LBSPR engine is tested elsewhere (skip).

test_that("mse_measure creates an object of the right class", {

  m <- mse_measure("Statu quo", "statuquo")
  expect_s3_class(m, "mse_measure")
  expect_equal(m$type, "statuquo")

})


test_that(".apply_measure: effort reduction reduces F/M", {

  m   <- mse_measure("-30%", "effort", effort_reduction = 0.3)
  mod <- .apply_measure(m, base_FM = 1.5, base_SL50 = 20, base_SL95 = 22)

  expect_equal(mod$FM, 1.5 * 0.7, tolerance = 1e-8)
  expect_equal(mod$SL50, 20)          # selectivity unchanged

})


test_that(".apply_measure: biological rest = proportional reduction", {

  m   <- mse_measure("Rest 3 months", "repos_biologique", closure_months = 3)
  mod <- .apply_measure(m, base_FM = 2, base_SL50 = 20, base_SL95 = 22)

  expect_equal(mod$FM, 2 * (1 - 3/12), tolerance = 1e-8)

})


test_that(".apply_measure: MPA reduces F according to the protected fraction", {

  m   <- mse_measure("MPA 20%", "amp", mpa_fraction = 0.2)
  mod <- .apply_measure(m, base_FM = 1, base_SL50 = 20, base_SL95 = 22)

  expect_equal(mod$FM, 0.8, tolerance = 1e-8)

})


test_that(".apply_measure: minimum size modifies selectivity", {

  m   <- mse_measure("Min size 25", "taille_min", Lc = 25)
  mod <- .apply_measure(m, base_FM = 1.5, base_SL50 = 20, base_SL95 = 22)

  expect_equal(mod$SL50, 25)
  expect_equal(mod$SL95, 25 * 1.1, tolerance = 1e-8)
  expect_equal(mod$FM, 1.5)           # effort unchanged

})


test_that(".apply_measure: taille_min without Lc raises an error", {

  m <- mse_measure("Without Lc", "taille_min")
  expect_error(.apply_measure(m, 1.5, 20, 22))

})

