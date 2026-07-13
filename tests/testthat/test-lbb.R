# Unit tests — lbb.R module
# lbb_reference_points and prepare_lbb_data do not depend on the LBB package.

test_that("lbb_reference_points applies Froese's formulas", {

  rp <- lbb_reference_points(Linf = 45, MK = 1.5, FM = 1)

  # Lopt/Linf = 3 / (3 + M/K)
  expect_equal(rp$Lopt_Linf, 3 / (3 + 1.5), tolerance = 1e-8)
  expect_equal(rp$Lopt, 45 * 3 / 4.5, tolerance = 1e-8)

  # Lc_opt/Linf = (2 + 3 F/M) / ((1 + F/M)(3 + M/K))
  expect_equal(rp$Lc_opt_Linf,
               (2 + 3 * 1) / ((1 + 1) * (3 + 1.5)),
               tolerance = 1e-8)

  # plausible bounds
  expect_true(rp$Lopt_Linf > 0 && rp$Lopt_Linf < 1)
  expect_true(rp$Lc_opt_Linf > 0 && rp$Lc_opt_Linf < 1)

})


test_that("prepare_lbb_data correctly aggregates frequencies", {

  x <- data.frame(
    annee  = c(2020, 2020, 2020, 2021),
    LCT    = c(20, 20, 22, 20),
    n      = c(3, 2, 4, 5)
  )

  d <- prepare_lbb_data(x, "LCT", "annee", "n", stock = "test")

  expect_true(all(c("Stock", "Year", "Length", "CatchNo") %in% names(d)))

  # 2020 / length 20 : 3 + 2 = 5
  v <- d$CatchNo[d$Year == 2020 & d$Length == 20]
  expect_equal(v, 5)

})


test_that("prepare_lbb_data counts 1 per row if freq_col is absent", {

  x <- data.frame(annee = c(2020, 2020, 2020),
                  L = c(10, 10, 12))
  d <- prepare_lbb_data(x, "L", "annee")

  expect_equal(d$CatchNo[d$Length == 10], 2)
  expect_equal(d$CatchNo[d$Length == 12], 1)

})
