# Tests unitaires — module lbspr.R
# Le coeur LBSPR (S4/LBSPRfit) est ignore si le package n'est pas installe ;
# la logique de statut SPR est testee independamment.

test_that(".lbspr_status correctly classifies SPR vs targets", {

  expect_match(.lbspr_status(0.45, 0.40, 0.20), "Healthy")
  expect_match(.lbspr_status(0.30, 0.40, 0.20), "Concerning")
  expect_match(.lbspr_status(0.10, 0.40, 0.20), "Overexploited")
  expect_true(is.na(.lbspr_status(NA_real_)))

})


test_that(".lbspr_status respects the exact bounds", {

  expect_match(.lbspr_status(0.40, 0.40, 0.20), "Healthy")        # equal to target
  expect_match(.lbspr_status(0.20, 0.40, 0.20), "Concerning") # equal to limit

})


test_that("lbspr_pars validates L95 > L50 (if LBSPR is present)", {

  skip_if_not_installed("LBSPR")

  expect_error(
    lbspr_pars(Linf = 45, L50 = 28, L95 = 22, MK = 1.5)
  )

})


test_that("lbspr_pars derives MK from M and K (if LBSPR is present)", {

  skip_if_not_installed("LBSPR")

  pr <- lbspr_pars(Linf = 45, L50 = 22, L95 = 28, M = 1.5, K = 1.0)
  expect_equal(pr@MK, 1.5, tolerance = 1e-8)

})