test_that("ypr_param assembles parameters and checks inputs", {

  p <- ypr_param(Linf = 130, K = 0.1, M = 0.28, a = 1e-5, b = 3,
                 Lc = 50, Lr = 35)
  expect_s3_class(p, "ypr_param")
  expect_equal(p$Linf, 130)
  expect_equal(p$M, 0.28)
  expect_equal(p$Lc, 50)

  ## Linf/K missing -> error
  expect_error(ypr_param(M = 0.3), "Linf and K are required")
  ## M missing -> error
  expect_error(ypr_param(Linf = 100, K = 0.2),
               "Natural mortality M is required")

})


test_that("run_ypr returns reference points (official TropFishR engine)", {

  skip_if_not_installed("TropFishR")

  data("hake", package = "TropFishR")
  hk <- hake
  hk$Lr <- 35

  yp <- run_ypr(hk, FM_change = seq(0, 3, 0.05), Lc_change = 50,
                curr.E = 0.73, curr.Lc = 50)

  expect_s3_class(yp, "FishStockYPR")
  rp <- yp$reference_points
  expect_true(all(c("F01", "Fmax", "E01", "Emax") %in% names(rp)))
  expect_true(is.finite(rp$F01[1]))
  expect_true(rp$E01[1] > 0 && rp$E01[1] < 1)
  ## F0.1 <= Fmax by construction
  expect_lte(rp$F01[1], rp$Fmax[1] + 1e-8)

})


test_that("compute_tac: ratio_E rule and precaution", {

  ## E target < E current -> TAC < current catch
  tac <- compute_tac(method = "ratio_E", catch_current = 5000,
                     E_current = 0.73, E_target = 0.40, buffer = 0)
  expect_s3_class(tac, "FishStockTAC")
  expect_equal(tac$TAC, 5000 * 0.40 / 0.73, tolerance = 1e-6)

  ## 10% precaution reduces TAC by 10%
  tac2 <- compute_tac(method = "ratio_E", catch_current = 5000,
                      E_current = 0.73, E_target = 0.40, buffer = 0.10)
  expect_equal(tac2$TAC, tac$TAC * 0.9, tolerance = 1e-6)

  ## biomass method
  tacb <- compute_tac(method = "biomass", biomass = 10000, F_target = 0.2)
  expect_equal(tacb$TAC, 2000)

  ## missing inputs -> error
  expect_error(compute_tac(method = "ratio_E", catch_current = 100),
               "requires")
  expect_error(compute_tac(method = "ratio_E", catch_current = 100,
                           E_current = 0, E_target = 0.4),
               "E_current must be")

})


test_that("compute_prorata aggregates and normalizes to 1", {

  df <- data.frame(
    annee = rep(2020:2022, each = 2),
    espece = rep(c("A", "B"), 3),
    capture_t = c(10, 30, 20, 20, 30, 10)
  )
  pr <- compute_prorata(df, segment_col = "espece",
                        catch_col = "capture_t")
  expect_equal(sum(pr$prorata), 1, tolerance = 1e-8)
  ## A total = 60, B total = 60 -> 0.5 / 0.5
  expect_equal(sort(pr$prorata), c(0.5, 0.5))

  ## filtering by year
  pr2 <- compute_prorata(df, segment_col = "espece",
                         catch_col = "capture_t",
                         year_col = "annee", years = 2020)
  ## 2020 : A=10, B=30 -> 0.25 / 0.75
  expect_equal(pr2$prorata[pr2$segment == "A"], 0.25)

})


test_that("allocate_quota distributes the TAC and preserves the sum", {

  pr <- data.frame(segment = c("glans", "pepo", "cymbium"),
                   prorata = c(0.5, 0.3, 0.2))
  q <- allocate_quota(tac = 1000, prorata = pr)
  expect_s3_class(q, "FishStockQuota")
  expect_equal(q$quotas$quota_t[q$quotas$segment == "glans"], 500)
  expect_equal(sum(q$quotas$quota_t), 1000, tolerance = 0.5)

  ## FishStockTAC object accepted as input
  tac <- compute_tac(method = "biomass", biomass = 5000, F_target = 0.2)
  q2 <- allocate_quota(tac = tac, prorata = pr)
  expect_equal(q2$TAC, 1000)

  ## unnormalized prorata -> warning + renormalization
  expect_warning(
    allocate_quota(tac = 1000,
                   prorata = data.frame(segment = c("A", "B"),
                                        prorata = c(2, 2))),
    "do not sum to 1"
  )

  ## named vector accepted
  q3 <- allocate_quota(tac = 1000, prorata = c(A = 0.6, B = 0.4))
  expect_equal(q3$quotas$quota_t[q3$quotas$segment == "A"], 600)

})