# INTEGRATION tests — verify the scientific core on simulated data
# with KNOWN parameters. Skipped on CRAN (computation time) but executed
# by devtools::test().

# ---------------------------------------------------------------------------
# Simulator: annual cohorts following a known von Bertalanffy,
# sampled monthly, with mortality and length-at-age noise.
# ---------------------------------------------------------------------------
simulate_lfq_data <- function(Linf = 45, K = 1.2, t0 = 0,
                              Z = 1.8, cv = 0.08,
                              n_years = 3, recruit_month = 6,
                              n0 = 4000, seed = 42) {

  set.seed(seed)
  out <- list()

  # sampling dates: every month
  samples <- expand.grid(year = seq_len(n_years), month = 1:12)

  for (i in seq_len(nrow(samples))) {

    yr <- samples$year[i]; mo <- samples$month[i]
    t_now <- yr + (mo - 0.5) / 12          # current time (years)

    # cohorts recruited each year in month 'recruit_month'
    births <- seq_len(n_years + 1) - 1 + (recruit_month - 0.5) / 12

    for (b in births) {

      age <- t_now - b
      if (age <= 0.05 || age > 4) next

      n <- round(n0 * exp(-Z * age))       # exponential decay
      if (n < 5) next

      mu <- Linf * (1 - exp(-K * (age - t0)))
      if (mu <= 0) next

      len <- stats::rnorm(n, mean = mu, sd = cv * mu)
      len <- len[len > 8 & len < Linf * 1.05]   # rough selectivity
      if (!length(len)) next

      out[[length(out) + 1]] <- data.frame(
        Length = len, Year = yr, Month = mo
      )
    }
  }

  do.call(rbind, out)
}


test_that("the simulator produces usable data", {

  d <- simulate_lfq_data()

  expect_s3_class(d, "data.frame")
  expect_true(all(c("Length", "Year", "Month") %in% names(d)))
  expect_gt(nrow(d), 1000)
  expect_true(all(d$Length > 0))
  expect_lt(max(d$Length), 45 * 1.1)          # nothing above Linf
})


test_that("ELEFAN_GA recovers the simulated growth parameters", {

  skip_on_cran()                       # stochastic computation, ~1 min
  skip_if_not_installed("TropFishR")

  TRUE_LINF <- 45
  TRUE_K    <- 1.2

  d   <- simulate_lfq_data(Linf = TRUE_LINF, K = TRUE_K)
  lfq <- prepare_tropfish(d, bin_size = 1, species = "sim",
                          stock_name = "sim", length_unit = "cm")

  set.seed(123)
  fit <- run_elefan_ga(
    lfq,
    control = list(
      low_par      = list(Linf = 35, K = 0.5, t_anchor = 0),
      up_par       = list(Linf = 55, K = 2.0, t_anchor = 1),
      popSize      = 40,
      maxiter      = 30,
      seasonalised = FALSE
    )
  )

  p <- extract_growth_parameters(fit)

  expect_s3_class(fit, "FishStockELEFAN")
  expect_true(is.finite(p$Linf) && is.finite(p$K))

  # ELEFAN is noisy: we tolerate +/- 20% on Linf and +/- 50% on K
  # (K is notoriously less well identified than Linf).
  expect_equal(p$Linf, TRUE_LINF, tolerance = 0.20)
  expect_equal(p$K,    TRUE_K,    tolerance = 0.50)

  # phi' must remain in a plausible range for this (Linf, K) pair
  phi_true <- log10(TRUE_K) + 2 * log10(TRUE_LINF)
  expect_equal(p$phiL, phi_true, tolerance = 0.15)
})


test_that("estimate_M_all is consistent with theory", {

  skip_if_not_installed("TropFishR")

  m <- estimate_M_all(Linf = 45, K = 1.2, temp = 27)

  # all estimates are positive and finite
  expect_true(all(m$table$M > 0))
  expect_true(all(is.finite(m$table$M)))

  # the consensus is bounded by the min and max of the methods
  expect_gte(m$consensus[["geomean"]], min(m$table$M))
  expect_lte(m$consensus[["geomean"]], max(m$table$M))

  # M/K must remain in the biologically plausible range (~1 to 2.5)
  MK <- m$consensus[["geomean"]] / 1.2
  expect_gt(MK, 0.5)
  expect_lt(MK, 4)
})


test_that("Froese's reference points satisfy Lc_opt < Lopt < Linf", {

  for (MK in c(1.0, 1.5, 2.0)) {
    for (FM in c(0.5, 1, 2)) {

      rp <- lbb_reference_points(Linf = 45, MK = MK, FM = FM)

      expect_lt(rp$Lopt, 45)              # Lopt always below Linf
      expect_lt(rp$Lc_opt, rp$Lopt)       # capture before the optimal size
      expect_gt(rp$Lc_opt, 0)
    }
  }
})


test_that("reducing effort increases the SPR (expected monotonicity)", {

  skip_if_not_installed("LBSPR")

  pr <- lbspr_pars(Linf = 45, L50 = 22, L95 = 28, M = 1.5, K = 1.2)

  scen <- list(
    mse_measure("Status quo",   "statuquo"),
    mse_measure("-30%",        "effort", effort_reduction = 0.30),
    mse_measure("-60%",        "effort", effort_reduction = 0.60)
  )

  res <- run_mse_equilibrium(pr, scen, base_FM = 2)
  spr <- res$comparison$SPR

  expect_equal(length(spr), 3)
  expect_true(all(is.finite(spr)))

  # SPR must increase as effort decreases
  expect_true(spr[2] > spr[1])
  expect_true(spr[3] > spr[2])

  # SPR remains bounded in [0, 1]
  expect_true(all(spr >= 0 & spr <= 1))
})


test_that("the allometry recovers the simulated exponent b", {

  set.seed(7)
  L <- stats::runif(800, 10, 40)
  a_true <- 0.008; b_true <- 3.10
  W <- a_true * L^b_true * exp(stats::rnorm(length(L), 0, 0.05))

  lw <- fit_length_weight(data.frame(L = L, W = W), "L", "W")

  expect_equal(lw$b, b_true, tolerance = 0.03)
  expect_gt(lw$r2, 0.98)

  # b > 3 => positive allometry detected
  expect_match(lw$allometry, "positive")
})