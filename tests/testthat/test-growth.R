## Growth protocol tests: Powell -> K-Scan -> ELEFAN
## (surface de reponse) -> courbe VBGF a travers les modes.

# Builds a small lfq object from the embedded data.
.make_lfq <- function() {
  skip_if_not_installed("TropFishR")
  freq <- subset(cymbium_freq, species == "Cymbium cymbium")
  freq$Length <- freq$LCQ
  freq$Year   <- as.integer(format(freq$Date, "%Y"))
  freq$Month  <- as.integer(format(freq$Date, "%m"))
  prepare_tropfish(freq, bin_size = 2, species = "Cymbium cymbium",
                   length_unit = "cm")
}

test_that("run_kscan estimates K with fixed Linf (numeric or Powell object)", {
  lfq <- .make_lfq()

  ks <- run_kscan(lfq, Linf = 42, K_range = seq(0.2, 1.0, 0.2))
  expect_s3_class(ks, "FishStockKScan")
  expect_equal(ks$Linf, 42)
  expect_true(is.finite(ks$K) && ks$K > 0)
  ## the score profile covers the scanned range
  expect_true(all(c("K", "score") %in% names(ks$scan)))
  expect_equal(nrow(ks$scan), 5L)
  ## the optimal K corresponds to the maximum score of the profile
  expect_equal(ks$K, ks$scan$K[which.max(ks$scan$score)], tolerance = 1e-6)
})

test_that("run_kscan accepts a FishStockPowell object", {
  lfq <- .make_lfq()
  pw  <- run_powell(lfq, reg_int = c(2, 12))
  ks  <- run_kscan(lfq, Linf = pw, K_range = seq(0.2, 1.0, 0.2))
  expect_s3_class(ks, "FishStockKScan")
  expect_equal(ks$Linf, pw$result$Linf_est, tolerance = 1e-6)
})

test_that("run_kscan rejects an invalid Linf", {
  lfq <- .make_lfq()
  expect_error(run_kscan(lfq, Linf = "abc"))
  expect_error(run_kscan(lfq, Linf = 42, K_range = 0.5))  # a single value
})

test_that("plot_response_surface produces a ggplot (2D and 1D)", {
  lfq <- .make_lfq()

  el <- run_elefan(lfq, Linf_range = seq(36, 44, 4),
                   K_range = seq(0.2, 1.0, 0.4), MA = 5)
  p2 <- plot_response_surface(el)
  expect_s3_class(p2, "ggplot")

  ks <- run_kscan(lfq, Linf = 42, K_range = seq(0.2, 1.0, 0.2))
  p1 <- plot_response_surface(ks)
  expect_s3_class(p1, "ggplot")
})

test_that("plot_growth_modes works on ELEFAN, K-Scan and list", {
  lfq <- .make_lfq()

  ## from a list of parameters
  pm_list <- plot_growth_modes(lfq, list(Linf = 42, K = 0.45), MA = 5)
  expect_s3_class(pm_list, "ggplot")

  ## from a K-Scan (flat fields)
  ks <- run_kscan(lfq, Linf = 42, K_range = seq(0.2, 1.0, 0.2))
  pm_ks <- plot_growth_modes(lfq, ks, MA = 5)
  expect_s3_class(pm_ks, "ggplot")

  ## from an ELEFAN (Linf/K under $result$par)
  el <- run_elefan(lfq, Linf_range = seq(36, 44, 4),
                   K_range = seq(0.2, 1.0, 0.4), MA = 5)
  pm_el <- plot_growth_modes(lfq, el, MA = 5)
  expect_s3_class(pm_el, "ggplot")
})

test_that("plot_growth_modes uses the best model from a full analysis", {
  lfq <- .make_lfq()

  ## Full analysis (ELEFAN grid only, fast): we check that
  ## .extract_growth_par correctly reads Linf/K from the best retained model ($best),
  ## and not systematically the first fitted model.
  an <- run_growth_analysis(lfq, methods = "elefan",
                            Linf_range = c(36, 46), K_range = c(0.2, 1.0),
                            powell = FALSE, verbose = FALSE)
  expect_s3_class(an, "FishStockGrowthAnalysis")

  gp <- stockflow:::.extract_growth_par(an)
  expect_equal(gp$Linf, an$best$Linf, tolerance = 1e-6)
  expect_equal(gp$K,    an$best$K,    tolerance = 1e-6)

  pm <- plot_growth_modes(lfq, an, MA = 5)
  expect_s3_class(pm, "ggplot")
})

test_that("the plot S3 methods return the object invisibly", {
  lfq <- .make_lfq()
  ks  <- run_kscan(lfq, Linf = 42, K_range = seq(0.2, 1.0, 0.2))
  expect_invisible(plot(ks))
})

test_that("estimate_t0 matches Pauly's empirical formula", {
  ## Pauly (1979): log10(-t0) = -0.3922 - 0.2752 log10(Linf) - 1.0380 log10(K)
  expected <- -(10 ^ (-0.3922 - 0.2752 * log10(42) - 1.0380 * log10(0.45)))
  expect_equal(estimate_t0(Linf = 42, K = 0.45), expected, tolerance = 1e-9)
  ## t0 is expected to be (slightly) negative
  expect_lt(estimate_t0(Linf = 30, K = 0.3), 0)
  ## input guards
  expect_error(estimate_t0(Linf = -1, K = 0.4), "Linf")
  expect_error(estimate_t0(Linf = 42, K = 0),   "K")
  expect_error(estimate_t0(Linf = c(40, 42), K = 0.4), "Linf")
})
