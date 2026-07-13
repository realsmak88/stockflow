# Unit tests — report.R module
# rmarkdown rendering is skipped if rmarkdown/pandoc are absent;
# aggregation and the summary table are tested without dependency.

test_that("collect_results creates a FishStockResults object", {

  res <- collect_results(meta = list(species = "Penaeus notialis"))

  expect_s3_class(res, "FishStockResults")
  expect_null(res$growth)
  expect_equal(res$meta$species, "Penaeus notialis")

})


test_that("summary_table returns an empty table without results", {

  res <- collect_results()
  st  <- summary_table(res)

  expect_s3_class(st, "data.frame")
  expect_equal(nrow(st), 0)

})


test_that("summary_table extracts available indicators", {

  # minimal objects mimicking module outputs
  growth <- list(best = data.frame(method = "ELEFAN_GA", Linf = 45,
                                   K = 1.2, phiL = 3.3))

  allometry <- structure(
    list(a = 0.005, b = 2.9, r2 = 0.97, n = 500,
         b_ci = c(2.85, 2.95), allometry = "negative allometry (b < 3)"),
    class = "FishStockLW")

  mortality <- structure(
    list(Z = list(Z = 1.8),
         M_consensus = 1.2,
         summary = data.frame(Z = 1.8, M = 1.2, F = 0.6, E = 0.333,
                              FM = 0.5, status = "Underexploitation (E < 0.5)",
                              stringsAsFactors = FALSE)),
    class = "FishStockMortality")

  res <- collect_results(growth = growth, allometry = allometry,
                         mortality = mortality,
                         maturity = list(L50 = 24, L95 = 30))

  st <- summary_table(res)

  expect_true(nrow(st) > 5)
  expect_true("Linf" %in% st$Indicator)
  expect_true("b (length-weight)" %in% st$Indicator)
  expect_true("M (consensus)" %in% st$Indicator)
  expect_true("L50" %in% st$Indicator)
  expect_true(all(c("Indicator", "Value", "Source") %in% names(st)))

})


test_that("fishstock_report rejects an object of the wrong class", {

  expect_error(fishstock_report(list(a = 1)))

})


test_that("the template is properly provided with the package", {

  tpl <- system.file("rmarkdown", "rapport_stockflow.Rmd",
                     package = "stockflow")
  skip_if(tpl == "", "package not installed (load_all mode)")
  expect_true(file.exists(tpl))

})