test_that("read_lengths reads the demo dataset", {

  demo <- system.file("extdata", "penaeus_lengths_demo.csv",
                      package = "stockflow")
  testthat::skip_if(demo == "", "Demo file not found")

  lengths <- read_lengths(demo)

  expect_s3_class(lengths, "data.frame")
  expect_true("length" %in% names(lengths))
  expect_gt(nrow(lengths), 0)
})

test_that("read_catches skips cleanly if the data is missing", {

  f <- "data/raw/catches.csv"
  testthat::skip_if_not(file.exists(f), "File catches.csv missing (not provided)")

  catch <- read_catches(f)
  expect_true("catch" %in% names(catch))
})