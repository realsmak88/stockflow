test_that("spawning_season computes the monthly proportion of mature individuals", {

  bio <- data.frame(
    Date = as.Date(c("2025-01-10", "2025-01-20",
                     "2025-06-05", "2025-06-15", "2025-06-25")),
    maturity = c("M", "IM", "M", "M", "IM")
  )
  sp <- spawning_season(bio, date_col = "Date", maturity_col = "maturity",
                        mature_codes = "M")
  expect_s3_class(sp, "data.frame")
  expect_equal(nrow(sp), 12)
  ## January: 1 mature / 2 -> 0.5
  expect_equal(sp$prop_mature[sp$month == 1], 0.5)
  ## June: 2 matures / 3
  expect_equal(sp$prop_mature[sp$month == 6], 2 / 3, tolerance = 1e-8)
  ## month with no data -> NA
  expect_true(is.na(sp$prop_mature[sp$month == 3]))

  ## missing maturity codes -> error
  expect_error(spawning_season(bio, "Date", "maturity"),
               "mature_codes")

})


test_that("recruitment_pattern returns 12 months summing to 100 (TropFishR)", {

  skip_if_not_installed("TropFishR")

  ## small synthetic lfq object with several dates
  set.seed(1)
  dates <- as.Date(paste0("2025-", sprintf("%02d", 1:12), "-15"))
  midL  <- seq(2, 40, 2)
  catch <- sapply(1:12, function(m) round(stats::rpois(length(midL), 20)))
  lfq <- list(midLengths = midL, catch = catch, dates = dates,
              Linf = 42, K = 0.45, t0 = 0)

  rc <- recruitment_pattern(lfq)
  expect_s3_class(rc, "data.frame")
  expect_equal(nrow(rc), 12)
  expect_equal(sum(rc$per_recruits), 100, tolerance = 1e-6)

})


test_that("optimal_closure finds the maximal circular window", {

  ## spawning signal peaking in December-January (tests the circular behavior)
  sp <- data.frame(month = 1:12, prop_mature = 0)
  sp$prop_mature[c(12, 1, 2)] <- c(0.9, 0.9, 0.9)

  cl <- optimal_closure(spawning = sp, duration = 3)
  expect_s3_class(cl, "FishStockClosure")
  ## the 3-month window must cover Dec-Jan-Feb (start = 12)
  expect_equal(cl$best$start_month, 12)
  expect_setequal(cl$best$months, c(12, 1, 2))
  ## almost all the score is protected
  expect_gt(cl$best$protected, 0.99)

})


test_that("optimal_closure compares durations and requires a signal", {

  rc <- data.frame(month = 1:12, per_recruits = 0)
  rc$per_recruits[6:8] <- c(30, 40, 30)

  cl <- optimal_closure(recruitment = rc, duration = NULL, max_duration = 5)
  expect_false(is.null(cl$by_duration))
  expect_equal(nrow(cl$by_duration), 5)
  ## the protected share grows with duration
  expect_true(all(diff(cl$by_duration$protected) >= -1e-9))
  ## 3-month window around the peak = Jun-Jul-Aug
  b3 <- cl$by_duration[cl$by_duration$duration == 3, ]
  expect_equal(b3$protected, 1, tolerance = 1e-8)

  ## no signal -> error
  expect_error(optimal_closure(), "At least one signal")
  ## duration out of bounds -> error
  expect_error(optimal_closure(recruitment = rc, duration = 12),
               "between 1 and 11")

})