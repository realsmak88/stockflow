# Unit tests — Froese backend script (lbb.R)
# Do NOT run JAGS: only check file preparation
# and the non-destructive patching of the script.

test_that("write_lbb_id writes an ID file with the correct columns", {

  tmp <- tempfile(fileext = ".csv")
  id  <- write_lbb_id(stock = "Penaeus", lf_file = "Penaeus_LF.csv",
                      species = "Penaeus notialis", mm = TRUE,
                      MK_prior = 1.5, Linf_prior = 45, file = tmp)

  expect_true(file.exists(tmp))
  need <- c("Stock", "Species", "File", "mm.user", "GausSel",
            "MergeLF", "Pile", "MK.user", "Linf.user", "Lm50")
  expect_true(all(need %in% names(id)))
  expect_equal(id$MK.user, 1.5)
  expect_true(id$mm.user)

})


test_that("write_lbb_inputs writes LF + ID and returns the paths", {

  wd <- file.path(tempdir(), paste0("lbbtest_", as.integer(runif(1, 1, 1e6))))

  dd <- data.frame(
    Stock   = "Penaeus",
    Year    = rep(c(2025, 2026), each = 3),
    Length  = c(18, 20, 22, 18, 20, 22),
    CatchNo = c(10, 25, 12, 8, 30, 15)
  )

  inp <- write_lbb_inputs(dd, stock = "Penaeus", workdir = wd,
                          mm = TRUE, MK_prior = 1.5)

  expect_true(file.exists(file.path(wd, inp$lf_file)))
  expect_true(file.exists(file.path(wd, inp$id_file)))

  lf <- utils::read.csv(file.path(wd, inp$lf_file))
  expect_true(all(c("Year", "Length", "CatchNo") %in% names(lf)))

  unlink(wd, recursive = TRUE)

})


test_that(".patch_lbb_script neutralizes rm() and injects Stock/ID.File", {

  # mini script mimicking the header of the Froese script
  src <- tempfile(fileext = ".R")
  writeLines(c(
    "rm(list=ls(all=TRUE))",
    'Stock       <-  "Sardinella aurita"',
    'ID.File     <-  "Aurita_ID.csv"',
    'output_dir <- "LBB_Results"',
    "x <- 1"
  ), src)

  patched <- .patch_lbb_script(src, "Penaeus",
                               "Penaeus_ID.csv", "LBB_Results")
  txt <- readLines(patched)

  expect_true(any(grepl("disabled by stockflow", txt)))       # rm() neutralized
  expect_true(any(grepl('Stock <- "Penaeus"', txt, fixed = TRUE)))
  expect_true(any(grepl('ID.File <- "Penaeus_ID.csv"', txt, fixed = TRUE)))
  expect_false(any(grepl("^rm\\(list", txt)))

})