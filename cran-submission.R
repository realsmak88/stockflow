###############################################################
#
# stockflow — Multi-platform CRAN submission script
#
# To be run from the ROOT of the package, on YOUR machine
# (network connection required; these checks run on
#  remote servers).
#
# Recommended order:
#   0. Prerequisites      (only once)
#   1. Local checks       (fast, everything must be green before sending)
#   2. URL check          (the GitHub URL must exist → push repo first)
#   3. win-builder        (Windows R-devel + R-release; result by email ~30 min)
#   4. R-hub v2            (Linux/Windows/macOS via GitHub Actions; repo required)
#   5. CRAN submission    (only when 1–4 are clean)
#
###############################################################

## ---- 0. Prerequisites (to install once) -----------------------------------
# install.packages(c("devtools", "rhub", "urlchecker",
#                     "spelling", "revdepcheck"))
# Windows toolchain not required locally: win-builder compiles server-side.

stopifnot(file.exists("DESCRIPTION"))
pkg <- read.dcf("DESCRIPTION", fields = "Package")[1, 1]
message("Package: ", pkg, " ", read.dcf("DESCRIPTION", fields = "Version")[1, 1])


## ---- 1. Local checks -----------------------------------------------------
# 1a. Full R CMD check --as-cran (vignette included, pandoc required)
#     -> MUST give 0 error. The only tolerated NOTE/WARNING are
#        environmental (qpdf, "unable to verify current time",
#        URL 404 as long as the GitHub repo is not yet public).
devtools::check(
  document = TRUE,          # regenerate man/ + NAMESPACE first
  args     = c("--as-cran"),
  build_args = character(),
  vignettes = TRUE
)

# 1b. Spelling (optional but recommended; set the language)
# spelling::spell_check_package()          # add false positives to
# spelling::update_wordlist()              #   inst/WORDLIST

# 1c. Test coverage (indicative target >= 70 %)
# covr::report(covr::package_coverage())


## ---- 2. URL check ----------------------------------------------------
# Requires that the public repository exist (otherwise the GitHub URL returns 404).
# Push first:  git push -u origin main
urlchecker::url_check()


## ---- 3. win-builder (Windows, CRAN server) ---------------------------------
# Results sent by email (Maintainer's address) within ~30 min.
devtools::check_win_devel()      # R-devel   (the strictest)
devtools::check_win_release()    # R-release


## ---- 4. R-hub v2 (multi-platform via GitHub Actions) ---------------------
# Prerequisite: public GitHub repository + `rhub::rhub_setup()` run once
# (adds the .github/workflows/rhub.yaml workflow and validates the account).
# rhub::rhub_setup()             # only once
# rhub::rhub_doctor()            # checks the configuration
rhub::rhub_check(
  platforms = c("linux", "windows", "macos", "macos-arm64")
)


## ---- 5. CRAN submission -----------------------------------------------------
# ONLY run when 1–4 are clean and cran-comments.md is up to date.
# devtools::submit_cran()

###############################################################
# After acceptance:
#   - git tag v0.1.0 && git push --tags
#   - create a GitHub "Release"
###############################################################

