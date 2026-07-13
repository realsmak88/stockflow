# stockflow — Publication Checklist (GitHub → CRAN)

Status as of July 13, 2026. Checked boxes are done; the others remain to
be handled.

## 1. Package structure

- [x] `R/` flat (subfolders are **not** supported by R; the former
      `R/prepare/` has been flattened)
- [x] No `library()` / `require()` in `R/` — everything in `package::function()`
- [x] No global variables (`<<-`)
- [x] Consistent S3 classes + `print` / `summary` / `plot` methods
- [x] Error handling via `tryCatch()`
- [x] `inst/extdata/` — demo datasets
- [x] `inst/scripts/LBB_ggplot.R` — Froese's official script (LBB backend)
- [x] `inst/rmarkdown/rapport_stockflow.Rmd` — report template
- [x] `.Rbuildignore` completed (`data/`, `output/`, `scripts/`, `_targets.R`,
      `LBB_Results/`, `renv/`, …)

## 2. DESCRIPTION (fixed)

- [x] **`tidyverse` removed from Imports** — meta-package, immediate CRAN
      rejection
- [x] **`LBB` removed from Suggests** — not on CRAN, would trigger
      "Packages suggested but not available"
- [x] `targets` moved to Suggests (pipeline tool, not a code dependency)
- [x] **`LazyData: true` (re)enabled**: 8 `data/*.rda` datasets
      (xz-compressed) provided per CRAN convention, documented in
      `R/data.R`, prepared by `data-raw/prepare_datasets.R`. The raw CSV
      analysis files have been moved to `data-raw/processed` and `data-raw/raw`
      (CRAN forbids subfolders in `data/`).
- [x] Imports: `FSA`, `TropFishR`, `ggplot2`, `graphics`, `janitor`, `methods`,
      `readr`, `readxl`, `stats`, `stringr`, `tibble`, `tools`, `utils`, `yaml`
- [x] Suggests: `car`, `fishmethods`, `knitr`, `LBSPR`, `MSEtool`, `R2jags`,
      `rmarkdown`, `targets`, `testthat`
- [x] `VignetteBuilder: knitr`, `Config/testthat/edition: 3`, `URL`, `BugReports`

## 3. Documentation (roxygen2)

- [x] `@export` added to the 27 remaining **public** functions
      (`read_maturity`, `read_catches`, `validate_stock`, `run_elefan_ga`,
      `compare_growth_models`, …). Internal helpers (`detect_separator`,
      `check_missing`, `get_Linf`, …) remain **unexported** — this is intentional.
- [x] `R/zzz-globals.R` — `utils::globalVariables()` to avoid NOTES
      "no visible binding for global variable" (aes/ggplot2)
- [x] **Complete roxygen documentation of the 27 public functions** (import: 10,
      validation: 4, growth: 13) — title, description, exhaustive `@param`,
      `@return`, `@examples` (`\dontrun{}` when a dependency or file is
      required), `@references`, `@seealso`.
- [x] S3 methods (`print.*`, `summary.*`, `plot.*`): `@export` alone, which is
      correct — roxygen emits `S3method()` and no `.Rd` file is required.
- [x] Aliases (`estimate_Z`, `calculate_F_E`) documented via `@rdname`.
- [ ] Re-run `devtools::document()` to regenerate `NAMESPACE` + `man/`

## 4. Tests

- [x] Standard structure: `tests/testthat.R` + `tests/testthat/test-*.R`
      (flat tests directly in `tests/` would **not** have been run by
      `R CMD check`)
- [x] 8 test files: import, mortality, allometry, LBB, LBB-Froese,
      LBSPR, MSE, report — heavy parts are behind `skip_if_not_installed()`
- [ ] Aim for ≥ 70% coverage (`covr::package_coverage()`)

## 5. Vignette

- [x] `vignettes/tutoriel_stockflow.Rmd` — complete pipeline, chunks protected
      by `requireNamespace()`
- [x] `devtools::load_all(".")` replaced by `library(stockflow)`
      (required: `load_all()` fails when building the vignette)
- [x] **Interactive chunks neutralized**: `run_powell()` calls
      `TropFishR::powell_wetherall()`, which requires an **interactive**
      selection of `reg_int` and returns `NULL` under `knitr` → was causing
      `R CMD build` to fail. The chunk is now set to `eval=FALSE` with an
      explicit `reg_int`.
- [x] **Long chunks neutralized**: ELEFAN_GA/SA (1–3 min) set to `eval=FALSE`;
      typical Linf/K values are injected so the rest runs
      (CRAN strongly limits vignette build time).
- [x] `check_powell()` error message made actionable (mentions
      `reg_int` and the non-interactive context).

## 6. Reproducible pipeline

- [x] `_targets.R` — complete DAG: import → validation → LFQ → growth →
      mortality → allometry → maturity → LBB → LBSPR → MSE → report
- [ ] `renv::init()` then `renv::snapshot()` to lock package versions

## 7. `R CMD check` status

**R CMD check `--as-cran`** (R 4.5.3, aarch64-apple-darwin20, vignette compiled
with pandoc 3.10): **0 error**. The only remaining messages are environmental
or expected:

- WARNING `'qpdf' is needed…` — qpdf absent from the local environment (the
  package contains no PDF).
- NOTE `unable to verify current time` — no NTP service in the environment.
- NOTE URL `github.com/realsmak88/stockflow` (404) — will disappear once the
  public repository is pushed.

Without the incoming checks (`_R_CHECK_CRAN_INCOMING_=FALSE`), the check is
**Status: OK, 0/0/0** (50 checks OK).

Tests: **220 PASS, 0 FAIL, 2 SKIP** (data file not distributed;
ELEFAN_GA example skipped on CRAN).

Fixes that made this possible:

- **non-ASCII**: accented characters escaped as `\uXXXX` in strings (display
  remains accented) and transliterated in comments/roxygen.
- **undeclared `requireNamespace("LBB")`**: since LBB is not on CRAN, the
  "package LBB" backend (`run_lbb()`, `.check_lbb()`, `.get_lbb_fun()`) was
  removed. Only `run_lbb_froese()` remains (official script + `R2jags`).
- **`\usage` / arguments**: a pre-existing roxygen block was attached to the
  wrong function (`create_default_config` carried `import_data`'s `@param`s);
  missing `@param`s added (`t0`, `Lmean`, `use_fsa`, `species_config`, `...`).
- **`quantile`** → `stats::quantile()`.
- **roxygen**: `\code{pars@L50}` escaped as `pars@@L50` (`@` is a roxygen tag).

## 8. Remaining before CRAN submission

- [x] `README.md`: badges, installation, minimal example, mention of JAGS for LBB
- [x] `cran-comments.md` (up-to-date check results, tested platforms)
- [x] Multi-platform submission script: `cran-submission.R`
- [ ] Push the public GitHub repository (`git push -u origin main`) → clears
      the 404 NOTE on the URL
- [ ] Remote multi-platform checks (see `cran-submission.R`):

```r
urlchecker::url_check()        # after the push (URL must exist)
devtools::check_win_devel()    # Windows R-devel (result by email)
devtools::check_win_release()  # Windows R-release
rhub::rhub_check()             # Linux / Windows / macOS
devtools::submit_cran()        # only once everything is green
```

- [ ] (optional) `covr::package_coverage()` — aim for ≥ 70%
- [ ] (optional) `renv::snapshot()` — lock package versions

## 9. Known points of attention

- **JAGS**: `run_lbb_froese()` requires `R2jags` + the JAGS software. All
  corresponding examples are in `\dontrun{}`.
- **Execution time**: ELEFAN_GA/SA is slow. Examples must stay in
  `\dontrun{}` or use minimal settings (CRAN limits to ~5 s/example).
- **`LBB`** off CRAN: document the GitHub installation in the README, or
  rely solely on the provided Froese script backend.
- Encoding: `Encoding: UTF-8` declared; avoid accented characters in the
  **code** (accented comments/docs are tolerated).