# cran-comments.md — stockflow 0.1.0

## Summary

First submission to CRAN.

stockflow provides a reproducible stock assessment pipeline for fisheries
based on length-frequency data. The package implements no new scientific
algorithm: it composes, via documented and robust *wrappers*, the official
implementations from TropFishR, FSA, LBSPR and MSEtool.

## Test results

### Local

- macOS 26.3.1, aarch64-apple-darwin20, R 4.5.3
- `R CMD check --as-cran`: **0 error | 1 warning | 2 notes**
  - WARNING: `'qpdf' is needed for checks on size reduction of PDFs`
    — limitation of the local environment (qpdf not installed), irrelevant
    since the package contains no PDF.
  - NOTE: `unable to verify current time` — no NTP service available in
    the local verification environment.
  - NOTE: URL `https://github.com/realsmak88/stockflow` returning 404 — the
    public repository will be created upon publication; the URL will then
    be valid.
- `testthat`: 220 PASS, 0 FAIL, 2 SKIP
  (skipped tests: data file not distributed; ELEFAN_GA example
  skipped on CRAN)

### Other platforms

- [ ] win-builder (R-devel)
- [ ] R-hub (Linux, Windows, macOS)

*(To be completed before submission; these checks require a network
connection.)*

## Dependencies

- **Imports**: all available on CRAN.
- **Suggests**: `car`, `fishmethods`, `knitr`, `LBSPR`, `MSEtool`, `R2jags`,
  `rmarkdown`, `targets`, `testthat` — all on CRAN.
- The package **loads and works without the suggested dependencies**:
  each call is protected by `requireNamespace()` and raises an
  explicit message if the package is missing. The corresponding tests use
  `skip_if_not_installed()`.

## Points that may prompt a reviewer question

- **Embedded third-party script.** `inst/scripts/LBB_ggplot.R` is the official
  script for the LBB method (Froese, R., Winker, H., Coro, G. et al. 2018,
  *ICES Journal of Marine Science* 75(6): 2004-2015), redistributed with
  the agreement of its license, and executed without modification of its core scientific
  content. It requires the external software **JAGS**; all examples that call
  it are wrapped in `\dontrun{}`.

- **Execution time.** ELEFAN methods (genetic algorithm, simulated
  annealing) are inherently costly. Their examples are wrapped in `\dontrun{}`
  and the corresponding vignette chunks are not evaluated at build time,
  in order to stay within CRAN's time limits.

- **Interactivity.** `TropFishR::powell_wetherall()` and
  `TropFishR::catchCurve()` require an interactive selection of the regression
  interval when `reg_int` is not provided. The wrappers document this
  behavior and raise an actionable message; no example or test
  depends on interaction.

## Miscellaneous

- The package is at the *experimental* lifecycle stage (badge in the README).
- No function writes to the user's directory without an
  explicit path; `run_lbb_froese()` writes by default to `tempdir()`.