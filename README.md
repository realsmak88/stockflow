# stockflow

<!-- badges: start -->
[![R-CMD-check](https://github.com/realsmak88/stockflow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/realsmak88/stockflow/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**stockflow** is a modular R package for **fish stock assessment based on
length-frequency data**, particularly suited to data-poor fisheries. It
provides robust wrappers around the field's reference tools (TropFishR, LBB,
LBSPR, MSEtool) within a coherent pipeline, without ever reimplementing their
algorithms.

## Pipeline

```
Import → Validation → LFQ preparation → Growth → Mortality
     → Allometry → Maturity → LBB / LBSPR → YPR / SPR
     → TAC + quota allocation
     → Biological rest period (closure) → Technical measures (Lc, mesh size)
     → HTML / PDF / DOCX report
```

The package covers the three objectives of a management plan: (i) biological
reference points → TAC → quotas per segment prorated by catches; (ii) optimal
windows and durations for biological rest periods; (iii) technical measures
(size at first capture, mesh size). A complete tutorial is provided as a
vignette (`vignette("tutoriel_stockflow", package = "stockflow")`) and eight
example (simulated) datasets are bundled (`data(package = "stockflow")`).

## Installation

From source (the package is not yet published):

```r
# install.packages("devtools")
devtools::install(".")     # from the repository root
# or, without installing, in development mode:
devtools::load_all(".")
```

Once the repository is published on GitHub:

```r
devtools::install_github("realsmak88/stockflow")
```

The package loads and works (import, validation, allometry) **without** the
heavy scientific dependencies. For the full pipeline (growth, LBB, LBSPR,
MSE), install the suggested packages:

```r
install.packages(c("TropFishR", "FSA", "LBSPR", "fishmethods", "MSEtool"))

# For LBB: Froese's official script is provided in inst/scripts/.
# It requires the R2jags package and the JAGS software (https://mcmc-jags.sourceforge.io).
install.packages("R2jags")
```

## Quick Start

```r
library(stockflow)

# 1. Importing a length file (automatic separator detection)
demo <- system.file("extdata", "penaeus_lengths_demo.csv", package = "stockflow")
lengths <- read_lengths(demo)

# 2. Quality validation
val <- validate_lengths(lengths, min_length = 0, max_length = 60)
print(val)

# 3. Length-weight allometry (no heavy dependency)
bio <- utils::read.csv(
  system.file("extdata", "penaeus_bio_demo.csv", package = "stockflow"))
lw <- fit_length_weight(bio, "length", "weight")
print(lw)   # a, b, r², isometry test

# 4. Preparing the LFQ object (requires TropFishR)
dd <- lengths
names(dd)[names(dd) == "length"] <- "Length"
names(dd)[names(dd) == "year"]   <- "Year"
names(dd)[names(dd) == "month"]  <- "Month"
lfq <- prepare_tropfish(dd, bin_size = 1, species = "Penaeus notialis")

# 5. Natural mortality: ~30 estimators + consensus
M <- estimate_M_all(Linf = 45, K = 1.2, temp = 27)
M$table       # one row per method (TropFishR, FSA::metaM, Jensen…)
M$consensus   # geometric mean / median
```

The complete tutorial, from the raw file to the assessment, is provided as a vignette:

```r
vignette("tutoriel_stockflow", package = "stockflow")
```

## Modules

| Module | File | Role |
|---|---|---|
| Import | `R/import.R` | Robust CSV/Excel import, harmonization, `FishStock` object |
| Validation | `R/validation.R` | Quality checks, overall score |
| LFQ preparation | `R/prepare_tropfish.R` | `lfq` object compatible with TropFishR |
| Growth | `R/growth.R` | Powell-Wetherall, ELEFAN (GA/SA), multi-criteria selection |
| Mortality | `R/mortality.R` | Z, M (multi-method), F, E |
| Allometry | `R/allometry.R` | Length-weight relationship |
| LBB | `R/lbb.R` | Length-Based Bayesian Biomass (Froese) |
| LBSPR | `R/lbspr.R` | Spawning Potential Ratio |
| MSE | `R/mse.R` | Simulation of management measures |

## Development

```r
devtools::load_all(".")   # loads the package in development mode
devtools::test()          # runs the test suite (testthat)
devtools::document()      # regenerates NAMESPACE + docs from roxygen tags
```

## License

GPL-3 — see [LICENSE.md](LICENSE.md).

## Author

Kamarel Ba (`bakamarel@gmail.com`).