# Tutorial — running the complete stockflow workflow

Usage guide for `scripts/full_workflow.R`, which evaluates the
**6 stocks** of an example EEZ from end to end and produces an HTML report
per stock plus a comparative summary.

---

## 1. What the script does

For each of the 6 stocks — octopus, white shrimp, and the 4 volutes
(*Cymbium cymbium*, *glans*, *marmoratum*, *pepo*) — it runs through:

```
loading → validation → LFQ object → growth (ELEFAN GA + SA, auto selection)
   → natural mortality M (~30 methods + consensus)
   → total mortality Z, then F and E
   → length-weight allometry → maturity ogive (L50, L95)
   → LBB (if JAGS) → LBSPR (SPR) → management scenarios → HTML report
```

A stock that fails does not interrupt the others: each step is isolated and
the failure is logged.

---

## 2. Prerequisites

```r
install.packages(c("TropFishR", "FSA", "LBSPR", "fishmethods",
                   "ggplot2", "rmarkdown", "knitr"))

# Optional — for the LBB step:
install.packages("R2jags")   # also requires the JAGS software
```

The cleaned data must be in `data/processed/`
(outputs of `R/prepare/01_clean_data.R`):

```
octopus_frequence_clean.csv   octopus_biologie_clean.csv
penaeus_frequence_clean.csv   penaeus_biologie_clean.csv
cymbium_frequence_clean.csv   cymbium_biologie_clean.csv
```

---

## 3. Running it

**The package does not need to be installed, nor published on GitHub.** Just
open the `stockflow.Rproj` project (the working directory is then set to the
package root) and load it from source:

```r
devtools::load_all(".")                      # loads stockflow from source
source("scripts/full_workflow.R")
```

The script takes care of `load_all()` itself if it detects that the package is
not installed — a simple `source()` is therefore enough.

If you prefer to install it once and for all:

```r
devtools::install(".")   # then: library(stockflow)
```

> ⚠️ `devtools::install_github("realsmak88/stockflow")` **cannot work as
> long as the repository has not been pushed to GitHub** (404 error).
> See §9 for publishing.

**Check the working directory** — the script reads `data/processed/` and
writes to `output/workflow/`, using relative paths:

```r
getwd()   # must end with /stockflow
```

Expect **15 to 40 minutes** depending on the machine (ELEFAN is the costly
step). For a quick first trial, set in BLOCK 1:

```r
ELEFAN_POP  <- 30
ELEFAN_ITER <- 20
```

---

## 4. The only parameter that requires your judgment: `reg_int`

This is the critical point of the workflow, and it needs to be looked at.

**What it is.** `reg_int` gives the two bounds of the **descending and
fully recruited** part of the linearized catch curve. It is on this
portion that the slope gives **Z**. TropFishR normally asks you to click these
bounds on the graph; the script fixes them so it can run without interaction.

**Why it matters.** Z determines F, then E = F/Z, hence **the stock's
status**. A poorly placed `reg_int` can turn a stock from "healthy" into
"overexploited." This is not a cosmetic setting.

**Procedure (normal back-and-forth, to be done once per stock)**

1. Run the script as is.
2. Open `output/workflow/<stock>_diagnostic_catchcurve.png`.
3. Identify the **peak** of the curve (full recruitment), then the
   **descending linear** portion right after. Discard the right-hand tail,
   often noisy due to low numbers.
4. Enter the indices in BLOCK 1:

```r
REG_INT <- list(
  "Octopus vulgaris"   = c(8, 18),   # ← adjust based on the graph
  "Penaeus notialis"   = c(15, 30),
  ...
)
```

5. Re-run.

The values provided are **plausible starting points**, not truths.

---

## 5. Reading the outputs

Everything ends up in `output/workflow/`:

| File | Content |
|---|---|
| `rapport_<stock>.html` | Full report: growth, M by method, Z/F/E, allometry, SPR, management scenarios, recommendation |
| `synthese_stocks.csv` | Comparative table of the 6 stocks (Linf, K, M, Z, F, E, L50, SPR, status) |
| `synthese_exploitation.png` | Exploitation rate E per stock, with the E = 0.5 reference |
| `<stock>_LFQ.png` | Restructured size frequencies (cohorts) |
| `<stock>_diagnostic_catchcurve.png` | **To be checked to set `reg_int`** |
| `lbb/<stock>/LBB_Results/` | LBB outputs (if JAGS): ggplot2 figures + annual table |

**Interpretation guidelines**

- `E = F/Z`: **E > 0.5** suggests overexploitation (Gulland optimum).
- `SPR`: **≥ 0.40** healthy · **0.20–0.40** concerning · **< 0.20** overexploited.
- Scenarios: the report automatically designates the one that **maximizes
  yield while keeping SPR above the target**. If none achieves this, it
  raises an alert — that is a result in itself.

---

## 6. Two alerts the script will send you

**`!! Linf or K hits a bound`**
ELEFAN converged on the edge of the search space: the estimate is
not reliable, it is constrained by your bounds rather than by the data.
Widen `Linf_min/Linf_max` or `K_min/K_max` for that stock in BLOCK 1.

**`!! poorly constrained ogive`**
Less than 15% (or more than 85%) mature in the sample: the logistic
does not have enough contrast, **L50 is unreliable**. This is the known case
of shrimp (≈ 7% mature). Two options: redo the ogive on **females
only**, or use a literature value (*P. notialis*: L50 ≈ 20–24 mm
cephalothorax length). Do not take the raw L50 without checking it.

---

## 7. Limitations to keep in mind

These are limitations of the **data**, not of the code — better to state
them than to discover them in a management committee.

- **Sampling duration (~1 year).** ELEFAN tracks cohort progression
  *over time*. With only one year, **K is fragile** (Linf holds up better).
  LBB, which works by year, has only 2 points. The results are
  **orders of magnitude**, to be consolidated with 2–3 years.
- **Volutes.** The 4 species are evaluated **separately**: their Linf range
  from about 15 cm (*marmoratum*) to 38 cm (*glans*). Grouping them into a
  single LFQ would create artificial modes that ELEFAN would read as growth.
- **Cephalopods and gastropods.** ELEFAN and the VBGF were designed for
  fish. Octopus (non-asymptotic growth, lifespan ~1–2 years) and
  volutes are poorly suited to it. Estimates must be cross-checked with
  regional literature.
- **Shrimp MLS.** Do not derive the minimum size from the raw L50 (see §6).

---

## 8. Moving on to management

The workflow stops at the per-stock diagnostic and equilibrium scenarios. For
**TACs, quotas by segment/fleet, and spatio-temporal closures**,
use the management scripts:

```r
source("R/analysis/04_stock_assessment.R")   # reference points (CMSY/JABBA)
source("R/analysis/05_management.R")         # TAC, quotas, closures, MLS
```

And for the closed-loop simulation under uncertainty (testing a harvest
control rule *before* adopting it):

```r
source("R/analysis/06_mse_openMSE.R")
```

---

## 9. Publishing the package on GitHub (optional)

`install_github()` returns a **404 error** as long as the repository does not
exist. Nothing is broken: the package is complete locally. To publish it:

```bash
# 1. Create the repository on github.com (account realsmak88), name: stockflow
#    → DO NOT initialize it with a README (yours already exists)

# 2. From the package root:
git init
git add .
git commit -m "stockflow 0.1.0 — complete pipeline, clean R CMD check"
git branch -M main
git remote add origin https://github.com/realsmak88/stockflow.git
git push -u origin main
```

Warning before pushing: `data/` contains your **raw survey data**. If it
should not be public, add it to `.gitignore`:

```
data/
output/
LBB_Results/
```

Only then will `devtools::install_github("realsmak88/stockflow")`
work.