# =============================================================================
# stockflow :: 04_stock_assessment.R
# Phase 2D — Data-limited stock assessment & biological reference points
#   Main method: CMSY++ / JABBA (Schaefer/Pella surplus production)
#   Inputs: total annual catches + standardized CPUE (03).
#   Outputs: r, k, MSY, Bmsy, Fmsy, B/Bmsy, F/Fmsy, Kobe.
#
# NB. Two equivalent pathways are provided: (A) CMSY (datalimited2), robust and
# fast; (B) JABBA (Bayesian, JAGS) for full uncertainty. SPiCT optional.
# =============================================================================
suppressPackageStartupMessages({library(tidyverse)})

proc <- "data/processed"; tab <- "output/tables"; fig <- "output/figures"

catch <- read_csv(file.path(proc,"total_catches_annual.csv"), show_col_types = FALSE)

# resilience -> prior on r (Froese et al. 2017, CMSY)
resil <- tibble::tribble(
  ~espece,              ~r_lo, ~r_hi, ~resilience,
  "Octopus vulgaris",   0.60,  1.50,  "High",
  "Penaeus notialis",   0.60,  1.50,  "High",
  "Cymbium spp",        0.20,  0.80,  "Medium")

# ---------------------------------------------------------------------------
# PATHWAY A — CMSY (datalimited2)  # remotes::install_github("datalimited/datalimited")
# ---------------------------------------------------------------------------
run_cmsy <- function(sp) {
  g  <- catch |> filter(espece == sp) |> arrange(annee)
  pr <- resil |> filter(espece == sp)
  cp <- file.path(tab, sprintf("cpue_std_%s.csv", tolower(word(sp,1))))
  idx <- if (file.exists(cp)) read_csv(cp, show_col_types = FALSE) else NULL

  fit <- datalimited2::cmsy2(
    year = g$annee, catch = g$capture_t,
    resilience = pr$resilience,
    # if CPUE available -> BSM (Bayesian Schaefer) via bsm() below
    verbose = FALSE)
  # ref. points
  rp <- fit$ref_pts
  tibble(stock = sp,
         r = rp["r","est"], k = rp["k","est"], MSY = rp["msy","est"],
         Bmsy = rp["k","est"]/2, Fmsy = rp["r","est"]/2,
         B_Bmsy = tail(fit$ref_ts$bbmsy, 1),
         F_Fmsy = tail(fit$ref_ts$ffmsy, 1))
}

# ---------------------------------------------------------------------------
# PATHWAY B — JABBA (Bayesian, recommended when CPUE is reliable)   library(JABBA)
# ---------------------------------------------------------------------------
run_jabba <- function(sp) {
  g   <- catch |> filter(espece == sp) |> arrange(annee)
  cdf <- data.frame(Year = g$annee, Total = g$capture_t)
  cp  <- file.path(tab, sprintf("cpue_std_%s.csv", tolower(word(sp,1))))
  idx <- read_csv(cp, show_col_types = FALSE)
  idf <- data.frame(Year = g$annee) |>
    dplyr::left_join(dplyr::transmute(idx, Year = annee, CPUE = index_scaled), by = "Year")
  pr  <- resil |> filter(espece == sp)

  jbinput <- JABBA::build_jabba(
    catch = cdf, cpue = idf[,c("Year","CPUE")],
    model.type = "Schaefer",
    r.prior = c(mean(c(pr$r_lo, pr$r_hi)), 0.4),   # lognormal (mean, cv)
    K.prior = c(20*max(cdf$Total), 1),
    psi.prior = c(0.9, 0.1),                        # B1/K ~ lightly exploited at start
    igamma = c(4, 0.01))
  fit <- JABBA::fit_jabba(jbinput, quickmcmc = TRUE, ni = 30000, nb = 5000, nc = 3)
  JABBA::jabba_plots(fit, output.dir = fig)
  fit
}

# ---- run CMSY for the 3 stocks --------------------------------------
brp <- purrr::map_dfr(c("Octopus vulgaris","Penaeus notialis","Cymbium spp"),
                      purrr::possibly(run_cmsy, otherwise = NULL))
write_csv(brp, file.path(tab, "stock_assessment_BRP.csv"))
print(brp)

# ---- Kobe (ggplot) ----------------------------------------------------------
kobe <- ggplot(brp, aes(B_Bmsy, F_Fmsy)) +
  annotate("rect", xmin=0, xmax=1,  ymin=1, ymax=Inf, fill="#FF7F7F", alpha=.35) +
  annotate("rect", xmin=0, xmax=1,  ymin=0, ymax=1,   fill="#FFF59D", alpha=.5) +
  annotate("rect", xmin=1, xmax=Inf,ymin=1, ymax=Inf, fill="#FFF59D", alpha=.5) +
  annotate("rect", xmin=1, xmax=Inf,ymin=0, ymax=1,   fill="#90EE90", alpha=.4) +
  geom_hline(yintercept=1, lty=2) + geom_vline(xintercept=1, lty=2) +
  geom_point(size=4) + ggrepel::geom_text_repel(aes(label = word(stock,1))) +
  labs(x=expression(B/B[MSY]), y=expression(F/F[MSY]),
       title="Kobe plot — stock status") + theme_bw()
ggsave(file.path(fig,"fig8_kobe.png"), kobe, width=7, height=6, dpi=130)

message("\n=== Phase 2D completed: BRP in output/tables/stock_assessment_BRP.csv ===")
# For full uncertainty: uncomment -> fit <- run_jabba("Octopus vulgaris")

