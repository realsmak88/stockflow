###############################################################################
#
#  COMPLETE STOCK ASSESSMENT WORKFLOW
#  Octopus (Octopus vulgaris) · White shrimp (Penaeus notialis) ·
#  Volutes (Cymbium cymbium, glans, marmoratum, pepo)
#
#  Runs the entire chain, without interaction:
#    import → validation → LFQ → growth (ELEFAN) → mortality (Z, M, F, E)
#    → allometry → maturity → LBB → LBSPR → management measures → HTML report
#
#  USAGE:  from the root of the stockflow package
#             source("scripts/full_workflow.R")
#
#  All settings are in BLOCK 1 below. The only parameter that
#  requires your judgment is `reg_int` (see the note at the top of the block).
#
###############################################################################

# --- Loading stockflow ------------------------------------------------
# The package does not need to be installed: if it isn't, it is loaded
# directly from source (development mode).
if (requireNamespace("stockflow", quietly = TRUE)) {
  suppressPackageStartupMessages(library(stockflow))
  message("stockflow loaded (installed package).")
} else if (requireNamespace("pkgload", quietly = TRUE) &&
           file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
  message("stockflow loaded from source (load_all).")
} else {
  stop(
    "stockflow not found.\n",
    "  - From the package root:  devtools::load_all(\".\")\n",
    "  - Or to install it once:  devtools::install(\".\")",
    call. = FALSE
  )
}

t_start <- Sys.time()


###############################################################################
#  BLOCK 1 — CONFIGURATION
###############################################################################

DATA_DIR <- "data/processed"     # outputs from cleaning (01_clean_data.R)
OUT_DIR  <- "output/workflow"    # reports, figures, tables
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TEMP_EAU   <- 25      # deg C, example average temperature (Pauly method)
SPR_CIBLE  <- 0.40    # management target
SPR_LIMITE <- 0.20    # management limit
BASE_FM    <- 1.5     # reference F/M (presumed current state)

# --- ELEFAN: computational cost -------------------------------------------------
# "Production" settings. For a quick first try, switch to
# ELEFAN_POP = 30 / ELEFAN_ITER = 20 (noisier results).
ELEFAN_POP  <- 60
ELEFAN_ITER <- 50

set.seed(42)   # ELEFAN is stochastic: the seed makes the run reproducible


# --- reg_int: regression interval of the catch curve --------------
#
# `reg_int` = the two bounds (size class indices) of the DESCENDING and
# FULLY RECRUITED part of the linearized catch curve.
# Z is derived from its slope, hence F, E, and the stock status.
#
# NULL  ->  AUTOMATIC SELECTION (auto_reg_int): the best interval is
#           chosen among the candidates, rejecting those that give a Z <= M
#           (negative F, biologically impossible).
#
# You can always impose your own bounds after looking at
# `output/workflow/<stock>_diagnostic_catchcurve.png`, e.g.:
#     "Octopus vulgaris" = c(8, 18)
#
# ⚠️ Automation reproduces a visual judgment, it does not replace it:
#    check the curves before releasing a Z.
#
REG_INT <- list(
  "Octopus vulgaris"   = NULL,
  "Penaeus notialis"   = NULL,
  "Cymbium cymbium"    = NULL,
  "Cymbium glans"      = NULL,
  "Cymbium marmoratum" = NULL,
  "Cymbium pepo"       = NULL
)


# --- Stock definitions ---------------------------------------------------
# Each row = one stock assessed separately.
#
# SCIENTIFIC NOTE: the 4 Cymbium are assessed SEPARATELY. Grouping them
# into a single LFQ would be a mistake: they have very different Linf
# (marmoratum ~15 cm vs glans ~38 cm) and mixing them creates artificial
# cohort modes that ELEFAN would interpret as growth.
#
STOCKS <- data.frame(
  stock      = c("Octopus vulgaris", "Penaeus notialis",
                 "Cymbium cymbium", "Cymbium glans",
                 "Cymbium marmoratum", "Cymbium pepo"),
  fichier    = c("octopus", "penaeus",
                 "cymbium", "cymbium", "cymbium", "cymbium"),
  col_long   = c("LM", "LCT", "LCQ", "LCQ", "LCQ", "LCQ"),
  unite      = c("cm", "mm", "cm", "cm", "cm", "cm"),
  bin        = c(1,    1,    1,    1,    0.5,  1),
  # ELEFAN search bounds (min/max Linf, min/max K)
  Linf_min   = c(20,   38,   16,   32,   12,   24),
  Linf_max   = c(32,   55,   24,   45,   20,   34),
  K_min      = c(0.4,  0.5,  0.15, 0.10, 0.20, 0.15),
  K_max      = c(1.6,  2.0,  0.80, 0.60, 0.90, 0.70),
  mat_type   = c("num", "num", "txt", "txt", "txt", "txt"),
  stringsAsFactors = FALSE
)


###############################################################################
#  BLOCK 2 — TOOLS
###############################################################################

log_etape <- function(...) cat(sprintf("   [%s] %s\n",
                                       format(Sys.time(), "%H:%M:%S"),
                                       paste0(...)))

# Runs a step while isolating errors: a stock that fails must not
# bring down the other five.
essayer <- function(libelle, expr) {
  res <- tryCatch(expr, error = function(e) {
    log_etape("FAILURE ", libelle, " : ", conditionMessage(e))
    NULL
  })
  if (!is.null(res)) log_etape("ok   ", libelle)
  res
}

# Loads a stock's data and puts it in the expected format
charger_stock <- function(cfg) {

  freq <- utils::read.csv(file.path(DATA_DIR,
                paste0(cfg$fichier, "_frequence_clean.csv")))
  bio  <- utils::read.csv(file.path(DATA_DIR,
                paste0(cfg$fichier, "_biologie_clean.csv")))

  # filter the species (useful for the 4 Cymbium sharing a file)
  if ("species" %in% names(freq)) freq <- freq[freq$"species" == cfg$stock, ]
  if ("species" %in% names(bio))  bio  <- bio[bio$"species"  == cfg$stock, ]

  # lengths: Length / Year / Month (prepare_tropfish format)
  d <- data.frame(
    Length = suppressWarnings(as.numeric(freq[[cfg$col_long]])),
    Year   = as.integer(format(as.Date(freq$Date), "%Y")),
    Month  = as.integer(format(as.Date(freq$Date), "%m"))
  )
  d <- d[is.finite(d$Length) & !is.na(d$Year) & !is.na(d$Month), ]

  # biology: length, weight, maturity
  b <- data.frame(
    length   = suppressWarnings(as.numeric(bio[[cfg$col_long]])),
    weight   = suppressWarnings(as.numeric(bio$weight)),
    maturite = as.character(bio$"maturity")
  )
  b <- b[is.finite(b$length), ]

  # mature: stages >= 3 (octopus/shrimp); M or P (volutes)
  b$mature <- if (cfg$mat_type == "num") {
    as.integer(suppressWarnings(as.numeric(b$maturite)) >= 3)
  } else {
    as.integer(toupper(trimws(b$maturite)) %in% c("M", "P"))
  }

  list(lfq_data = d, bio = b)
}

# Logistic maturity ogive -> L50, L95
ogive_maturite <- function(b) {
  d <- b[is.finite(b$length) & !is.na(b$mature), ]
  if (length(unique(d$mature)) < 2) return(NULL)   # no contrast
  g  <- stats::glm(mature ~ length, data = d, family = stats::binomial())
  cf <- stats::coef(g)
  list(L50 = as.numeric(-cf[1] / cf[2]),
       L95 = as.numeric((log(0.95 / 0.05) - cf[1]) / cf[2]),
       n   = nrow(d),
       prop_mature = mean(d$mature))
}


###############################################################################
#  BLOCK 3 — MAIN LOOP OVER STOCKS
###############################################################################

resultats  <- list()
synthese   <- list()

for (i in seq_len(nrow(STOCKS))) {

  cfg <- STOCKS[i, ]
  sp  <- cfg$stock

  cat("\n", strrep("=", 74), "\n", sep = "")
  cat("  STOCK ", i, "/", nrow(STOCKS), " : ", sp, "\n", sep = "")
  cat(strrep("=", 74), "\n", sep = "")

  # --- 3.1 Loading ------------------------------------------------------
  dat <- essayer("data loading", charger_stock(cfg))
  if (is.null(dat) || nrow(dat$lfq_data) < 100) {
    log_etape("stock skipped (insufficient sample size)"); next
  }
  log_etape(nrow(dat$lfq_data), " lengths · ", nrow(dat$bio), " biology individuals")

  # --- 3.2 Validation ------------------------------------------------------
  val <- essayer("validation", {
    v <- validate_lengths(
      data.frame(length = dat$lfq_data$Length),
      max_length = cfg$Linf_max * 1.2
    )
    v
  })

  # --- 3.3 LFQ object -------------------------------------------------------
  lfq <- essayer("LFQ preparation", {
    l <- prepare_tropfish(dat$lfq_data,
                          bin_size    = cfg$bin,
                          species     = sp,
                          stock_name  = sp,
                          length_unit = cfg$unite,
                          Lmin        = NULL)   # no truncation
    TropFishR::lfqRestructure(l, MA = 5, addl.sqrt = FALSE)
  })
  if (is.null(lfq)) next

  grDevices::png(file.path(OUT_DIR, sprintf("%s_LFQ.png", gsub(" ", "_", sp))),
                 width = 1000, height = 600)
  try(graphics::plot(lfq, Fname = "rcounts", date.axis = "modern"), silent = TRUE)
  grDevices::dev.off()

  # --- 3.4 Growth: ELEFAN GA + SA, automatic selection ---------------
  growth <- essayer("growth (ELEFAN GA + SA)", {
    run_growth_analysis(
      lfq,
      # powell = FALSE: powell_wetherall() requires an INTERACTIVE selection of
      # reg_int; without it, it displays "Please choose the minimum and maximum
      # point in the graph" and returns an unreliable result. Z is estimated
      # properly further below, with an explicit reg_int.
      powell     = FALSE,
      methods    = c("elefan_ga", "elefan_sa"),
      Linf_range = c(cfg$Linf_min, cfg$Linf_max),
      K_range    = c(cfg$K_min,    cfg$K_max),
      criterion  = "Rn_max",
      control = list(
        low_par = list(Linf = cfg$Linf_min, K = cfg$K_min, t_anchor = 0),
        up_par  = list(Linf = cfg$Linf_max, K = cfg$K_max, t_anchor = 1),
        popSize = ELEFAN_POP, maxiter = ELEFAN_ITER, seasonalised = FALSE
      ),
      verbose = FALSE
    )
  })
  if (is.null(growth)) next

  Linf <- growth$best$Linf
  K    <- growth$best$K
  log_etape(sprintf("Linf = %.2f %s · K = %.3f /year · phi' = %.3f",
                    Linf, cfg$unite, K, growth$best$phiL))

  # Warn if the estimate hits a search bound
  if (isTRUE(growth$best$boundary_flag))
    log_etape("!! Linf or K hits a bound: WIDEN Linf_range/K_range")

  # --- 3.5 Natural mortality: ~30 methods + consensus -------------------
  M_all <- essayer("natural mortality M", {
    estimate_M_all(Linf = Linf, K = K, temp = TEMP_EAU)
  })
  M_star <- if (!is.null(M_all)) M_all$consensus[["geomean"]] else NA_real_
  if (!is.null(M_all))
    log_etape(sprintf("M (consensus, %d methods) = %.3f /year · M/K = %.2f",
                      nrow(M_all$table), M_star, M_star / K))

  # --- 3.6 Total mortality Z + F + E --------------------------------------
  ri <- REG_INT[[sp]]

  grDevices::png(file.path(OUT_DIR,
                 sprintf("%s_diagnostic_catchcurve.png", gsub(" ", "_", sp))),
                 width = 1000, height = 650)
  mort <- essayer("total mortality Z (catch curve)", {
    run_mortality(lfq,
                  growth_model = growth$models[[1]],
                  temp    = TEMP_EAU,
                  reg_int = ri,
                  plot    = TRUE)
  })
  grDevices::dev.off()

  if (!is.null(mort)) {

    ri_used <- mort$Z$reg_int
    log_etape(sprintf("reg_int = c(%d, %d)%s",
                      ri_used[1], ri_used[2],
                      if (isTRUE(mort$Z$reg_int_auto)) " (auto)" else " (imposed)"))

    log_etape(sprintf("Z = %.3f · F = %.3f · E = %.3f — %s",
                      mort$Z$Z, mort$summary$F, mort$summary$E,
                      mort$summary$statut))

    if (is.finite(mort$summary$F) && mort$summary$F < 0)
      log_etape("!! NEGATIVE F (Z < M): reg_int poorly placed, result unusable")
  }

  # --- 3.7 Length-weight allometry -----------------------------------------
  lw <- essayer("length-weight allometry", {
    fit_length_weight(dat$bio, "length", "weight")
  })
  if (!is.null(lw))
    log_etape(sprintf("W = %.5g · L^%.3f (r2 = %.3f) — %s",
                      lw$a, lw$b, lw$r2, lw$allometry))

  # --- 3.8 Maturity --------------------------------------------------------
  mat <- essayer("maturity ogive", ogive_maturite(dat$bio))
  if (!is.null(mat)) {
    log_etape(sprintf("L50 = %.2f · L95 = %.2f (%.0f %% mature in the sample)",
                      mat$L50, mat$L95, 100 * mat$prop_mature))
    if (mat$prop_mature < 0.15 || mat$prop_mature > 0.85)
      log_etape("!! poorly constrained ogive (too little contrast): L50 unreliable")
  }

  # --- 3.9 LBB (Froese script) — only if JAGS is installed -------------
  lbb <- NULL
  if (requireNamespace("R2jags", quietly = TRUE) && !is.null(mat)) {
    lbb <- essayer("LBB (JAGS)", {

      # FROESE SCRIPT UNIT CONVENTION (a trap):
      #   · length file: ALWAYS in mm
      #   · priors from the ID file: in the working unit (mm.user)
      # Stocks in cm are therefore converted to mm here, while priors
      # (Linf, Lm50) remain in cm with mm = FALSE.
      fac <- if (cfg$unite == "cm") 10 else 1

      dl <- prepare_lbb_data(
        data.frame(length = dat$lfq_data$Length * fac,
                   year   = dat$lfq_data$Year),
        "length", "year", stock = gsub(" ", "_", sp))

      run_lbb_froese(dl,
                     stock      = gsub(" ", "_", sp),
                     workdir    = file.path(OUT_DIR, "lbb", gsub(" ", "_", sp)),
                     mm         = (cfg$unite == "mm"),
                     MK_prior   = M_star / K,
                     Linf_prior = Linf,     # in the working unit
                     Lm50       = mat$L50)
    })
  } else {
    log_etape("LBB skipped (R2jags/JAGS missing or maturity unavailable)")
  }

  # --- 3.10 LBSPR ----------------------------------------------------------
  lbspr <- NULL; pars <- NULL
  if (requireNamespace("LBSPR", quietly = TRUE) && !is.null(mat) &&
      is.finite(mat$L50) && is.finite(mat$L95) && mat$L95 > mat$L50) {

    pars <- essayer("LBSPR — parameters", {
      lbspr_pars(Linf = Linf, L50 = mat$L50, L95 = mat$L95,
                 M = M_star, K = K, species = sp, bin_width = cfg$bin)
    })

    if (!is.null(pars)) {
      lbspr <- essayer("LBSPR — fitting", {
        ln <- lbspr_lengths(
          data.frame(length = dat$lfq_data$Length, year = dat$lfq_data$Year),
          pars, length_col = "length", year_col = "year", bin_width = cfg$bin)
        run_lbspr(pars, ln, spr_target = SPR_CIBLE, spr_limit = SPR_LIMITE)
      })
      if (!is.null(lbspr)) {
        last <- lbspr$summary[nrow(lbspr$summary), ]
        log_etape(sprintf("SPR = %.3f — %s", last$SPR, last$statut))
      }
    }
  }

  # --- 3.11 Management measures --------------------------------------------
  mse <- NULL
  if (!is.null(pars)) {
    mse <- essayer("management measures (equilibrium MSE)", {
      scen <- list(
        mse_measure("Status quo",       "statuquo"),
        mse_measure("-20 % effort",     "effort",           effort_reduction = 0.20),
        mse_measure("-40 % effort",     "effort",           effort_reduction = 0.40),
        mse_measure("2-month closure",  "repos_biologique", closure_months   = 2),
        mse_measure("3-month closure",  "repos_biologique", closure_months   = 3),
        mse_measure("MPA 20 %",         "amp",              mpa_fraction     = 0.20),
        mse_measure("MLS = L50",        "taille_min",       Lc = mat$L50),
        mse_optimal_size(Linf = Linf, MK = M_star / K)
      )
      run_mse_equilibrium(pars, scen, base_FM = BASE_FM,
                          spr_target = SPR_CIBLE, spr_limit = SPR_LIMITE)
    })
  }

  # --- 3.12 HTML report ---------------------------------------------------
  res <- collect_results(
    growth = growth, mortality = mort, allometry = lw,
    maturity = mat, lbb = lbb, lbspr = lbspr, mse = mse,
    meta = list(
      espece    = sp,
      stock     = "Example EEZ",
      periode   = paste(range(dat$lfq_data$Year), collapse = "-"),
      n_mesures = nrow(dat$lfq_data),
      unite     = cfg$unite
    )
  )
  resultats[[sp]] <- res

  essayer("HTML report", {
    fishstock_report(
      res,
      file   = file.path(OUT_DIR, sprintf("rapport_%s.html", gsub(" ", "_", sp))),
      format = "html",
      title  = paste("Stock assessment —", sp),
      author = "Kamarel Ba — CRODT"
    )
  })

  # --- 3.13 Summary row ---------------------------------------------
  synthese[[sp]] <- data.frame(
    Stock  = sp,
    Linf   = round(Linf, 2),
    K      = round(K, 3),
    M      = round(M_star, 3),
    Z      = if (!is.null(mort)) round(mort$Z$Z, 3)        else NA_real_,
    F_     = if (!is.null(mort)) round(mort$summary$F, 3)  else NA_real_,
    E      = if (!is.null(mort)) round(mort$summary$E, 3)  else NA_real_,
    L50    = if (!is.null(mat))  round(mat$L50, 2)         else NA_real_,
    SPR    = if (!is.null(lbspr)) round(lbspr$summary$SPR[nrow(lbspr$summary)], 3) else NA_real_,
    Statut = if (!is.null(mort)) mort$summary$statut       else NA_character_,
    stringsAsFactors = FALSE
  )
}


###############################################################################
#  BLOCK 4 — CROSS-STOCK SUMMARY
###############################################################################

cat("\n", strrep("=", 74), "\n", sep = "")
cat("  SUMMARY — ", length(synthese), " stocks assessed\n", sep = "")
cat(strrep("=", 74), "\n\n", sep = "")

if (length(synthese)) {

  tab <- do.call(rbind, synthese)
  rownames(tab) <- NULL
  print(tab, row.names = FALSE)

  utils::write.csv(tab, file.path(OUT_DIR, "synthese_stocks.csv"),
                   row.names = FALSE)

  # Exploitation diagram: E vs target 0.5 (Gulland)
  if (requireNamespace("ggplot2", quietly = TRUE) && any(is.finite(tab$E))) {
    p <- ggplot2::ggplot(tab[is.finite(tab$E), ],
                         ggplot2::aes(x = stats::reorder(Stock, E), y = E)) +
      ggplot2::geom_col(fill = "#2c7fb8") +
      ggplot2::geom_hline(yintercept = 0.5, linetype = 2, colour = "#c0392b") +
      ggplot2::annotate("text", x = 0.7, y = 0.52,
                        label = "E = 0.5 (Gulland optimum)",
                        hjust = 0, size = 3, colour = "#c0392b") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "Exploitation rate by stock",
                    subtitle = "E > 0.5: overexploitation",
                    x = NULL, y = "E = F/Z") +
      ggplot2::theme_bw()
    ggplot2::ggsave(file.path(OUT_DIR, "synthese_exploitation.png"),
                    p, width = 8, height = 5, dpi = 150)
  }

  cat("\nOutputs written to:", normalizePath(OUT_DIR), "\n")
  cat("  · rapport_<stock>.html            one report per stock\n")
  cat("  · synthese_stocks.csv             comparative table\n")
  cat("  · synthese_exploitation.png       exploitation rate\n")
  cat("  · <stock>_diagnostic_catchcurve.png  ← CHECK reg_int on it\n")

} else {
  cat("No stock could be assessed. Check DATA_DIR and the data.\n")
}

cat(sprintf("\nTotal duration: %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

###############################################################################
#  IMPORTANT REMINDERS
#
#  1. reg_int — The default values are STARTING POINTS. Open the
#     *_diagnostic_catchcurve.png files and adjust BLOCK 1. Z (hence F and E,
#     hence stock status) depends directly on this choice.
#
#  2. boundary_flag — If the log shows "hits a bound", ELEFAN converged
#     on the edge of the search space: the estimate is not reliable.
#     Widen Linf_min/max or K_min/max for this stock.
#
#  3. Sampling duration — This data covers ~1 year. ELEFAN tracks
#     cohorts over time: with only one year, K estimates are
#     fragile, and LBB (which works by year) only has 2 points. Results
#     are ORDERS OF MAGNITUDE, to be consolidated with 2-3 years of data.
#
#  4. Volutes — The 4 species are assessed separately (very different Linf).
#     Do not group them into the same LFQ.
###############################################################################