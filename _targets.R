###############################################################
#
# stockflow — reproducible {targets} pipeline
#
# Full chain: import -> validation -> LFQ -> growth ->
# mortality -> allometry -> maturity -> LBB -> LBSPR -> management
# measures -> automatic report.
#
# Usage:
#   targets::tar_make()          # runs only what has changed
#   targets::tar_visnetwork()    # interactive DAG
#   targets::tar_read(brp)       # read a target
#
# The parameters (paths, priors, ELEFAN bounds) are centralized
# in config.yml.
###############################################################

library(targets)   # required by targets to define the pipeline

tar_option_set(
  packages = c("stockflow"),
  format   = "rds",
  error    = "continue"   # a failed target does not stop everything
)

# ---- Configuration -----------------------------------------------------
cfg <- if (file.exists("config.yml")) {
  yaml::read_yaml("config.yml")
} else {
  list()
}

get_cfg <- function(path, default) {
  x <- cfg
  for (p in strsplit(path, "\\.")[[1]]) {
    if (is.null(x[[p]])) return(default)
    x <- x[[p]]
  }
  x
}

FILE_LENGTHS <- get_cfg("data.lengths", "inst/extdata/penaeus_lengths_demo.csv")
FILE_BIO     <- get_cfg("data.biology", "inst/extdata/penaeus_bio_demo.csv")
SPECIES      <- get_cfg("species.name", "Penaeus notialis")
BIN          <- get_cfg("sampling.bin_size", 1)
LINF_RANGE   <- c(get_cfg("biology.Linf_min", 40), get_cfg("biology.Linf_max", 55))
K_RANGE      <- c(get_cfg("biology.K_min", 0.4),  get_cfg("biology.K_max", 2.0))
TEMP         <- get_cfg("biology.temp", 27)
SPR_TARGET   <- get_cfg("management.SPR_target", 0.40)
CATCH_CURRENT <- get_cfg("management.catch_current", NA)   # recent catch (t) for the TAC
TAC_BUFFER    <- get_cfg("management.tac_buffer", 0.10)     # precautionary coefficient
FILE_PRORATA  <- get_cfg("management.prorata_file", NA)     # historical catch CSV (segment)
PRORATA_SEG   <- get_cfg("management.prorata_segment", "espece")
PRORATA_CATCH <- get_cfg("management.prorata_catch", "capture_t")
CLOSURE_DUR   <- get_cfg("management.closure_duration", NULL)  # NULL = compare 1..6 months
MATURE_CODES  <- get_cfg("biology.mature_codes", NULL)         # "mature" maturity codes


list(

  # ---- 1. Import (tracking the source file) -----------------------------
  tar_target(file_lengths, FILE_LENGTHS, format = "file"),
  tar_target(file_bio,     FILE_BIO,     format = "file"),

  tar_target(lengths, read_lengths(file_lengths)),
  tar_target(bio,     utils::read.csv(file_bio)),

  # ---- 2. Validation ---------------------------------------------------
  tar_target(validation,
             validate_lengths(lengths, species_config = list(),
                              min_length = 0, max_length = 60)),

  # ---- 3. LFQ object (TropFishR) ----------------------------------------
  tar_target(lfq, {
    d <- lengths
    names(d)[names(d) == "length"] <- "Length"
    names(d)[names(d) == "year"]   <- "Year"
    names(d)[names(d) == "month"]  <- "Month"
    prepare_tropfish(d, bin_size = BIN, species = SPECIES,
                     stock_name = SPECIES, length_unit = "cm")
  }),

  # ---- 4. Growth (ELEFAN GA + SA, automatic selection) -----------
  tar_target(growth,
             run_growth_analysis(
               lfq,
               methods    = c("elefan_ga", "elefan_sa"),
               Linf_range = LINF_RANGE,
               K_range    = K_RANGE,
               criterion  = "Rn_max",
               control = list(
                 low_par = list(Linf = LINF_RANGE[1], K = K_RANGE[1], t_anchor = 0),
                 up_par  = list(Linf = LINF_RANGE[2], K = K_RANGE[2], t_anchor = 1),
                 popSize = 60, maxiter = 50, seasonalised = FALSE
               ),
               verbose = FALSE)),

  tar_target(growth_par, growth$best),

  # ---- 5. Natural mortality (all methods + consensus) ------------
  tar_target(M_all,
             estimate_M_all(Linf = growth_par$Linf,
                            K    = growth_par$K,
                            temp = TEMP)),

  tar_target(M_star, M_all$consensus[["geomean"]]),

  # ---- 6. Length-weight allometry --------------------------------------
  tar_target(allometry, fit_length_weight(bio, "length", "weight")),

  # ---- 7. maturity (L50 / L95) -----------------------------------------
  tar_target(maturity, {
    b <- bio
    b$mature <- as.integer(b$maturity >= 3)
    g  <- stats::glm(mature ~ length, data = b, family = stats::binomial())
    cf <- stats::coef(g)
    list(L50 = as.numeric(-cf[1] / cf[2]),
         L95 = as.numeric((log(0.95 / 0.05) - cf[1]) / cf[2]))
  }),

  # ---- 8. LBB (Froese backend script; requires R2jags + JAGS) --------
  tar_target(lbb_data,
             prepare_lbb_data(lengths, "length", "year",
                              stock = gsub(" ", "_", SPECIES))),

  tar_target(lbb, {
    if (!requireNamespace("R2jags", quietly = TRUE)) NULL
    else run_lbb_froese(lbb_data,
                        mm = TRUE,
                        MK_prior   = M_star / growth_par$K,
                        Linf_prior = growth_par$Linf,
                        Lm50       = maturity$L50)
  }),

  # ---- 9. LBSPR ---------------------------------------------------------
  tar_target(lbspr_par,
             lbspr_pars(Linf = growth_par$Linf,
                        L50  = maturity$L50,
                        L95  = maturity$L95,
                        M    = M_star,
                        K    = growth_par$K,
                        species = SPECIES)),

  tar_target(lbspr_len,
             lbspr_lengths(lengths, lbspr_par,
                           length_col = "length", year_col = "year",
                           bin_width = BIN)),

  tar_target(lbspr,
             run_lbspr(lbspr_par, lbspr_len,
                       spr_target = SPR_TARGET, spr_limit = 0.20)),

  # ---- 10. Management measures (equilibrium MSE) ------------------------
  tar_target(scenarios, list(
    mse_measure("Statu quo",     "statuquo"),
    mse_measure("-30% effort",   "effort",           effort_reduction = 0.30),
    mse_measure("Repos 3 mois",  "repos_biologique", closure_months   = 3),
    mse_measure("AMP 20%",       "amp",              mpa_fraction     = 0.20),
    mse_measure("Taille min L50", "taille_min",      Lc = maturity$L50),
    mse_optimal_size(Linf = growth_par$Linf,
                     MK   = M_star / growth_par$K)
  )),

  tar_target(mse,
             run_mse_equilibrium(lbspr_par, scenarios,
                                 base_FM = 1.5,
                                 spr_target = SPR_TARGET)),

  # ---- 11. YPR (yield per recruit) + reference points ------------
  tar_target(ypr_input,
             ypr_param(Linf = growth_par$Linf, K = growth_par$K,
                       M = M_star, a = allometry$a, b = allometry$b,
                       Lc = maturity$L50, Lr = maturity$L50 * 0.5)),

  tar_target(ypr,
             tryCatch(run_ypr(ypr_input,
                              FM_change = seq(0, 3, 0.05),
                              curr.E = 0.5, curr.Lc = maturity$L50),
                      error = function(e) NULL)),

  # ---- 12. TAC + quota allocation (objective i) -------------------
  tar_target(tac, {
     e_target <- if (!is.null(ypr) && !is.null(ypr$reference_points))
      ypr$reference_points$E01[1] else NA
    e_curr <- if (!is.null(ypr) && !is.null(ypr$currents))
      ypr$currents$curr.E[1] else NA
    if (is.na(CATCH_CURRENT) || is.na(e_target) || is.na(e_curr)) NULL
    else compute_tac(method = "ratio_E", catch_current = CATCH_CURRENT,
                     E_current = e_curr, E_target = e_target,
                     buffer = TAC_BUFFER, label = SPECIES)
  }),

  tar_target(quota, {
    if (is.null(tac) || is.na(FILE_PRORATA) || !file.exists(FILE_PRORATA))
      NULL
    else {
      cap <- utils::read.csv(FILE_PRORATA)
      pr  <- compute_prorata(cap, segment_col = PRORATA_SEG,
                             catch_col = PRORATA_CATCH)
      allocate_quota(tac = tac, prorata = pr)
    }
  }),

  # ---- 13. Optimal biological rest period (objective ii) ---------------------
  tar_target(rec_pattern,
             tryCatch({
               l <- lfq
               l$Linf <- growth_par$Linf; l$K <- growth_par$K; l$t0 <- 0
               recruitment_pattern(l)
             }, error = function(e) NULL)),

  tar_target(spawning, {
    if (is.null(MATURE_CODES) || !"maturity" %in% names(bio)) NULL
    else tryCatch(
      spawning_season(bio, date_col = "date", maturity_col = "maturity",
                      mature_codes = MATURE_CODES),
      error = function(e) NULL)
  }),

  tar_target(closure, {
    if (is.null(rec_pattern) && is.null(spawning)) NULL
    else optimal_closure(spawning = spawning, recruitment = rec_pattern,
                         duration = CLOSURE_DUR)
  }),

  # ---- 14. Technical measures: management sizes (objective iii) -----
  tar_target(size_advice,
             recommend_sizes(candidates = seq(maturity$L50 * 0.7,
                                              growth_par$Linf * 0.7,
                                              length.out = 6),
                             Linf = growth_par$Linf,
                             L50_maturity = maturity$L50,
                             MK = M_star / growth_par$K)),

  # ---- 15. Aggregation + report ----------------------------------------
  tar_target(results,
             collect_results(
               growth      = growth,
               allometry   = allometry,
               maturity    = maturity,
               lbb         = lbb,
               lbspr       = lbspr,
               mse         = mse,
               ypr         = ypr,
               tac         = tac,
               quota       = quota,
               closure     = closure,
               size_advice = size_advice,
               meta = list(species = SPECIES,
                           n_ind   = nrow(lengths),
                           period  = paste(range(lengths$year), collapse = "-")))),

  tar_target(report,
             fishstock_report(results,
                              file   = "output/rapport_stockflow.html",
                              format = "html",
                              title  = paste("Stock assessment -", SPECIES),
                              author = "stockflow"),
             format = "file")

)