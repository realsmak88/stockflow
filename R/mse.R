###############################################################
#
# stockflow
#
# Module : mse.R
#
# Simulation and comparison of management measures.
#
# Two engines, wrappers of OFFICIAL functions:
#   (A) Per-recruit equilibrium (length-based):
#         LBSPR::LBSPRsim()      -> SPR, relative yield, SSB
#         TropFishR::predict_mod() -> Thompson-Bell (absolute yield)
#   (B) Stochastic closed loop (optional):
#         MSEtool::runMSE()      -> openMSE MSE
#
# Each management measure (biological rest period, minimum size,
# mesh size, optimal size, effort reduction, MPA, quota) is
# translated into selectivity modifiers (SL50/SL95) and/or
# the F/M ratio. No scientific algorithm is rewritten.
#
###############################################################

#==============================================================
# Dependencies
#==============================================================

.check_lbspr_mse <- function() {

  if (!requireNamespace("LBSPR", quietly = TRUE))
    stop("The 'LBSPR' package must be installed (equilibrium engine).",
         call. = FALSE)

  if (!requireNamespace("methods", quietly = TRUE))
    stop("The 'methods' package is required (S4 objects).", call. = FALSE)

  invisible(TRUE)

}


#==============================================================
# Definition of a management measure
#==============================================================

#' Definition of a management measure (scenario)
#'
#' Translates a management measure into modifiers applied to selectivity
#' (\code{SL50}, \code{SL95}) and/or to the \eqn{F/M} ratio, relative to a
#' reference state. Spatio-temporal closures and MPAs are represented,
#' to first order, as a proportional reduction in \eqn{F}.
#'
#' @param label Scenario name.
#' @param type Type of measure: \code{"statuquo"}, \code{"effort"},
#'   \code{"repos_biologique"}, \code{"amp"}, \code{"taille_min"},
#'   \code{"maillage"}, \code{"taille_optimale"}, \code{"quota"}.
#' @param effort_reduction Fraction of effort reduction (0-1), for
#'   \code{type = "effort"}.
#' @param closure_months Closure duration in months (0-12), for
#'   \code{type = "repos_biologique"}.
#' @param mpa_fraction Fraction of the zone/stock protected (0-1), for
#'   \code{type = "amp"}.
#' @param Lc New size at first capture (SL50), for
#'   \code{"taille_min"} / \code{"maillage"}.
#' @param SL95 New SL95 (default \code{1.1 * Lc}).
#' @param TAC Target catch (t), for \code{type = "quota"} (evaluated in closed
#'   loop, see \code{run_mse_closed_loop()}).
#'
#' @return List of class \code{mse_measure}.
#' @export
mse_measure <- function(label,
                        type = c("statuquo", "effort", "repos_biologique",
                                 "amp", "taille_min", "maillage",
                                 "taille_optimale", "quota"),
                        effort_reduction = 0,
                        closure_months = 0,
                        mpa_fraction = 0,
                        Lc = NULL,
                        SL95 = NULL,
                        TAC = NULL) {

  type <- match.arg(type)

  structure(
    list(label = label, type = type,
         effort_reduction = effort_reduction,
         closure_months = closure_months,
         mpa_fraction = mpa_fraction,
         Lc = Lc, SL95 = SL95, TAC = TAC),
    class = "mse_measure"
  )

}


## Translates a measure into (FM_mult, SL50, SL95) from the reference state.
.apply_measure <- function(m, base_FM, base_SL50, base_SL95) {

  FM_mult <- 1
  SL50 <- base_SL50
  SL95 <- base_SL95

  if (m$type == "effort")
    FM_mult <- 1 - m$effort_reduction

  if (m$type == "repos_biologique")
    FM_mult <- 1 - m$closure_months / 12          # 1st order: fishing time

  if (m$type == "amp")
    FM_mult <- 1 - m$mpa_fraction                 # 1st order: protected fraction

  if (m$type %in% c("taille_min", "maillage", "taille_optimale")) {
    if (is.null(m$Lc))
      stop("Measure '", m$label, "' : 'Lc' required.")
    SL50 <- m$Lc
    SL95 <- if (!is.null(m$SL95)) m$SL95 else m$Lc * 1.1
  }

  list(FM = base_FM * FM_mult, SL50 = SL50, SL95 = SL95,
       FM_mult = FM_mult)

}


#' "Optimal size" measure (Lc = Froese's Lopt)
#'
#' Constructs a mesh-size-type measure targeting the optimal size at first
#' capture \eqn{L_{opt}} derived from \eqn{L_{inf}} and \eqn{M/K}
#' (\code{lbb_reference_points()}).
#'
#' @param Linf,MK Parameters for \eqn{L_{opt} = 3/(3+M/K)\cdot L_{inf}}.
#' @param FM Current \eqn{F/M} ratio (for \eqn{L_{c,opt}}).
#' @param label Scenario name.
#' @return \code{mse_measure} object.
#' @export
mse_optimal_size <- function(Linf, MK, FM = 1, label = "Taille optimale") {

  rp <- lbb_reference_points(Linf = Linf, MK = MK, FM = FM)

  mse_measure(label, type = "taille_optimale",
              Lc = rp$Lc_opt, SL95 = rp$Lc_opt * 1.1)

}


#==============================================================
# Equilibrium engine: LBSPR::LBSPRsim
#==============================================================

.lbsim_indicators <- function(pars) {

  sim <- LBSPR::LBSPRsim(pars)

  getslot <- function(nm)
    tryCatch(as.numeric(methods::slot(sim, nm))[1],
             error = function(e) NA_real_)

  list(
    SPR   = getslot("SPR"),
    Yield = getslot("Yield"),   # relative yield (0-1)
    YPR   = getslot("YPR"),
    SSB   = getslot("SSB")
  )

}


#' Comparison of management measures at equilibrium (LBSPR)
#'
#' Evaluates, for each measure, the SPR, relative yield and relative spawning
#' biomass at equilibrium via \code{LBSPR::LBSPRsim()}, and diagnoses the
#' status relative to SPR targets.
#'
#' @param pars \code{LB_pars} object (see \code{lbspr_pars()}) describing the
#'   biology (Linf, L50, L95, MK).
#' @param measures List of \code{mse_measure} objects.
#' @param base_FM Reference \eqn{F/M} ratio (current state).
#' @param base_SL50,base_SL95 Reference selectivity (SL50 defaults to
#'   \code{pars@@L50}).
#' @param spr_target,spr_limit SPR targets (0.40 / 0.20).
#'
#' @return S3 object \code{FishStockMSE}: \code{comparison} (data.frame:
#'   scenario, FM, SL50, SPR, Yield_rel, SSB_rel, status) and \code{targets}.
#'
#' @references
#'   Hordyk, A. et al. (2015). ICES J. Mar. Sci. 72(1).
#'   Froese, R. et al. (2016). Fish and Fisheries 17(3).
#' @examples
#' \dontrun{
#'   pr <- lbspr_pars(Linf = 45, L50 = 22, L95 = 28, M = 1.5, K = 1.2)
#'   sc <- list(
#'     mse_measure("Status quo", "statuquo"),
#'     mse_measure("-30\% effort", "effort", effort_reduction = 0.3),
#'     mse_measure("3-month rest", "repos_biologique", closure_months = 3),
#'     mse_measure("Min size 25", "taille_min", Lc = 25)
#'   )
#'   res <- run_mse_equilibrium(pr, sc, base_FM = 1.5)
#'   res$comparison
#' }
#' @export
run_mse_equilibrium <- function(pars,
                                measures,
                                base_FM,
                                base_SL50 = NULL,
                                base_SL95 = NULL,
                                spr_target = 0.40,
                                spr_limit = 0.20) {

  .check_lbspr_mse()

  if (!methods::is(pars, "LB_pars"))
    stop("'pars' must be an LB_pars object (see lbspr_pars()).")

  if (is.null(base_SL50)) base_SL50 <- pars@L50
  if (is.null(base_SL95)) base_SL95 <- pars@L95


  rows <- lapply(measures, function(m) {

    if (!inherits(m, "mse_measure"))
      stop("Each element of 'measures' must be an mse_measure object.")

    mod <- .apply_measure(m, base_FM, base_SL50, base_SL95)

    p2 <- pars
    p2@SL50 <- mod$SL50
    p2@SL95 <- mod$SL95
    p2@FM   <- mod$FM

    ind <- tryCatch(.lbsim_indicators(p2),
                    error = function(e)
                      list(SPR = NA_real_, Yield = NA_real_,
                           YPR = NA_real_, SSB = NA_real_))

    data.frame(
      scenario  = m$label,
      type      = m$type,
      FM        = mod$FM,
      SL50      = mod$SL50,
      SPR       = ind$SPR,
      Yield_rel = ind$Yield,
      SSB_rel   = ind$SSB,
      statut    = .lbspr_status(ind$SPR, spr_target, spr_limit),
      stringsAsFactors = FALSE
    )

  })

  comparison <- do.call(rbind, rows)
  rownames(comparison) <- NULL


  structure(
    list(comparison = comparison,
         targets = c(SPR_target = spr_target, SPR_limit = spr_limit)),
    class = "FishStockMSE"
  )

}


#==============================================================
# Absolute equilibrium engine: Thompson-Bell (TropFishR)
#==============================================================

#' Thompson-Bell yield prediction (TropFishR)
#'
#' Wrapper around \code{TropFishR::predict_mod(type = "ThompBell")} to
#' evaluate the effect of a change in effort (\code{FM_change}) and/or size
#' at first capture (\code{Lc_change}) on yield and biomass.
#'
#' @param param List of TropFishR-compatible parameters (Linf, K, M, a, b,
#'   plus composition/selectivity depending on version).
#' @param FM_change Vector of F multipliers to explore.
#' @param Lc_change Vector of sizes at first capture to explore
#'   (optional, for a mesh-size x effort isopleth analysis).
#' @param ... Arguments passed to \code{TropFishR::predict_mod()}.
#' @return Raw output of \code{TropFishR::predict_mod()}.
#' @export
run_thompson_bell <- function(param,
                              FM_change = seq(0, 2, 0.1),
                              Lc_change = NULL,
                              ...) {

  if (!requireNamespace("TropFishR", quietly = TRUE))
    stop("The 'TropFishR' package must be installed.", call. = FALSE)

  args <- c(list(param, type = "ThompBell", FM_change = FM_change),
            if (!is.null(Lc_change)) list(Lc_change = Lc_change),
            list(...))

  tryCatch(
    do.call(TropFishR::predict_mod, args),
    error = function(e)
      stop("predict_mod : ", conditionMessage(e), call. = FALSE)
  )

}


#==============================================================
# openMSE bridge (stochastic closed loop) - optional
#==============================================================

#' Closed-loop MSE (openMSE / MSEtool) - optional bridge
#'
#' Minimal wrapper around \code{MSEtool::runMSE()} to test management
#' procedures (management procedures: quotas, TAC, effort) under
#' uncertainty. Requires the \code{MSEtool} package (openMSE).
#'
#' @param OM \code{OM} (operating model) MSEtool object.
#' @param MPs Vector of management procedures (\code{MSEtool::avail("MP")}).
#' @param ... Arguments passed to \code{MSEtool::runMSE()}.
#' @return MSEtool MSE object.
#' @references Carruthers, T.R., Hordyk, A.R. (2018). openMSE. MEE 9(12).
#' @export
run_mse_closed_loop <- function(OM,
                                MPs = c("AvC", "DCAC", "Itarget1",
                                        "matlenlim", "curE"),
                                ...) {

  if (!requireNamespace("MSEtool", quietly = TRUE))
    stop("The 'MSEtool' package (openMSE) must be installed.", call. = FALSE)

  tryCatch(
    MSEtool::runMSE(OM, MPs = MPs, ...),
    error = function(e)
      stop("runMSE : ", conditionMessage(e), call. = FALSE)
  )

}


#==============================================================
# S3 methods
#==============================================================

#' @export
print.FishStockMSE <- function(x, ...) {

  cat("\n=====================================\n")
  cat(" stockflow : Comparison of management measures (equilibrium)\n")
  cat("=====================================\n\n")

  cat(sprintf("SPR targets: target %.2f, limit %.2f\n\n",
              x$targets["SPR_target"], x$targets["SPR_limit"]))

  print(x$comparison, row.names = FALSE)

  invisible(x)

}


#' Relative yield vs SPR trade-off by scenario
#' @param x \code{FishStockMSE} object.
#' @param ... Ignored.
#' @export
plot.FishStockMSE <- function(x, ...) {

  df <- x$comparison

  if (!requireNamespace("ggplot2", quietly = TRUE)) {

    graphics::plot(df$SPR, df$Yield_rel, pch = 19,
                   xlab = "SPR", ylab = "Relative yield")
    graphics::text(df$SPR, df$Yield_rel, df$scenario, pos = 3, cex = 0.7)
    graphics::abline(v = x$targets["SPR_target"], lty = 2, col = "darkgreen")
    graphics::abline(v = x$targets["SPR_limit"],  lty = 3, col = "red")
    return(invisible(x))

  }

  p <- ggplot2::ggplot(df, ggplot2::aes(SPR, Yield_rel, label = scenario)) +
    ggplot2::annotate("rect",
                      xmin = 0, xmax = x$targets["SPR_limit"],
                      ymin = -Inf, ymax = Inf, fill = "#FF7F7F", alpha = .18) +
    ggplot2::geom_vline(xintercept = x$targets["SPR_target"],
                        linetype = 2, colour = "darkgreen") +
    ggplot2::geom_point(size = 3, colour = "#12507b") +
    ggplot2::geom_text(vjust = -0.8, size = 3) +
    ggplot2::labs(title = "Management measures: yield vs reproductive potential",
                  x = "SPR (reproductive potential)",
                  y = "Relative yield") +
    ggplot2::theme_bw()

  print(p)

  invisible(p)

}