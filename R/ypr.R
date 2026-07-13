###############################################################
#
# stockflow
#
# Module : ypr.R
#
# Yield-Per-Recruit (YPR) analysis and associated biological
# reference points.
#
# Wrappers around OFFICIAL functions only:
#   TropFishR::predict_mod(type = "ypr")       -> Beverton & Holt YPR
#   TropFishR::predict_mod(type = "ThompBell") -> Thompson & Bell
#
# The module does NOT rewrite ANY scientific algorithm. It:
#   (1) assembles the parameter list expected by predict_mod()
#       from the objects already produced by the pipeline
#       (FishStockMortality, FishStockLW, growth);
#   (2) runs the official model;
#   (3) extracts and compares the reference points (F0.1, Fmax,
#       E0.1, Emax, E0.5) that predict_mod() computes itself.
#
# Reference :
#   Beverton, R.J.H. & Holt, S.J. (1957). On the dynamics of
#     exploited fish populations. Fish. Invest. Ser. II, 19.
#   Sparre, P. & Venema, S.C. (1998). Introduction to tropical
#     fish stock assessment. FAO Fish. Tech. Pap. 306/1.
#   Mildenberger, T.K., Taylor, M.H. & Wolff, M. (2017).
#     TropFishR. Methods Ecol. Evol. 8(11), 1520-1527.
#
###############################################################


#==============================================================
# Dependencies
#==============================================================

.check_tropfishr <- function() {

  if (!requireNamespace("TropFishR", quietly = TRUE))
    stop("The 'TropFishR' package must be installed (per-recruit engine).",
         call. = FALSE)

  invisible(TRUE)

}


## Extracts Linf, K, t0 from a heterogeneous growth object (list,
## FishStockGrowth, FishStockVBGF, ELEFAN output...). Returns a
## named list; missing fields are NULL.
.ypr_growth_pars <- function(growth) {

  if (is.null(growth))
    return(list(Linf = NULL, K = NULL, t0 = NULL))

  ## A growth object from the package often exposes $par or $Linf/$K.
  cand <- growth
  if (!is.null(growth$par) && is.list(growth$par)) cand <- growth$par
  if (!is.null(growth$growth) && is.list(growth$growth)) cand <- growth$growth

  getv <- function(nm) {
    for (k in nm) if (!is.null(cand[[k]]) && is.finite(suppressWarnings(
      as.numeric(cand[[k]])[1]))) return(as.numeric(cand[[k]])[1])
    NULL
  }

  list(
    Linf = getv(c("Linf", "linf", "L_inf")),
    K    = getv(c("K", "k")),
    t0   = getv(c("t0", "t_0", "tzero"))
  )

}


#==============================================================
# Assembling parameters for predict_mod()
#==============================================================

#' Assemble the parameters for a per-recruit analysis
#'
#' Builds the parameter list expected by
#' \code{TropFishR::predict_mod()}, preferably from the objects already
#' produced by the \pkg{stockflow} pipeline (mortality, allometry,
#' growth). Any parameter passed explicitly takes precedence over the
#' value inferred from the objects.
#'
#' @param mortality \code{FishStockMortality} object (from
#'   \code{run_mortality()}) providing \eqn{M} and the growth
#'   parameters. Optional if \code{Linf}, \code{K} and \code{M} are provided.
#' @param allometry \code{FishStockLW} object (from
#'   \code{fit_length_weight()}) providing \eqn{a} and \eqn{b}. Optional if
#'   \code{a} and \code{b} are provided.
#' @param growth Growth object (ELEFAN, VBGF...) providing
#'   \eqn{L_\infty}, \eqn{K}, \eqn{t_0}. Optional.
#' @param Linf,K,t0,M Growth and natural mortality parameters. Take
#'   precedence over the objects if provided.
#' @param a,b Coefficients of the length-weight relationship \eqn{W = a L^b}.
#' @param Lc,Lr Length at first capture (\eqn{L_c}, knife-edge
#'   selectivity) and recruitment length (\eqn{L_r}). \code{Lc} is required
#'   for the length-based per-recruit engine.
#' @param bin_size Width of the length classes (cm) for the discretization
#'   in \code{predict_mod()} (default 1).
#'
#' @return Named list of class \code{ypr_param}, directly usable by
#'   \code{run_ypr()} and \code{run_thompson_bell()}.
#' @seealso \code{\link{run_ypr}}, \code{\link{run_thompson_bell}}
#' @export
ypr_param <- function(mortality = NULL,
                      allometry = NULL,
                      growth = NULL,
                      Linf = NULL, K = NULL, t0 = NULL, M = NULL,
                      a = NULL, b = NULL,
                      Lc = NULL, Lr = NULL,
                      bin_size = 1) {

  ## --- Growth: priority to mortality$growth object, then growth ---
  gp_mort <- if (!is.null(mortality)) .ypr_growth_pars(mortality$growth) else
    list(Linf = NULL, K = NULL, t0 = NULL)
  gp_grow <- .ypr_growth_pars(growth)

  pick <- function(x, ...) {
    for (v in list(x, ...)) if (!is.null(v) && is.finite(v)) return(v)
    NULL
  }

  Linf_v <- pick(Linf, gp_mort$Linf, gp_grow$Linf)
  K_v    <- pick(K,    gp_mort$K,    gp_grow$K)
  t0_v   <- pick(t0,   gp_mort$t0,   gp_grow$t0)

  ## --- Natural mortality: argument, then object consensus ---
  M_v <- pick(M,
              if (!is.null(mortality)) mortality$M_consensus else NULL)

  ## --- Allometry ---
  a_v <- pick(a, if (!is.null(allometry)) allometry$a else NULL)
  b_v <- pick(b, if (!is.null(allometry)) allometry$b else NULL)

  ## --- Checks ---
  if (is.null(Linf_v) || is.null(K_v))
    stop("Linf and K are required (via 'mortality'/'growth' or as arguments).",
         call. = FALSE)
  if (is.null(M_v))
    stop("Natural mortality M is required (via 'mortality' or M=).",
         call. = FALSE)

  param <- list(
    Linf = Linf_v, K = K_v, t0 = if (is.null(t0_v)) 0 else t0_v,
    M = M_v, a = a_v, b = b_v,
    Lr = Lr, Lc = Lc,
    bin_size = bin_size
  )
  param <- param[!vapply(param, is.null, logical(1))]

  structure(param, class = c("ypr_param", "list"))

}


#' @export
print.ypr_param <- function(x, ...) {

  cat("\n<ypr_param> parameters for the per-recruit analysis\n")
  cat(sprintf("  Linf = %.2f | K = %.3f | t0 = %.3f | M = %.3f\n",
              x$Linf %||% NA, x$K %||% NA,
              (x$t0 %||% NA), x$M %||% NA))
  if (!is.null(x$a) && !is.null(x$b))
    cat(sprintf("  a = %.4g | b = %.3f (W = a L^b)\n", x$a, x$b))
  if (!is.null(x$Lc)) cat(sprintf("  Lc = %.2f\n", x$Lc))
  invisible(x)

}

## small internal fallback operator (not exported)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a


#==============================================================
# Per-recruit engine: Beverton & Holt (TropFishR)
#==============================================================

#' Beverton & Holt per-recruit analysis (YPR)
#'
#' Wrapper around \code{TropFishR::predict_mod(type = "ypr")}. Explores a
#' range of fishing mortalities (\code{FM_change}, expressed as \eqn{F/M} or
#' directly as \eqn{F} depending on the parameterization) and, optionally, a
#' range of lengths at first capture (\code{Lc_change}) to produce a mesh
#' \eqn{\times} effort isopleth analysis.
#'
#' @param param \code{ypr_param} object (see \code{ypr_param()}) or a list
#'   compatible with \code{TropFishR::predict_mod()}. In length-based mode,
#'   the official function requires a recruitment length \code{Lr}: it is
#'   taken from \code{param$Lr} or, failing that, initialized to the smallest
#'   \code{Lc} explored.
#' @param FM_change Vector of effort values to explore.
#' @param Lc_change Vector of lengths at first capture. If \code{NULL} and
#'   \code{param$Lc} exists, this single value is used.
#' @param s_list Selectivity list for \code{predict_mod()} (e.g.
#'   \code{list(selecType = "trawl_ogive", L50 = , L75 = )}). If \code{NULL},
#'   a trawl ogive is built from the first \code{Lc}
#'   (\eqn{L_{50} = L_c}, \eqn{L_{75} = 1.08\,L_c}).
#' @param FM_relative Logical. The YPR engine of \code{predict_mod()} only
#'   accepts absolute \eqn{F} values; this parameter therefore defaults to
#'   \code{FALSE} and any \code{TRUE} value is reset to \code{FALSE} with a
#'   warning.
#' @param plot Logical: whether to plot the \code{predict_mod()} diagram.
#' @param ... Additional arguments passed to
#'   \code{TropFishR::predict_mod()} (e.g. \code{curr.E}, \code{curr.Lc},
#'   \code{Lmin}, \code{Lincr}).
#'
#' @return Object of class \code{FishStockYPR}: a list containing the raw
#'   output of \code{predict_mod()} (\code{$mod}), the extracted reference
#'   points (\code{$reference_points}), the current state (\code{$currents})
#'   and the engine type (\code{"ypr"}).
#' @references Beverton & Holt (1957) ; Sparre & Venema (1998), FAO 306/1.
#' @seealso \code{\link{ypr_reference_points}}, \code{\link{run_thompson_bell}}
#' @examples
#' \dontrun{
#' p  <- ypr_param(mortality = mort, allometry = lw, Lc = 12)
#' yp <- run_ypr(p, FM_change = seq(0, 3, 0.05))
#' yp$reference_points
#' plot(yp)
#' }
#' @export
run_ypr <- function(param,
                    FM_change = seq(0, 3, 0.05),
                    Lc_change = NULL,
                    s_list = NULL,
                    FM_relative = FALSE,
                    plot = FALSE,
                    ...) {

  .check_tropfishr()

  ## The YPR engine of predict_mod only accepts absolute F values.
  if (isTRUE(FM_relative)) {
    warning("The YPR engine only accepts absolute F values: ",
            "'FM_relative' forced to FALSE.", call. = FALSE)
    FM_relative <- FALSE
  }

  ## --- Determine Lc_change (length-based) ---
  if (is.null(Lc_change) && !is.null(param$Lc)) Lc_change <- param$Lc

  ## --- In length-based mode, predict_mod requires Lr (recruitment length) ---
  length_based <- !is.null(param$Linf) &&
    is.null(param$tr) && is.null(param$Winf)
  if (length_based && is.null(param$Lr)) {
    param$Lr <- if (!is.null(Lc_change)) min(Lc_change) else
      stop("In length-based mode, 'Lr' (recruitment length) is required in ",
           "param, or 'Lc_change' must be provided.", call. = FALSE)
  }

  ## --- Build a default selectivity ogive if needed ---
  if (is.null(s_list) && !is.null(Lc_change)) {
    Lc1 <- Lc_change[1]
    s_list <- list(selecType = "trawl_ogive", L50 = Lc1, L75 = Lc1 * 1.08)
  }

  args <- c(
    list(param, type = "ypr", FM_change = FM_change,
         FM_relative = FM_relative, plot = plot),
    if (!is.null(Lc_change)) list(Lc_change = Lc_change),
    if (!is.null(s_list))    list(s_list = s_list),
    list(...)
  )

  mod <- tryCatch(
    do.call(TropFishR::predict_mod, args),
    error = function(e)
      stop("predict_mod(type='ypr') : ", conditionMessage(e), call. = FALSE)
  )

  rp <- tryCatch(ypr_reference_points(mod),
                 error = function(e) NULL)

  structure(
    list(mod = mod, reference_points = rp,
         currents = mod$currents, engine = "ypr",
         FM_relative = FM_relative),
    class = "FishStockYPR"
  )

}


#==============================================================
# Extraction of per-recruit reference points
#==============================================================

#' Extract the per-recruit reference points
#'
#' Retrieves, from the output of \code{TropFishR::predict_mod()}, the
#' biological reference points already computed by the official function:
#' \eqn{F_{max}} (maximum yield per recruit), \eqn{F_{0.1}} (precautionary
#' management point: slope = 10 \% of the slope at the origin) and
#' \eqn{F_{0.5}} / \eqn{E_{0.5}} (effort reducing biomass per recruit to 50 \%
#' of the virgin state), with the corresponding exploitation rates
#' \eqn{E = F/Z}.
#'
#' @param mod Output of \code{TropFishR::predict_mod()} (or a
#'   \code{FishStockYPR} object).
#' @return \code{data.frame} of reference points (one row per \code{Lc}/
#'   \code{tc} value explored).
#' @references Gulland (1983) for \eqn{F_{0.1}} ; Sparre & Venema (1998).
#' @export
ypr_reference_points <- function(mod) {

  if (inherits(mod, "FishStockYPR")) mod <- mod$mod

  ## predict_mod() stores the reference points in $df_Es
  df <- mod$df_Es
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0)
    stop("No reference points ($df_Es) in the predict_mod() output.",
         call. = FALSE)

  ## predict_mod may duplicate the row when Lc appears both in param
  ## and in Lc_change; we de-duplicate.
  df <- unique(df)
  tibble::as_tibble(df)

}


#==============================================================
# Comparison of the two official engines
#==============================================================

#' Compare the per-recruit engines (Beverton-Holt vs Thompson-Bell)
#'
#' Runs the two OFFICIAL per-recruit models from \pkg{TropFishR} on the same
#' set of parameters and gathers their reference points into a single table,
#' in keeping with the package's philosophy (\dQuote{implement all the
#' methods, then compare}).
#'
#' @param param \code{ypr_param} object or compatible list.
#' @param FM_change Vector of effort values to explore.
#' @param ... Arguments passed to both wrappers.
#' @return Object of class \code{FishStockYPRcompare}: a list containing both
#'   fits (\code{$ypr}, \code{$thompbell}) and a combined table
#'   (\code{$comparison}).
#' @seealso \code{\link{run_ypr}}, \code{\link{run_thompson_bell}}
#' @export
compare_per_recruit <- function(param,
                                FM_change = seq(0, 3, 0.05),
                                ...) {

  .check_tropfishr()

  yp <- tryCatch(run_ypr(param, FM_change = FM_change, ...),
                 error = function(e) {
                   warning("YPR : ", conditionMessage(e), call. = FALSE)
                   NULL
                 })

  tb <- tryCatch(run_thompson_bell(param, FM_change = FM_change, ...),
                 error = function(e) {
                   warning("Thompson-Bell : ", conditionMessage(e),
                           call. = FALSE)
                   NULL
                 })

  rows <- list()

  if (!is.null(yp) && !is.null(yp$reference_points)) {
    d <- yp$reference_points
    rows$ypr <- data.frame(
      moteur = "Beverton-Holt (YPR)",
      Fmax = .first_num(d, c("Fmax", "F_max")),
      F01  = .first_num(d, c("F01", "F0.1", "F_01")),
      F05  = .first_num(d, c("F05", "F0.5", "F_05")),
      Emax = .first_num(d, c("Emax", "E_max")),
      E01  = .first_num(d, c("E01", "E0.1")),
      E05  = .first_num(d, c("E05", "E0.5")),
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(tb) && !is.null(tb$df_Es)) {
    d <- tb$df_Es
    rows$tb <- data.frame(
      moteur = "Thompson-Bell",
      Fmax = .first_num(d, c("Fmax", "F_max")),
      F01  = .first_num(d, c("F01", "F0.1", "F_01")),
      F05  = .first_num(d, c("F05", "F0.5", "F_05")),
      Emax = .first_num(d, c("Emax", "E_max")),
      E01  = .first_num(d, c("E01", "E0.1")),
      E05  = .first_num(d, c("E05", "E0.5")),
      stringsAsFactors = FALSE
    )
  }

  comparison <- if (length(rows)) do.call(rbind, rows) else NULL

  structure(
    list(ypr = yp, thompbell = tb, comparison = comparison),
    class = "FishStockYPRcompare"
  )

}


## Returns the first numeric column found (1st value) among 'cands'.
.first_num <- function(df, cands) {
  for (nm in cands) if (nm %in% names(df))
    return(suppressWarnings(as.numeric(df[[nm]])[1]))
  NA_real_
}


#==============================================================
# S3 methods
#==============================================================

#' @export
print.FishStockYPR <- function(x, ...) {

  cat("\n=====================================\n")
  cat(" stockflow : per-recruit analysis (Beverton & Holt)\n")
  cat("=====================================\n\n")

  if (!is.null(x$reference_points)) {
    cat("Reference points (predict_mod) :\n")
    print(as.data.frame(x$reference_points), row.names = FALSE)
  } else {
    cat("Reference points unavailable.\n")
  }

  invisible(x)

}

#' @export
summary.FishStockYPR <- function(object, ...) {
  print(object)
  invisible(object$reference_points)
}

#' Yield-per-recruit diagram
#' @param x \code{FishStockYPR} object.
#' @param ... Ignored.
#' @export
plot.FishStockYPR <- function(x, ...) {
  .check_tropfishr()
  ## predict_mod() knows how to plot itself via its object.
  graphics::plot(x$mod)
  invisible(x)
}

#' @export
print.FishStockYPRcompare <- function(x, ...) {

  cat("\n=====================================\n")
  cat(" stockflow : comparison of per-recruit engines\n")
  cat("=====================================\n\n")

  if (!is.null(x$comparison))
    print(x$comparison, row.names = FALSE)
  else
    cat("No comparison available (both engines failed).\n")

  invisible(x)

}