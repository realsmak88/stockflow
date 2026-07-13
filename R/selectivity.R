###############################################################
#
# stockflow
#
# Module: selectivity.R
#
# Technical management measures: gear selectivity,
# mesh size, size at first capture (Lc), minimum landing
# size (MLS) and optimal size.
#
# Objective (iii) of the management plan.
#
# Wrappers around OFFICIAL TropFishR functions only:
#   select_Millar()  -> fitting the selectivity of gillnets
#                       (experimental mesh size data)
#   select_ogive()   -> retention curve (trawl ogive,
#                       knife edge)
#
# The reference biological sizes (L_mat, L_opt, Lc_opt)
# come from modules already present (allometry/maturity,
# lbb_reference_points). No scientific algorithm is rewritten.
#
# References:
#   Millar, R.B. & Holst, R. (1997). Estimation of gillnet and hook
#     selectivity using log-linear models. ICES J. Mar. Sci. 54.
#   Sparre, P. & Venema, S.C. (1998). FAO Fish. Tech. Pap. 306/1.
#   Froese, R. et al. (2016). Minimizing the impact of fishing.
#     Fish and Fisheries 17(3).
#
###############################################################


#==============================================================
# Dependencies
#==============================================================

.check_tropfishr_sel <- function() {
  if (!requireNamespace("TropFishR", quietly = TRUE))
    stop("The 'TropFishR' package must be installed (selectivity).",
         call. = FALSE)
  invisible(TRUE)
}


#==============================================================
# 1. Gillnet selectivity (select_Millar)
#==============================================================

#' Fit the selectivity of a gillnet (Millar, TropFishR)
#'
#' Wrapper around \code{TropFishR::select_Millar()}. Fits a log-linear
#' selectivity model to experimental catch-per-mesh data
#' (size x mesh matrix), and derives the relationship between mesh size and
#' selected size.
#'
#' @param data List in the format expected by \code{select_Millar}:
#'   \code{midLengths} (size classes), \code{meshSizes} (mesh sizes) and
#'   \code{CatchPerNet_mat} (catch-per-net matrix).
#' @param x0 Initial values of the model (see \code{?select_Millar}).
#' @param rtype Type of retention curve: \code{"norm.loc"} (default),
#'   \code{"norm.sca"}, \code{"lognorm"}, \code{"binorm.sca"}, \code{"gamma"}.
#' @param rel.power Relative fishing power per mesh size (optional).
#' @param plot Logical: plot the \code{select_Millar} diagnostic.
#'
#' @return Object of class \code{FishStockSelectivity}: list containing the
#'   raw output of \code{select_Millar()} (\code{$fit}), the gear type
#'   (\code{"gillnet"}) and the main estimates (\code{$estimates}).
#' @references Millar & Holst (1997); Sparre & Venema (1998).
#' @seealso \code{\link{gillnet_mesh_for_Lc}}, \code{\link{retention_curve}}
#' @examples
#' \dontrun{
#' data(gillnet, package = "TropFishR")
#' fit <- fit_gillnet_selectivity(gillnet, x0 = c(60, 4))
#' }
#' @export
fit_gillnet_selectivity <- function(data,
                                    x0 = NULL,
                                    rtype = "norm.loc",
                                    rel.power = NULL,
                                    plot = FALSE) {

  .check_tropfishr_sel()

  fit <- tryCatch(
    TropFishR::select_Millar(data, x0 = x0, rtype = rtype,
                             rel.power = rel.power, plot = plot),
    error = function(e)
      stop("select_Millar : ", conditionMessage(e), call. = FALSE)
  )

  structure(
    list(fit = fit, gear = "gillnet",
         rtype = rtype, estimates = fit$estimates),
    class = "FishStockSelectivity"
  )

}


#' Target mesh size for a given size at first capture
#'
#' Based on a fitted gillnet selectivity model
#' (\code{fit_gillnet_selectivity()}), computes the mesh size that centers the
#' selection curve on a target size \code{Lc}. For models of type
#' \code{norm.loc} / \code{lognorm}, the modal size is proportional to the
#' mesh size (\eqn{mode = k \cdot mesh}), which provides an invertible
#' linear relationship.
#'
#' @param selectivity Object \code{FishStockSelectivity} (gear \code{"gillnet"}).
#' @param Lc Target size(s) at first capture.
#' @param mesh_ref,mode_ref Reference mesh size and modal size defining
#'   the proportionality. By default, derived from the first experimental
#'   mesh size and the estimated modal size.
#'
#' @return \code{data.frame} \code{Lc, mesh} (recommended mesh size per size).
#' @seealso \code{\link{fit_gillnet_selectivity}}
#' @export
gillnet_mesh_for_Lc <- function(selectivity, Lc,
                                mesh_ref = NULL, mode_ref = NULL) {

  if (!inherits(selectivity, "FishStockSelectivity") ||
      selectivity$gear != "gillnet")
    stop("'selectivity' must be a gillnet fit.",
         call. = FALSE)

  fit <- selectivity$fit

  ## Modal size at the first experimental mesh size (mode1) and associated mesh size.
  if (is.null(mode_ref)) mode_ref <- fit$estimates[1, 1]
  if (is.null(mesh_ref)) mesh_ref <- fit$meshSizes[1]

  if (!is.finite(mode_ref) || !is.finite(mesh_ref) || mesh_ref <= 0)
    stop("Invalid mesh/mode reference.", call. = FALSE)

  k <- mode_ref / mesh_ref               # modal size per unit of mesh size
  data.frame(Lc = Lc, mesh = Lc / k)

}


#==============================================================
# 2. Retention curve (select_ogive)
#==============================================================

#' Retention curve of a gear (trawl ogive, TropFishR)
#'
#' Wrapper around \code{TropFishR::select_ogive()}. Builds the proportion
#' retained per size class for a trawl ogive (\code{trawl_ogive},
#' defined by \eqn{L_{50}} and \eqn{L_{75}}) or a knife-edge selectivity
#' (\code{knife_edge}, defined by \eqn{L_c}).
#'
#' @param Lt Vector of sizes at which to evaluate retention.
#' @param type \code{"trawl_ogive"} (default) or \code{"knife_edge"}.
#' @param L50,L75 Sizes at 50 \% and 75 \% retention (\code{trawl_ogive}).
#' @param Lc Cutoff size (\code{knife_edge}).
#'
#' @return \code{data.frame} \code{Lt, retention}.
#' @references Sparre & Venema (1998), FAO 306/1.
#' @seealso \code{\link{fit_gillnet_selectivity}}, \code{\link{recommend_sizes}}
#' @examples
#' \dontrun{
#' retention_curve(Lt = seq(10, 50, 2), L50 = 28, L75 = 32)
#' }
#' @export
retention_curve <- function(Lt,
                            type = c("trawl_ogive", "knife_edge"),
                            L50 = NULL, L75 = NULL, Lc = NULL) {

  .check_tropfishr_sel()
  type <- match.arg(type)

  s_list <- if (type == "trawl_ogive") {
    if (is.null(L50) || is.null(L75))
      stop("'trawl_ogive' requires L50 and L75.", call. = FALSE)
    list(selecType = "trawl_ogive", L50 = L50, L75 = L75)
  } else {
    if (is.null(Lc)) stop("'knife_edge' requires Lc.", call. = FALSE)
    list(selecType = "knife_edge", L50 = Lc)
  }

  ret <- tryCatch(
    TropFishR::select_ogive(s_list, Lt = Lt),
    error = function(e)
      stop("select_ogive : ", conditionMessage(e), call. = FALSE)
  )

  data.frame(Lt = Lt, retention = as.numeric(ret))

}


#==============================================================
# 3. Recommendation of management sizes (Lc, MLS, Lopt)
#==============================================================

#' Recommend management sizes (Lc, MLS, optimal size)
#'
#' Compares one or more candidate sizes at first capture against the
#' reference biological sizes of the stock: size at maturity \eqn{L_{50}},
#' optimal size \eqn{L_{opt}} and optimal size at first capture
#' \eqn{L_{c,opt}} (Froese; \code{lbb_reference_points()}). For each
#' candidate, a transparent verdict indicates whether it protects
#' reproduction (\eqn{L_c \ge L_{50}}) and whether it is close to the
#' optimal size.
#'
#' @param candidates Vector of sizes at first capture to evaluate (cm).
#' @param Linf Asymptotic length (cm).
#' @param L50_maturity Size at 50 \% maturity (cm). Optional but
#'   recommended: without it, the reproduction protection criterion is
#'   not evaluated.
#' @param MK Ratio \eqn{M/K} (default 1.5) for \eqn{L_{opt}} and
#'   \eqn{L_{c,opt}}.
#' @param FM Current ratio \eqn{F/M} (default 1) for \eqn{L_{c,opt}}.
#' @param tolerance Relative tolerance around \eqn{L_{opt}} to judge a
#'   size \dQuote{close to the optimum} (default 0.1, i.e. +/- 10 \%).
#'
#' @return Object of class \code{FishStockSizeAdvice}: list with the
#'   reference sizes (\code{$reference}) and the candidate evaluation table
#'   (\code{$advice}).
#' @references Froese et al. (2016); Beverton (1992) for \eqn{L_{opt}}.
#' @seealso \code{\link{lbb_reference_points}}, \code{\link{retention_curve}}
#' @examples
#' \dontrun{
#' recommend_sizes(candidates = c(10, 12, 14), Linf = 42,
#'                 L50_maturity = 13, MK = 1.4)
#' }
#' @export
recommend_sizes <- function(candidates,
                            Linf,
                            L50_maturity = NULL,
                            MK = 1.5,
                            FM = 1,
                            tolerance = 0.1) {

  if (!is.numeric(candidates) || length(candidates) == 0)
    stop("'candidates' must be a non-empty numeric vector.",
         call. = FALSE)
  if (!is.numeric(Linf) || Linf <= 0)
    stop("'Linf' must be a positive number.", call. = FALSE)

  ## Reference sizes: Lopt and Lc_opt via the package's official function
  ref <- lbb_reference_points(Linf = Linf, MK = MK, FM = FM)
  Lopt   <- ref$Lopt
  Lc_opt <- ref$Lc_opt

  lo <- Lopt * (1 - tolerance)
  hi <- Lopt * (1 + tolerance)

  advice <- data.frame(
    Lc = candidates,
    protege_reproduction = if (!is.null(L50_maturity))
      candidates >= L50_maturity else NA,
    ratio_Lc_L50 = if (!is.null(L50_maturity))
      round(candidates / L50_maturity, 3) else NA_real_,
    proche_Lopt = candidates >= lo & candidates <= hi,
    ratio_Lc_Lopt = round(candidates / Lopt, 3),
    stringsAsFactors = FALSE
  )

  ## Summary verdict
  advice$verdict <- vapply(seq_len(nrow(advice)), function(i) {
    prot <- advice$protege_reproduction[i]
    prox <- advice$proche_Lopt[i]
    if (!is.na(prot) && !prot) return("Insufficient (Lc < L50 maturity)")
    if (isTRUE(prox)) return("Optimal (close to Lopt)")
    if (candidates[i] < lo) return("Acceptable but below the optimum")
    "Above the optimum"
  }, character(1))

  structure(
    list(
      reference = list(Linf = Linf, L50_maturity = L50_maturity,
                       Lopt = Lopt, Lc_opt = Lc_opt, MK = MK, FM = FM),
      advice = advice,
      tolerance = tolerance
    ),
    class = "FishStockSizeAdvice"
  )

}


#==============================================================
# S3 methods
#==============================================================

#' @export
print.FishStockSelectivity <- function(x, ...) {
  cat("\n=====================================\n")
  cat(" stockflow: gear selectivity (", x$gear, ")\n", sep = "")
  cat("=====================================\n\n")
  cat("Curve type:", x$rtype, "\n\n")
  cat("Estimates:\n")
  print(x$estimates)
  invisible(x)
}

#' Gear selectivity diagram
#' @param x Object \code{FishStockSelectivity}.
#' @param ... Ignored.
#' @export
plot.FishStockSelectivity <- function(x, ...) {
  .check_tropfishr_sel()
  graphics::plot(x$fit)
  invisible(x)
}

#' @export
print.FishStockSizeAdvice <- function(x, ...) {
  cat("\n=====================================\n")
  cat(" stockflow: recommended management sizes\n")
  cat("=====================================\n\n")
  r <- x$reference
  cat(sprintf("Linf = %.1f cm | Lopt = %.1f cm | Lc_opt = %.1f cm\n",
              r$Linf, r$Lopt, r$Lc_opt))
  if (!is.null(r$L50_maturity))
    cat(sprintf("L50 maturity = %.1f cm\n", r$L50_maturity))
  cat("\nEvaluation of candidate sizes:\n")
  print(x$advice, row.names = FALSE)
  invisible(x)
}