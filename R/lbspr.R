###############################################################
#
# stockflow
#
# Module : lbspr.R
#
# LBSPR - Length-Based Spawning Potential Ratio
# (Hordyk, Ono, Sainsbury, Loneragan, Prince 2015)
#
# Principle: robust wrapper around the OFFICIAL 'LBSPR' package.
# The estimator (LBSPRfit) is never rewritten; we handle the
# construction of S4 objects (LB_pars, LB_lengths), the call and
# the extraction of indicators (SPR, SL50, SL95, F/M) + the
# diagnostic relative to management targets.
#
###############################################################

#==============================================================
# Dependencies
#==============================================================

.check_lbspr <- function() {

  if (!requireNamespace("LBSPR", quietly = TRUE))
    stop("The 'LBSPR' package must be installed.", call. = FALSE)

  if (!requireNamespace("methods", quietly = TRUE))
    stop("The 'methods' package is required (S4 objects).", call. = FALSE)

  invisible(TRUE)

}


#==============================================================
# Construction of biological parameters (LB_pars)
#==============================================================

#' Construction of an LB_pars object (LBSPR)
#'
#' Assembles the biological parameters needed for LBSPR into an S4 object
#' \code{LB_pars}. The ratio \eqn{M/K} can be provided directly
#' (\code{MK}) or derived from \code{M} and \code{K}.
#'
#' @param Linf Asymptotic length.
#' @param L50,L95 Lengths at 50 and 95 percent maturity.
#' @param MK Ratio \eqn{M/K}. Takes priority if provided.
#' @param M,K Natural mortality and growth coefficient (used if
#'   \code{MK} is \code{NULL}).
#' @param species Species name (label).
#' @param L_units Length unit (\code{"cm"} by default).
#' @param bin_width Width of the size classes.
#' @param Walpha,Wbeta Length-weight coefficients (optional; \code{Wbeta}
#'   also serves as the default fecundity exponent).
#'
#' @return S4 object \code{LB_pars} (LBSPR package).
#' @references Hordyk, A. et al. (2015). ICES J. Mar. Sci. 72(1): 217-231.
#' @export
lbspr_pars <- function(Linf,
                       L50,
                       L95,
                       MK = NULL,
                       M = NULL,
                       K = NULL,
                       species = "stock",
                       L_units = "cm",
                       bin_width = NULL,
                       Walpha = NULL,
                       Wbeta = NULL) {

  .check_lbspr()

  if (is.null(MK)) {

    if (is.null(M) || is.null(K))
      stop("Provide 'MK', or 'M' and 'K'.")

    MK <- M / K

  }

  if (L95 <= L50)
    stop("L95 must be greater than L50.")


  pars <- methods::new("LB_pars", verbose = FALSE)

  pars@Species <- species
  pars@Linf    <- as.numeric(Linf)
  pars@L50     <- as.numeric(L50)
  pars@L95     <- as.numeric(L95)
  pars@MK      <- as.numeric(MK)
  pars@L_units <- L_units

  if (!is.null(bin_width)) pars@BinWidth <- bin_width
  if (!is.null(Walpha))    pars@Walpha   <- Walpha
  if (!is.null(Wbeta)) {
    pars@Wbeta <- Wbeta
    pars@FecB  <- Wbeta
  }

  pars

}


#==============================================================
# Construction of size compositions (LB_lengths)
#==============================================================

#' Construction of an LB_lengths object (LBSPR)
#'
#' Creates an S4 object \code{LB_lengths} either from a frequency matrix
#' (size classes x years), or from individual lengths
#' grouped into classes.
#'
#' @param x Either a matrix/data.frame of frequencies (rows = classes,
#'   columns = years), or a data.frame of individual lengths.
#' @param pars Object \code{LB_pars} (provides \code{BinWidth}, \code{Linf}).
#' @param length_col,year_col Columns (individual lengths mode).
#' @param LMids Class midpoints (matrix mode). If \code{NULL},
#'   derived from row names.
#' @param years Vector of years (matrix mode).
#' @param bin_width Class width (individual lengths mode); by
#'   default that of \code{pars} or 1.
#'
#' @return S4 object \code{LB_lengths}.
#' @export
lbspr_lengths <- function(x,
                          pars,
                          length_col = NULL,
                          year_col = NULL,
                          LMids = NULL,
                          years = NULL,
                          bin_width = NULL) {

  .check_lbspr()

  if (is.null(bin_width))
    bin_width <- if (length(pars@BinWidth)) pars@BinWidth else 1


  ## ---- individual lengths mode ------------------------------
  if (!is.null(length_col) && !is.null(year_col)) {

    if (!all(c(length_col, year_col) %in% names(x)))
      stop("Length/year columns not found.")

    L <- suppressWarnings(as.numeric(x[[length_col]]))
    Y <- x[[year_col]]

    ok <- is.finite(L) & L > 0 & !is.na(Y)
    L <- L[ok]; Y <- Y[ok]

    ## Classes must exceed Linf: LBSPRfit rejects a last class
    ## smaller than the asymptotic size ("Maximum length bin can't be
    ## smaller than asymptotic size"). The OBSERVED maximum may well be
    ## below Linf (sampling, selectivity): so we extend up to
    ## 1.15 x Linf.
    Lmax <- max(L, na.rm = TRUE)

    if (length(pars@Linf) && is.finite(pars@Linf))
      Lmax <- max(Lmax, 1.15 * pars@Linf)

    breaks <- seq(0, ceiling(Lmax / bin_width) * bin_width + bin_width,
                  by = bin_width)
    mids   <- breaks[-length(breaks)] + bin_width / 2

    years  <- sort(unique(Y))

    mat <- vapply(years, function(yr) {
      as.numeric(table(cut(L[Y == yr], breaks = breaks,
                           include.lowest = TRUE)))
    }, numeric(length(mids)))

    LMids <- mids

  } else {

    ## ---- frequency matrix mode ------------------------------
    mat <- as.matrix(x)

    if (is.null(LMids)) {

      rn <- suppressWarnings(as.numeric(rownames(mat)))

      if (all(is.finite(rn)))
        LMids <- rn
      else
        stop("Provide 'LMids' (class midpoints).")

    }

    if (is.null(years))
      years <- if (!is.null(colnames(mat)))
        suppressWarnings(as.numeric(colnames(mat)))
      else seq_len(ncol(mat))

  }


  lbl <- methods::new("LB_lengths")
  lbl@LMids  <- as.numeric(LMids)
  lbl@LData  <- matrix(as.numeric(mat), nrow = length(LMids))
  lbl@Years  <- as.numeric(years)
  lbl@NYears <- length(years)

  lbl

}


#==============================================================
# LBSPR fitting
#==============================================================

#' LBSPR fitting (SPR, selectivity, F/M)
#'
#' Wrapper around \code{LBSPR::LBSPRfit()}. Estimates by year the SPR, the
#' selectivity (SL50, SL95) and the ratio \eqn{F/M} from the size
#' compositions and the biological parameters.
#'
#' @param pars Object \code{LB_pars} (see \code{lbspr_pars()}).
#' @param lengths Object \code{LB_lengths} (see \code{lbspr_lengths()}).
#' @param spr_target,spr_limit Management targets for SPR (defaults 0.40 and
#'   0.20).
#' @param ... Arguments passed to \code{LBSPR::LBSPRfit()}.
#'
#' @return S3 object \code{FishStockLBSPR}: \code{fit} (LB_obj), \code{summary}
#'   (yearly data.frame: Year, SL50, SL95, FM, SPR, status) and
#'   \code{targets}.
#'
#' @references Hordyk, A. et al. (2015). ICES J. Mar. Sci. 72(1).
#' @examples
#' \dontrun{
#'   pr <- lbspr_pars(Linf = 45, L50 = 22, L95 = 28, M = 1.5, K = 1.2)
#'   ln <- lbspr_lengths(freq, pr, length_col = "LCT", year_col = "annee")
#'   fit <- run_lbspr(pr, ln)
#'   fit$summary
#' }
#' @export
run_lbspr <- function(pars,
                      lengths,
                      spr_target = 0.40,
                      spr_limit = 0.20,
                      ...) {

  .check_lbspr()

  if (!methods::is(pars, "LB_pars"))
    stop("'pars' must be an LB_pars object (see lbspr_pars()).")

  if (!methods::is(lengths, "LB_lengths"))
    stop("'lengths' must be an LB_lengths object (see lbspr_lengths()).")


  fit <- tryCatch(

    LBSPR::LBSPRfit(pars, lengths, ...),

    error = function(e)
      stop("LBSPRfit: ", conditionMessage(e), call. = FALSE)

  )


  smy <- extract_lbspr(fit, spr_target, spr_limit)


  structure(

    list(
      fit     = fit,
      summary = smy,
      targets = c(SPR_target = spr_target, SPR_limit = spr_limit)
    ),

    class = "FishStockLBSPR"

  )

}


#==============================================================
# Extraction of LBSPR estimates
#==============================================================

#' Extraction of yearly LBSPR estimates
#'
#' @param fit Object \code{LB_obj} returned by \code{LBSPR::LBSPRfit()}.
#' @param spr_target,spr_limit Targets for the status diagnostic.
#' @return data.frame: Year, SL50, SL95, FM, SPR, status.
#' @export
extract_lbspr <- function(fit, spr_target = 0.40, spr_limit = 0.20) {

  ## @Ests: smoothed estimates (SL50, SL95, FM, SPR); otherwise raw slots
  ests <- tryCatch(methods::slot(fit, "Ests"), error = function(e) NULL)

  years <- tryCatch(methods::slot(fit, "Years"), error = function(e) NULL)

  if (!is.null(ests) && is.matrix(ests) && ncol(ests) >= 4) {

    df <- data.frame(
      Year = if (!is.null(years)) years else seq_len(nrow(ests)),
      SL50 = ests[, "SL50"],
      SL95 = ests[, "SL95"],
      FM   = ests[, "FM"],
      SPR  = ests[, "SPR"],
      stringsAsFactors = FALSE
    )

  } else {

    ## fallback: raw vector slots
    getslot <- function(nm)
      tryCatch(as.numeric(methods::slot(fit, nm)),
               error = function(e) NA_real_)

    df <- data.frame(
      Year = if (!is.null(years)) years else NA_real_,
      SL50 = getslot("SL50"),
      SL95 = getslot("SL95"),
      FM   = getslot("FM"),
      SPR  = getslot("SPR"),
      stringsAsFactors = FALSE
    )

  }


  df$statut <- vapply(df$SPR,
                      .lbspr_status,
                      character(1),
                      target = spr_target,
                      limit = spr_limit)

  rownames(df) <- NULL

  df

}


## Status classification based on SPR (testable helper)
.lbspr_status <- function(spr, target = 0.40, limit = 0.20) {

  if (is.na(spr)) return(NA_character_)

  if (spr >= target) "Healthy (SPR >= target)"
  else if (spr >= limit) "Concerning (limit <= SPR < target)"
  else "Overexploited (SPR < limit)"

}


#==============================================================
# Equilibrium simulation (SPR ~ F/M or ~ Lc curve)
#==============================================================

#' Equilibrium LBSPR simulation
#'
#' Wrapper for \code{LBSPR::LBSPRsim()}: SPR, relative yield and expected
#' size composition for a set of parameters (useful for exploring the effect
#' of a length at first capture or a level of F/M).
#'
#' @param pars Complete \code{LB_pars} object (with \code{SL50}, \code{SL95},
#'   \code{FM} or \code{SPR} filled in).
#' @param ... Arguments passed to \code{LBSPR::LBSPRsim()}.
#' @return Simulation object \code{LB_obj}.
#' @export
lbspr_sim <- function(pars, ...) {

  .check_lbspr()

  tryCatch(
    LBSPR::LBSPRsim(pars, ...),
    error = function(e)
      stop("LBSPRsim: ", conditionMessage(e), call. = FALSE)
  )

}


#==============================================================
# S3 methods
#==============================================================

#' @export
print.FishStockLBSPR <- function(x, ...) {

  cat("\n=====================================\n")
  cat(" stockflow : LBSPR\n")
  cat("=====================================\n\n")

  cat(sprintf("Targets: SPR >= %.2f (target), >= %.2f (limit)\n\n",
              x$targets["SPR_target"], x$targets["SPR_limit"]))

  print(utils::tail(x$summary, 6), row.names = FALSE)

  invisible(x)

}


#' @export
summary.FishStockLBSPR <- function(object, ...) {

  print(object)
  invisible(object$summary)

}


#' Plot of LBSPR estimates (SPR, SL50/SL95, F/M over time)
#' @param x Object \code{FishStockLBSPR}.
#' @param ... Ignored.
#' @export
plot.FishStockLBSPR <- function(x, ...) {

  ## priority given to the official LBSPR plot if a single multi-year
  done <- tryCatch({
    LBSPR::plotEsts(x$fit)
    TRUE
  }, error = function(e) FALSE)

  if (done) return(invisible(x))


  if (!requireNamespace("ggplot2", quietly = TRUE)) {

    with(x$summary, {
      graphics::plot(Year, SPR, type = "b", ylim = c(0, 1),
                     ylab = "SPR", xlab = "year")
      graphics::abline(h = x$targets["SPR_target"], lty = 2, col = "darkgreen")
      graphics::abline(h = x$targets["SPR_limit"],  lty = 3, col = "red")
    })

    return(invisible(x))

  }


  df <- x$summary

  p <- ggplot2::ggplot(df, ggplot2::aes(Year, SPR)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0,
                                      ymax = x$targets["SPR_limit"]),
                         fill = "#FF7F7F", alpha = .2) +
    ggplot2::geom_line(colour = "#12507b") +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_hline(yintercept = x$targets["SPR_target"],
                        linetype = 2, colour = "darkgreen") +
    ggplot2::geom_hline(yintercept = x$targets["SPR_limit"],
                        linetype = 3, colour = "red") +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(title = "LBSPR: spawning potential ratio (SPR)",
                  x = NULL, y = "SPR") +
    ggplot2::theme_bw()

  print(p)

  invisible(p)

}