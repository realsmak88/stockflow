###############################################################
#
# stockflow
#
# Module : allometry.R
#
# Length-weight relationship (allometry) W = a * L^b
# Wrappers around stats::lm() (log-log) and the OFFICIAL
# functions of FSA (FSA::hoCoef for the isometry test).
#
###############################################################

#==============================================================
# Fitting the length-weight relationship
#==============================================================

#' Length-weight relationship (allometry) W = a * L^b
#'
#' Fits the length-weight relationship by log-log linear regression and tests
#' isometry (\eqn{b = 3}) via the official function \code{FSA::hoCoef()}.
#' An inter-group comparison (ANCOVA) is provided when \code{group}
#' is specified.
#'
#' @param data data.frame containing the lengths and weights.
#' @param length_col,weight_col Names of the length and weight columns.
#' @param group (optional) Name of a group column (sex, species, zone)
#'   to compare relationships.
#' @param log_base Logarithmic base (\code{10} by default, fisheries
#'   convention; \code{exp(1)} possible).
#'
#' @return S3 object \code{FishStockLW}: list with \code{a}, \code{b}, the CI of
#'   \code{b}, the isometry test, the \eqn{r^2}, the sample size, the
#'   \code{lm} model and, where applicable, the group comparison.
#'
#' @references
#'   Le Cren, E.D. (1951). J. Anim. Ecol. 20(2).
#'   Froese, R. (2006). J. Appl. Ichthyol. 22.
#'   Ogle, D.H. (2016). Introductory Fisheries Analyses with R.
#'
#' @examples
#' \dontrun{
#'   lw <- fit_length_weight(bio, "LM", "weight")
#'   lw$a; lw$b
#' }
#' @export
fit_length_weight <- function(data,
                              length_col,
                              weight_col,
                              group = NULL,
                              log_base = 10) {

  if (!all(c(length_col, weight_col) %in% names(data)))
    stop("Length/weight columns not found in 'data'.")


  d <- data.frame(
    L = suppressWarnings(as.numeric(data[[length_col]])),
    W = suppressWarnings(as.numeric(data[[weight_col]]))
  )

  if (!is.null(group))
    d$grp <- as.factor(data[[group]])

  d <- d[is.finite(d$L) & is.finite(d$W) & d$L > 0 & d$W > 0, ]

  if (nrow(d) < 10)
    stop("Insufficient sample size (<10) to fit the length-weight relationship.")


  logf <- function(x) log(x, base = log_base)

  d$logL <- logf(d$L)
  d$logW <- logf(d$W)


  ## ---- main fit (all groups combined) --------------
  mdl <- stats::lm(logW ~ logL, data = d)

  co  <- stats::coef(mdl)
  b   <- unname(co[2])
  a   <- log_base ^ unname(co[1])
  r2  <- summary(mdl)$r.squared
  ci  <- stats::confint(mdl)["logL", ]


  ## ---- isometry test b = 3 via FSA (official) ------------------
  iso <- NULL

  if (requireNamespace("FSA", quietly = TRUE)) {

    iso <- tryCatch(
      FSA::hoCoef(mdl, term = 2, bo = 3),
      error = function(e) NULL
    )

  }

  ## fallback: manual t test if FSA unavailable
  if (is.null(iso)) {

    se <- summary(mdl)$coefficients["logL", "Std. Error"]
    tval <- (b - 3) / se
    pval <- 2 * stats::pt(-abs(tval), df = stats::df.residual(mdl))
    iso <- data.frame(Ho.value = 3, b = b, t = tval, p.value = pval)

  }

  allometry_type <-
    if (iso$p.value[1] >= 0.05) "isometric (b = 3)"
  else if (b > 3) "positive allometry (b > 3)"
  else "negative allometry (b < 3)"


  ## ---- group comparison (ANCOVA) ----------------------------
  group_test <- NULL

  if (!is.null(group) && nlevels(d$grp) > 1) {

    m_full <- stats::lm(logW ~ logL * grp, data = d)
    m_add  <- stats::lm(logW ~ logL + grp, data = d)
    m_base <- stats::lm(logW ~ logL,       data = d)

    group_test <- list(
      pentes     = stats::anova(m_add, m_full),   # difference in b
      ordonnees  = stats::anova(m_base, m_add)     # difference in a
    )

  }


  structure(

    list(
      a            = a,
      b            = b,
      b_ci         = as.numeric(ci),
      r2           = r2,
      n            = nrow(d),
      isometry     = iso,
      allometry    = allometry_type,
      model        = mdl,
      group_test   = group_test,
      log_base     = log_base
    ),

    class = "FishStockLW"

  )

}


#' @export
print.FishStockLW <- function(x, ...) {

  cat("\n")
  cat("--------------------------------------\n")
  cat("Length-weight relationship  W = a * L^b\n")
  cat("--------------------------------------\n\n")

  cat(sprintf("a         : %.5g\n", x$a))
  cat(sprintf("b         : %.4f\n", x$b))
  cat(sprintf("95%% CI (b): [%.3f ; %.3f]\n", x$b_ci[1], x$b_ci[2]))
  cat(sprintf("r2        : %.4f\n", x$r2))
  cat(sprintf("n         : %d\n", x$n))
  cat(sprintf("Allometry: %s\n", x$allometry))

  if (!is.null(x$group_test))
    cat("\n(Group comparison available in $group_test)\n")

  invisible(x)

}


#==============================================================
# Weight predictor from length
#==============================================================

#' weight predicted from length (from a FishStockLW object)
#'
#' @param object \code{FishStockLW} object.
#' @param length Vector of lengths.
#' @param ... Ignored.
#' @return Vector of predicted weights.
#' @export
predict.FishStockLW <- function(object, length, ...) {

  object$a * length ^ object$b

}