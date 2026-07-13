###############################################################
#
# stockflow
#
# Module : tac.R
#
# Determination of the Total Allowable Catch (TAC) from
# biological reference points, and allocation of quotas by
# segment / fleet in proportion to historical catches.
#
# Objective (i) of the management plan :
#   reference points -> TAC -> allocation by segment.
#
# This module does NOT invent any assessment algorithm: it applies
# standard, transparent control rules (harvest control rules) to
# the outputs of the official modules (mortality, YPR,
# LBSPR, Thompson-Bell) already produced by the pipeline.
#
# References :
#   Gulland, J.A. (1971). The fish resources of the ocean.
#   Restrepo, V.R. et al. (1998). NOAA Tech. Memo. NMFS-F/SPO-31
#     (harvest control rules, reference F/M and E).
#   FAO (1997). Fisheries management. FAO Tech. Guidelines 4.
#
###############################################################


#==============================================================
# 1. TAC calculation
#==============================================================

#' Calculate a Total Allowable Catch (TAC)
#'
#' Translates an exploitation reference point into an authorized catch, using
#' a transparent control rule. Three complementary approaches are
#' offered depending on the data available:
#'
#' \describe{
#'   \item{\code{"ratio_E"}}{Rule based on the exploitation rate. The TAC is the
#'     recent catch scaled by the ratio between the
#'     target exploitation rate (\eqn{E_{cible}}, e.g. \eqn{E_{0.1}} or a target
#'     derived from LBSPR) and the current rate \eqn{E_{cur}} :
#'     \eqn{TAC = C_{cur} \cdot E_{cible} / E_{cur}}. Suitable for data-limited
#'     stocks where only \eqn{E}, \eqn{F/M} and catch are known.}
#'   \item{\code{"ratio_F"}}{Same but based on fishing mortality :
#'     \eqn{TAC = C_{cur} \cdot F_{cible} / F_{cur}}.}
#'   \item{\code{"biomass"}}{Rule based on biomass : \eqn{TAC = F_{cible} \cdot B}
#'     (with \eqn{B} the estimated exploitable biomass). Use only if an
#'     absolute biomass estimate is available.}
#' }
#'
#' @param method Calculation approach : \code{"ratio_E"} (default), \code{"ratio_F"}
#'   or \code{"biomass"}.
#' @param catch_current Recent reference catch \eqn{C_{cur}} (t), for the
#'   ratio-based methods.
#' @param E_current,E_target Current and target exploitation rates (\code{ratio_E}).
#' @param F_current,F_target Current and target fishing mortalities
#'   (\code{ratio_F} ; \code{F_target} also used by \code{biomass}).
#' @param biomass Exploitable biomass \eqn{B} (t), for \code{biomass}.
#' @param buffer Precautionary fraction applied to the TAC (0-1 ; the TAC is
#'   multiplied by \code{1 - buffer}). Default 0 (no reduction).
#' @param label Stock/scenario label.
#'
#' @return Object of class \code{FishStockTAC} : list containing the raw TAC, the
#'   TAC after precaution (\code{$TAC}), the method and the inputs.
#' @references Restrepo et al. (1998) ; FAO (1997) Tech. Guidelines 4.
#' @examples
#' \dontrun{
#' # From the pipeline outputs (current E = mortality$summary,
#' # target E = YPR reference point)
#' compute_tac(method = "ratio_E", catch_current = 5000,
#'             E_current = 0.62, E_target = 0.40, buffer = 0.10)
#' }
#' @export
compute_tac <- function(method = c("ratio_E", "ratio_F", "biomass"),
                        catch_current = NULL,
                        E_current = NULL, E_target = NULL,
                        F_current = NULL, F_target = NULL,
                        biomass = NULL,
                        buffer = 0,
                        label = "stock") {

  method <- match.arg(method)

  if (!is.numeric(buffer) || buffer < 0 || buffer >= 1)
    stop("'buffer' must be in [0, 1).", call. = FALSE)

  tac_raw <- switch(
    method,

    ratio_E = {
      if (is.null(catch_current) || is.null(E_current) || is.null(E_target))
        stop("method='ratio_E' requires catch_current, E_current, E_target.",
             call. = FALSE)
      if (E_current <= 0) stop("E_current must be > 0.", call. = FALSE)
      catch_current * (E_target / E_current)
    },

    ratio_F = {
      if (is.null(catch_current) || is.null(F_current) || is.null(F_target))
        stop("method='ratio_F' requires catch_current, F_current, F_target.",
             call. = FALSE)
      if (F_current <= 0) stop("F_current must be > 0.", call. = FALSE)
      catch_current * (F_target / F_current)
    },

    biomass = {
      if (is.null(biomass) || is.null(F_target))
        stop("method='biomass' requires biomass and F_target.", call. = FALSE)
      F_target * biomass
    }
  )

  tac <- tac_raw * (1 - buffer)

  structure(
    list(
      label   = label,
      method  = method,
      TAC_raw = tac_raw,
      buffer  = buffer,
      TAC     = tac,
      inputs  = list(catch_current = catch_current,
                     E_current = E_current, E_target = E_target,
                     F_current = F_current, F_target = F_target,
                     biomass = biomass)
    ),
    class = "FishStockTAC"
  )

}


#' @export
print.FishStockTAC <- function(x, ...) {

  cat("\n=====================================\n")
  cat(" stockflow : Total Allowable Catch (TAC)\n")
  cat("=====================================\n\n")
  cat(sprintf("Stock    : %s\n", x$label))
  cat(sprintf("Method   : %s\n", x$method))
  cat(sprintf("Raw TAC  : %.1f t\n", x$TAC_raw))
  if (x$buffer > 0)
    cat(sprintf("Precaution : -%.0f%%\n", 100 * x$buffer))
  cat(sprintf("Final TAC: %.1f t\n", x$TAC))
  invisible(x)

}


#==============================================================
# 2. Historical catch proportions
#==============================================================

#' Calculate proportions from historical catches
#'
#' Aggregates a table of historical catches (long format) by segment
#' (\code{espece}, \code{engin}, \code{armement}, \code{port}...) over a range
#' of years, and returns the relative share of each segment. Serves as an
#' allocation key for \code{allocate_quota()}.
#'
#' @param catches \code{data.frame} of catches in long format, with at least
#'   a segment column and a catch column.
#' @param segment_col Name of the segment column (default \code{"espece"}).
#' @param catch_col Name of the catch column (default \code{"capture_t"}).
#' @param year_col Name of the year column (optional ; required if
#'   \code{years} is provided).
#' @param years Vector of years to keep (optional ; otherwise all).
#'
#' @return \code{data.frame} \code{segment, capture_totale, prorata} sorted by
#'   decreasing proportion (\code{prorata} sums to 1).
#' @examples
#' \dontrun{
#' compute_prorata(captures, segment_col = "espece",
#'                 catch_col = "capture_t", year_col = "annee",
#'                 years = 2017:2023)
#' }
#' @export
compute_prorata <- function(catches,
                            segment_col = "espece",
                            catch_col = "capture_t",
                            year_col = NULL,
                            years = NULL) {

  if (!is.data.frame(catches))
    stop("'catches' must be a data.frame.", call. = FALSE)
  if (!segment_col %in% names(catches))
    stop("Segment column '", segment_col, "' missing.", call. = FALSE)
  if (!catch_col %in% names(catches))
    stop("Catch column '", catch_col, "' missing.", call. = FALSE)

  d <- catches

  if (!is.null(years)) {
    if (is.null(year_col) || !year_col %in% names(d))
      stop("'years' provided but 'year_col' missing or not found.",
           call. = FALSE)
    d <- d[d[[year_col]] %in% years, , drop = FALSE]
    if (nrow(d) == 0)
      stop("No rows for the requested years.", call. = FALSE)
  }

  seg <- as.character(d[[segment_col]])
  cap <- suppressWarnings(as.numeric(d[[catch_col]]))
  ok  <- !is.na(seg) & !is.na(cap)
  seg <- seg[ok]; cap <- cap[ok]

  agg <- tapply(cap, seg, sum, na.rm = TRUE)
  total <- sum(agg)
  if (total <= 0) stop("Sum of catches is zero or negative.", call. = FALSE)

  out <- data.frame(
    segment        = names(agg),
    capture_totale = as.numeric(agg),
    prorata        = as.numeric(agg) / total,
    stringsAsFactors = FALSE
  )
  out[order(-out$prorata), , drop = FALSE]

}


#==============================================================
# 3. Allocation of the TAC into quotas
#==============================================================

#' Allocate a TAC into quotas by segment / fleet
#'
#' Distributes a Total Allowable Catch among segments (species, gear,
#' fleets, ports...) according to an allocation key. The key can be:
#' a named vector of proportions, a \code{data.frame} of proportions (produced
#' by \code{compute_prorata()}), or a \code{data.frame} of historical catches
#' aggregated on the fly.
#'
#' @param tac Total TAC to allocate : a number (t) or a
#'   \code{FishStockTAC} object.
#' @param prorata Allocation key. Either a named numeric vector (the names
#'   are the segments), or a \code{data.frame} containing the columns
#'   \code{segment} and \code{prorata} (see \code{compute_prorata()}).
#' @param min_quota Optional floor quota per segment (t) ; segments
#'   below this threshold are raised and the excess is recovered proportionally
#'   from the others (conservative distribution of the total).
#' @param round_digits Number of decimal places for rounding the quotas (default 1).
#'
#' @return Object of class \code{FishStockQuota} : list with the TAC, the
#'   \code{$quotas} table (\code{segment, prorata, quota_t}) and a sum check.
#' @seealso \code{\link{compute_tac}}, \code{\link{compute_prorata}}
#' @examples
#' \dontrun{
#' pr <- compute_prorata(captures, years = 2017:2023)
#' allocate_quota(tac = 4500, prorata = pr)
#' }
#' @export
allocate_quota <- function(tac,
                           prorata,
                           min_quota = NULL,
                           round_digits = 1) {

  tac_val <- if (inherits(tac, "FishStockTAC")) tac$TAC else as.numeric(tac)[1]
  if (!is.finite(tac_val) || tac_val <= 0)
    stop("Invalid TAC (must be a positive number).", call. = FALSE)

  ## --- Normalize 'prorata' into (segment, prorata) ---
  if (is.data.frame(prorata)) {
    if (!all(c("segment", "prorata") %in% names(prorata)))
      stop("The 'prorata' data.frame must contain 'segment' and 'prorata'.",
           call. = FALSE)
    seg <- as.character(prorata$segment)
    p   <- as.numeric(prorata$prorata)
  } else if (is.numeric(prorata) && !is.null(names(prorata))) {
    seg <- names(prorata)
    p   <- as.numeric(prorata)
  } else {
    stop("'prorata' must be a named numeric vector or a data.frame ",
         "(segment, prorata).", call. = FALSE)
  }

  if (any(!is.finite(p)) || any(p < 0))
    stop("Invalid proportions (non-finite or negative values).",
         call. = FALSE)

  ## Defensive renormalization (the proportions must sum to 1).
  s <- sum(p)
  if (s <= 0) stop("Sum of proportions is zero.", call. = FALSE)
  if (abs(s - 1) > 1e-6)
    warning("The proportions do not sum to 1 (sum = ", round(s, 4),
            ") : renormalization applied.", call. = FALSE)
  p <- p / s

  quota <- tac_val * p

  ## --- Optional floor with sum conservation ---
  if (!is.null(min_quota) && is.finite(min_quota) && min_quota > 0) {
    below <- quota < min_quota
    if (any(below) && !all(below)) {
      deficit <- sum(min_quota - quota[below])
      quota[below] <- min_quota
      ## recover the deficit proportionally from segments above the floor
      donors <- !below
      w <- quota[donors] / sum(quota[donors])
      quota[donors] <- quota[donors] - deficit * w
    }
  }

  quotas <- data.frame(
    segment = seg,
    prorata = round(p, 4),
    quota_t = round(quota, round_digits),
    stringsAsFactors = FALSE
  )
  quotas <- quotas[order(-quotas$quota_t), , drop = FALSE]

  structure(
    list(
      TAC        = tac_val,
      quotas     = quotas,
      somme_quotas = sum(quotas$quota_t),
      ecart_somme  = sum(quotas$quota_t) - tac_val
    ),
    class = "FishStockQuota"
  )

}


#' @export
print.FishStockQuota <- function(x, ...) {

  cat("\n=====================================\n")
  cat(" stockflow : Quota allocation\n")
  cat("=====================================\n\n")
  cat(sprintf("Total TAC       : %.1f t\n", x$TAC))
  cat(sprintf("Sum of quotas   : %.1f t (difference %.2f t)\n\n",
              x$somme_quotas, x$ecart_somme))
  print(x$quotas, row.names = FALSE)
  invisible(x)

}