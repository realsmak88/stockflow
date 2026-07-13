###############################################################
#
# stockflow
#
# Module : closure.R
#
# Determination of the ideal window and duration for biological rest
# (spatio-temporal closure).
#
# Objective (ii) of the management plan.
#
# Principle: three monthly signals locate the period to protect,
# each derived from OFFICIAL FUNCTIONS or transparent
# aggregations; no scientific algorithm is rewritten.
#
#   (1) spawning season: monthly proportion of mature individuals
#       (the user declares the maturity codes considered
#        as "mature").
#   (2) Recruitment pulse: TropFishR::recruitment() reconstructs
#       the monthly recruitment pattern from length
#       frequencies and growth (VBGF).
#   (3) Monthly closure score: weighted and normalized combination
#       of the two signals; the optimal window of duration d
#       is the sliding window (circular, since the year is cyclical)
#       that maximizes the cumulative protected score.
#
# References :
#   Sparre, P. & Venema, S.C. (1998). Introduction to tropical fish
#     stock assessment. FAO Fish. Tech. Pap. 306/1 (recruitment).
#   Cochrane, K.L. (ed.) (2002). A fishery manager's guidebook. FAO.
#
###############################################################


#==============================================================
# Internal utilities
#==============================================================

.month_labels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

## Normalizes a vector of 12 values to sum 1 (seasonality pattern).
## Returns a uniform vector if the sum is zero.
.norm_share <- function(x) {
  x[!is.finite(x)] <- 0
  s <- sum(x)
  if (s <= 0) rep(1 / length(x), length(x)) else x / s
}

## Converts a heterogeneous date column into a month number (1-12).
.to_month <- function(x) {
  if (inherits(x, "Date") || inherits(x, "POSIXt"))
    return(as.integer(format(x, "%m")))
  ## attempt string parsing
  d <- suppressWarnings(as.Date(x))
  if (all(is.na(d))) {
    ## maybe already a month number
    m <- suppressWarnings(as.integer(x))
    if (all(is.na(m) | (m >= 1 & m <= 12), na.rm = TRUE)) return(m)
    stop("Unable to extract the month from the date column.",
         call. = FALSE)
  }
  as.integer(format(d, "%m"))
}


#==============================================================
# 1. Spawning season (monthly proportion of mature individuals)
#==============================================================

#' Monthly proportion of mature individuals (spawning signal)
#'
#' Computes, month by month, the fraction of individuals at the mature stage in a
#' set of individual biological data. The peak of this proportion signals the
#' breeding season, a priority target for biological rest.
#'
#' @param bio \code{data.frame} of individual biological data.
#' @param date_col Name of the date (or month) column.
#' @param maturity_col Name of the maturity stage column.
#' @param mature_codes Vector of values of \code{maturity_col} considered
#'   as mature (e.g. \code{c("M", "P")} for immature/mature/spawning, or
#'   \code{c(2, 3, 4)} for a 1-4 scale).
#'
#' @return \code{data.frame} \code{month, month_label, n, n_mature, prop_mature}
#'   (12 rows, missing months filled with 0).
#' @examples
#' \dontrun{
#' spawning_season(bio, date_col = "Date", maturity_col = "maturity",
#'                 mature_codes = c(2, 3, 4))
#' }
#' @export
spawning_season <- function(bio,
                            date_col = "Date",
                            maturity_col = "maturity",
                            mature_codes) {

  if (!is.data.frame(bio))
    stop("'bio' must be a data.frame.", call. = FALSE)
  if (!date_col %in% names(bio))
    stop("Date column '", date_col, "' missing.", call. = FALSE)
  if (!maturity_col %in% names(bio))
    stop("Maturity column '", maturity_col, "' missing.", call. = FALSE)
  if (missing(mature_codes) || length(mature_codes) == 0)
    stop("'mature_codes' is required (maturity codes considered mature).",
         call. = FALSE)

  m   <- .to_month(bio[[date_col]])
  mat <- as.character(bio[[maturity_col]]) %in% as.character(mature_codes)
  ok  <- !is.na(m)
  m   <- m[ok]; mat <- mat[ok]

  n_tot <- tapply(rep(1L, length(m)), factor(m, levels = 1:12), sum)
  n_mat <- tapply(as.integer(mat),    factor(m, levels = 1:12), sum)
  n_tot[is.na(n_tot)] <- 0
  n_mat[is.na(n_mat)] <- 0

  data.frame(
    month       = 1:12,
    month_label = .month_labels,
    n           = as.integer(n_tot),
    n_mature    = as.integer(n_mat),
    prop_mature = ifelse(n_tot > 0, as.numeric(n_mat) / as.numeric(n_tot), NA_real_),
    stringsAsFactors = FALSE
  )

}


#==============================================================
# 2. Monthly recruitment pattern (TropFishR::recruitment wrapper)
#==============================================================

#' Monthly recruitment pattern (TropFishR)
#'
#' Wrapper around \code{TropFishR::recruitment()} that reconstructs, from
#' an \code{lfq} object (length frequencies + growth parameters), the
#' percentage of recruits per month. The peak(s) indicate the arrival of
#' juveniles, a second target for biological rest.
#'
#' @param lfq \code{lfq} object from \pkg{TropFishR} containing at least
#'   \code{midLengths}, \code{catch}, \code{Linf}, \code{K} and \code{t0}.
#' @param tsample Sampling time (fraction of the year). If \code{NULL},
#'   inferred from the \code{lfq} dates when available.
#' @param ... Arguments passed to \code{TropFishR::recruitment()}.
#'
#' @return \code{data.frame} \code{month, month_label, per_recruits}
#'   (percentage of recruits per month, 12 rows).
#' @references Sparre & Venema (1998), FAO 306/1.
#' @export
recruitment_pattern <- function(lfq, tsample = NULL, ...) {

  if (!requireNamespace("TropFishR", quietly = TRUE))
    stop("The 'TropFishR' package must be installed.", call. = FALSE)

  if (is.null(tsample)) {
    if (!is.null(lfq$dates)) {
      d <- as.POSIXlt(lfq$dates)
      tsample <- d$yday / 365
    } else {
      stop("'tsample' is required (no dates in the lfq object).", call. = FALSE)
    }
  }

  rec <- tryCatch(
    TropFishR::recruitment(param = lfq, tsample = tsample, plot = FALSE, ...),
    error = function(e)
      stop("recruitment: ", conditionMessage(e), call. = FALSE)
  )

  ## per_recruits is a monthly vector (summing to 100), aligned with rec$months
  ## (month numbers 1-12) when t0 is provided.
  pr <- as.numeric(rec$per_recruits)
  mo <- rec$months

  agg <- numeric(12)
  if (!is.null(mo) && length(mo) == length(pr) &&
      all(mo %in% 1:12)) {
    ## direct alignment by month number
    for (i in seq_along(pr))
      agg[mo[i]] <- agg[mo[i]] + pr[i]
  } else if (length(pr) == 12) {
    ## fallback: assumes Jan-Dec order
    agg <- pr
  } else {
    stop("Unexpected recruitment() structure (per_recruits of length ",
         length(pr), ").", call. = FALSE)
  }

  data.frame(
    month        = 1:12,
    month_label  = .month_labels,
    per_recruits = agg,
    stringsAsFactors = FALSE
  )

}


#==============================================================
# 3. Optimal window and duration for biological rest
#==============================================================

#' Determine the optimal window and duration for biological rest
#'
#' Combines the monthly spawning and recruitment signals into a closure
#' score, then identifies the contiguous window of \code{duration} months that
#' maximizes the cumulative protected score. Since the year is cyclical, the
#' search is circular (a closure can overlap December-January). If
#' \code{duration = NULL}, all durations from 1 to 6 months are evaluated and
#' compared.
#'
#' @param spawning \code{data.frame} from \code{spawning_season()}
#'   (\code{prop_mature} column). Optional.
#' @param recruitment \code{data.frame} from \code{recruitment_pattern()}
#'   (\code{per_recruits} column). Optional. At least one of the two signals
#'   is required.
#' @param duration Closure duration in months (integer 1-11) or \code{NULL} to
#'   compare 1 to \code{max_duration} months.
#' @param max_duration Maximum duration evaluated when \code{duration = NULL}
#'   (default 6).
#' @param w_spawning,w_recruitment relative weight of the two signals in the score
#'   (normalized internally). Default 0.5 / 0.5.
#'
#' @return Object of class \code{FishStockClosure}: list containing the monthly
#'   score (\code{$monthly}), the best window for the requested duration
#'   (\code{$best}) and, where applicable, the duration comparison table
#'   (\code{$by_duration}).
#' @seealso \code{\link{spawning_season}}, \code{\link{recruitment_pattern}}
#' @examples
#' \dontrun{
#' sp <- spawning_season(bio, "Date", "maturity", mature_codes = c(2,3,4))
#' rc <- recruitment_pattern(lfq)
#' optimal_closure(spawning = sp, recruitment = rc, duration = 3)
#' }
#' @export
optimal_closure <- function(spawning = NULL,
                            recruitment = NULL,
                            duration = NULL,
                            max_duration = 6,
                            w_spawning = 0.5,
                            w_recruitment = 0.5) {

  if (is.null(spawning) && is.null(recruitment))
    stop("At least one signal (spawning or recruitment) is required.",
         call. = FALSE)

  ## --- Build normalized monthly shares (12) ---
  sp_share <- if (!is.null(spawning))
    .norm_share(spawning$prop_mature) else rep(0, 12)
  rc_share <- if (!is.null(recruitment))
    .norm_share(recruitment$per_recruits) else rep(0, 12)

  ## normalized weights (an absent signal -> weight shifted to the other)
  ws <- if (is.null(spawning)) 0 else w_spawning
  wr <- if (is.null(recruitment)) 0 else w_recruitment
  wtot <- ws + wr
  if (wtot <= 0) stop("Total weight is zero.", call. = FALSE)
  ws <- ws / wtot; wr <- wr / wtot

  score <- ws * sp_share + wr * rc_share
  score <- .norm_share(score)

  monthly <- data.frame(
    month        = 1:12,
    month_label  = .month_labels,
    spawning     = round(sp_share, 4),
    recruitment  = round(rc_share, 4),
    score        = round(score, 4),
    stringsAsFactors = FALSE
  )

  ## --- Search for the maximum circular window for a given duration ---
  best_window <- function(d) {
    d <- as.integer(d)
    if (d < 1 || d > 11)
      stop("'duration' must be between 1 and 11 months.", call. = FALSE)
    cov <- vapply(1:12, function(start) {
      idx <- ((start - 1 + 0:(d - 1)) %% 12) + 1
      sum(score[idx])
    }, numeric(1))
    start <- which.max(cov)
    idx <- ((start - 1 + 0:(d - 1)) %% 12) + 1
    list(
      duration     = d,
      start_month  = start,
      months       = idx,
      months_label = paste(.month_labels[idx], collapse = "-"),
      protected    = cov[start]     # fraction of the score protected
    )
  }

  by_duration <- NULL
  if (is.null(duration)) {
    max_duration <- min(as.integer(max_duration), 11)
    rows <- lapply(1:max_duration, function(d) {
      b <- best_window(d)
      data.frame(duration = b$duration, start_month = b$start_month,
                 months = b$months_label, protected = round(b$protected, 4),
                 stringsAsFactors = FALSE)
    })
    by_duration <- do.call(rbind, rows)
    ## "recommended" duration: smallest duration protecting >= 60% of the score,
    ## otherwise the one with the highest marginal gain (elbow).
    thr <- which(by_duration$protected >= 0.60)
    rec_d <- if (length(thr)) min(thr) else {
      gains <- diff(c(0, by_duration$protected))
      which.max(gains)
    }
    best <- best_window(by_duration$duration[rec_d])
  } else {
    best <- best_window(duration)
  }

  structure(
    list(
      monthly     = monthly,
      best        = best,
      by_duration = by_duration,
      weights     = c(spawning = ws, recruitment = wr)
    ),
    class = "FishStockClosure"
  )

}


#==============================================================
# S3 Methods
#==============================================================

#' @export
print.FishStockClosure <- function(x, ...) {

  cat("\n=====================================\n")
  cat(" stockflow: Optimal biological rest\n")
  cat("=====================================\n\n")

  cat("Monthly closure score (share protected per month):\n")
  print(x$monthly[, c("month_label", "spawning", "recruitment", "score")],
        row.names = FALSE)

  if (!is.null(x$by_duration)) {
    cat("\nDuration comparison:\n")
    print(x$by_duration, row.names = FALSE)
  }

  b <- x$best
  cat(sprintf(
    "\nRecommended window: %d months (%s), from month %d.\n",
    b$duration, b$months_label, b$start_month))
  cat(sprintf("Share of score protected: %.1f %%\n", 100 * b$protected))

  invisible(x)

}

#' Chart of the monthly closure score with optimal window
#' @param x \code{FishStockClosure} object.
#' @param ... Ignored.
#' @export
plot.FishStockClosure <- function(x, ...) {

  df <- x$monthly
  protected_months <- x$best$months

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cols <- ifelse(df$month %in% protected_months, "#c0392b", "grey70")
    graphics::barplot(df$score, names.arg = df$month_label, col = cols,
                      ylab = "Closure score",
                      main = "Biological rest: monthly score")
    return(invisible(x))
  }

  df$protege <- ifelse(df$month %in% protected_months,
                       "Proposed window", "Outside closure")
  df$month_label <- factor(df$month_label, levels = .month_labels)

  p <- ggplot2::ggplot(df, ggplot2::aes(month_label, score, fill = protege)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("Proposed window" = "#c0392b",
                                          "Outside closure" = "grey70")) +
    ggplot2::labs(
      title = "Biological rest: monthly score and optimal window",
      subtitle = sprintf("Window: %s (%d months, %.0f %% of score protected)",
                         x$best$months_label, x$best$duration,
                         100 * x$best$protected),
      x = NULL, y = "Closure score (monthly share)",
      fill = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  print(p)
  invisible(p)

}