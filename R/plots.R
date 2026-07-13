###############################################################
#
# stockflow
#
# Module : plots.R
#
# Summary plots for stock assessment from length
# frequencies:
#   - length histograms (by sex, by month)
#   - length-weight relationship (W = a L^b)
#   - von Bertalanffy growth curve (VBGF)
#   - maturity ogive (L50)
#   - catch curve (total mortality Z)
#   - sex-ratio (overall and by length class)
#   - condition index (Fulton K, Le Cren Kn)
#   - catch / CPUE series
#
# All functions return a 'ggplot' object (not printed),
# to remain composable and testable. Scientific fits
# rely on official functions
# (TropFishR::VBGF, stats::glm for the logistic ogive, etc.).
#
# Style inspired by earlier fisheries analysis work.
#
###############################################################


#==============================================================
# Theme and internal utilities
#==============================================================

## Common publication theme (classic, legend at bottom).
.theme_fish <- function(base = 12) {
  ggplot2::theme_classic(base_size = base) +
    ggplot2::theme(
      strip.background   = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.text         = ggplot2::element_text(face = "bold"),
      panel.grid.major.y = ggplot2::element_line(colour = "grey88"),
      plot.title         = ggplot2::element_text(face = "bold", size = base + 1),
      plot.subtitle      = ggplot2::element_text(colour = "grey30"),
      legend.position    = "bottom"
    )
}

## Sex palette consistent across all figures.
.pal_sexe <- c(F = "#c2185b", M = "#1565c0", "F" = "#c2185b", "M" = "#1565c0")

## Month abbreviations (consistent with closure.R).
.mois_abbr <- c("Jan","Fev","Mar","Avr","Mai","Jun",
                "Jui","Aou","Sep","Oct","Nov","Dec")

## Checks that a column is present, clear message otherwise.
.need_col <- function(data, col, arg) {
  if (!col %in% names(data))
    stop("Column '", col, "' missing (argument ", arg, ").", call. = FALSE)
  invisible(TRUE)
}

## Wilson binomial confidence interval (robust to small sample sizes
## and extreme proportions, unlike the Wald approximation).
.wilson_ci <- function(x, n, level = 0.95) {
  if (n == 0) return(c(NA_real_, NA_real_))
  z  <- stats::qnorm(1 - (1 - level) / 2)
  ph <- x / n
  den <- 1 + z^2 / n
  ctr <- (ph + z^2 / (2 * n)) / den
  hw  <- z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2)) / den
  c(max(0, ctr - hw), min(1, ctr + hw))
}


#==============================================================
# 1. Length frequency histogram
#==============================================================

#' Length frequency histogram
#'
#' Length distribution, optionally coloured by sex and/or faceted
#' by month (useful for visualizing recruitment and length structure).
#'
#' @param data data.frame of individual measurements.
#' @param length_col Name of the length column.
#' @param sex_col Name of the sex column (optional).
#' @param date_col Name of a date column to facet by month (optional).
#' @param binwidth Bin width (default: range / 40).
#' @param title,xlab Title and label of the length axis.
#'
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' plot_length_frequency(cymbium_bio, "LCQ", sex_col = "sex")
#' }
#' @export
plot_length_frequency <- function(data,
                                  length_col,
                                  sex_col = NULL,
                                  date_col = NULL,
                                  binwidth = NULL,
                                  title = "Length frequency distribution",
                                  xlab = NULL) {

  .need_col(data, length_col, "length_col")
  d <- data[is.finite(data[[length_col]]), , drop = FALSE]
  if (nrow(d) == 0) stop("No valid length.", call. = FALSE)

  if (is.null(binwidth))
    binwidth <- diff(range(d[[length_col]], na.rm = TRUE)) / 40
  if (is.null(xlab)) xlab <- length_col

  aes_hist <- if (!is.null(sex_col) && sex_col %in% names(d))
    ggplot2::aes(x = .data[[length_col]], fill = .data[[sex_col]])
  else
    ggplot2::aes(x = .data[[length_col]])

  p <- ggplot2::ggplot(d, aes_hist) +
    ggplot2::geom_histogram(binwidth = binwidth, colour = "white",
                            linewidth = 0.1, position = "stack")

  if (!is.null(sex_col) && sex_col %in% names(d))
    p <- p + ggplot2::scale_fill_manual(values = .pal_sexe, name = "sex",
                                        na.value = "grey70")

  if (!is.null(date_col) && date_col %in% names(d)) {
    d$.mois <- as.integer(format(as.Date(d[[date_col]]), "%m"))
    p$data <- d
    p <- p + ggplot2::facet_wrap(~ factor(.data[[".mois"]], levels = 1:12,
                                          labels = .mois_abbr),
                                 scales = "free_y")
  }

  p +
    ggplot2::labs(title = title, x = xlab, y = "Count") +
    .theme_fish()

}


#' Temporal evolution of length frequencies
#'
#' Visualizes the progression of length distributions over time, by
#' month or by year, as faceted histograms or overlaid density
#' curves. The density view highlights modal progression
#' (cohort shift), useful for spotting recruitment and
#' documenting length structure before/after a closure.
#'
#' @param data data.frame of individual measurements.
#' @param length_col Length column.
#' @param date_col Date column.
#' @param period Time scale: \code{"month"} (calendar months,
#'   aggregated across all years) or \code{"year"} (years).
#' @param type Representation: \code{"histogram"} (faceted by period) or
#'   \code{"density"} (overlaid curves, one per period).
#' @param binwidth Bin width for histograms (default: range / 40).
#' @param sex_col Sex column (optional, fills the histograms).
#' @param title,xlab Title and label of the length axis.
#'
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plot_length_frequency}}
#' @examples
#' \dontrun{
#' plot_length_evolution(cymbium_bio, "LCQ", "Date",
#'                       period = "month", type = "density")
#' plot_length_evolution(cymbium_bio, "LCQ", "Date",
#'                       period = "year", type = "histogram")
#' }
#' @export
plot_length_evolution <- function(data,
                                  length_col,
                                  date_col,
                                  period = c("month", "year"),
                                  type = c("histogram", "density"),
                                  binwidth = NULL,
                                  sex_col = NULL,
                                  title = NULL,
                                  xlab = NULL) {

  period <- match.arg(period)
  type   <- match.arg(type)
  .need_col(data, length_col, "length_col")
  .need_col(data, date_col, "date_col")

  d <- data[is.finite(data[[length_col]]), , drop = FALSE]
  dt <- as.Date(d[[date_col]])
  d  <- d[!is.na(dt), , drop = FALSE]
  if (nrow(d) == 0) stop("No valid dated length.", call. = FALSE)

  ## ordered period factor
  if (period == "month") {
    m <- as.integer(format(as.Date(d[[date_col]]), "%m"))
    d$.periode <- factor(.mois_abbr[m], levels = .mois_abbr)
    per_lab <- "Month"
  } else {
    y <- format(as.Date(d[[date_col]]), "%Y")
    d$.periode <- factor(y, levels = sort(unique(y)))
    per_lab <- "year"
  }
  d <- d[!is.na(d$.periode), , drop = FALSE]

  if (is.null(binwidth))
    binwidth <- diff(range(d[[length_col]], na.rm = TRUE)) / 40
  if (is.null(xlab)) xlab <- length_col
  if (is.null(title))
    title <- sprintf("Length evolution by %s",
                     if (period == "month") "month" else "year")

  if (type == "histogram") {
    aes_h <- if (!is.null(sex_col) && sex_col %in% names(d))
      ggplot2::aes(x = .data[[length_col]], fill = .data[[sex_col]])
    else ggplot2::aes(x = .data[[length_col]])

    p <- ggplot2::ggplot(d, aes_h) +
      ggplot2::geom_histogram(binwidth = binwidth, colour = "white",
                              linewidth = 0.1, position = "stack") +
      ggplot2::facet_wrap(ggplot2::vars(.periode), scales = "free_y")

    if (!is.null(sex_col) && sex_col %in% names(d))
      p <- p + ggplot2::scale_fill_manual(values = .pal_sexe, name = "sex",
                                          na.value = "grey70")

    return(p +
             ggplot2::labs(title = title, x = xlab, y = "Count") +
             .theme_fish())
  }

  ## type == "density": overlaid curves, one per period
  ggplot2::ggplot(d, ggplot2::aes(x = .data[[length_col]],
                                  colour = .periode)) +
    ggplot2::geom_density(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::scale_colour_viridis_d(name = per_lab, option = "C",
                                    end = 0.92) +
    ggplot2::labs(title = title, x = xlab, y = "Density") +
    .theme_fish()

}


#==============================================================
# 2. Length-weight relationship
#==============================================================

#' Length-weight relationship (W = a L^b)
#'
#' Length-weight scatter plot with power law fit \eqn{W = a L^b}
#' (linear regression on logarithms). Can be plotted in
#' natural scale (the fitted power curve is overlaid), in
#' log-log scale (the relationship becomes linear), or both side
#' by side.
#'
#' @param data data.frame of individual measurements.
#' @param length_col,weight_col Length and weight columns.
#' @param scale Plot scale: \code{"both"} (default, natural
#'   scale and log-log side by side), \code{"power"} (natural scale
#'   with the \eqn{W = a L^b} curve) or \code{"log"} (linear log-log).
#' @param log_scale Deprecated: kept for compatibility. If provided,
#'   \code{TRUE} forces \code{scale = "log"} and \code{FALSE} \code{scale =
#'   "power"} (ignored if \code{scale} is explicitly given).
#' @param level Confidence level of the plotted band (default 0.95).
#' @param title Plot title.
#'
#' @return A \code{ggplot} object (or a \code{patchwork} assembly for
#'   \code{scale = "both"}). The fitted parameters (a, b with its CI, R2)
#'   are annotated and a confidence band surrounds the curve.
#' @seealso \code{\link{fit_length_weight}} for the full fit (FSA).
#' @examples
#' \dontrun{
#' plot_length_weight(cymbium_bio, "LCQ", "weight")               # both
#' plot_length_weight(cymbium_bio, "LCQ", "weight", scale = "power")
#' }
#' @export
plot_length_weight <- function(data,
                               length_col,
                               weight_col,
                               scale = c("both", "power", "log"),
                               log_scale = NULL,
                               level = 0.95,
                               title = "Length-weight relationship (W = a L^b)") {

  ## backward compatibility: the old log_scale argument takes priority
  ## only if 'scale' is not explicitly provided by the caller.
  if (missing(scale) && !is.null(log_scale))
    scale <- if (isTRUE(log_scale)) "log" else "power"
  scale <- match.arg(scale)

  .need_col(data, length_col, "length_col")
  .need_col(data, weight_col, "weight_col")

  d <- data[is.finite(data[[length_col]]) & is.finite(data[[weight_col]]) &
            data[[length_col]] > 0 & data[[weight_col]] > 0, , drop = FALSE]
  if (nrow(d) < 3) stop("Too few valid points for the fit.",
                        call. = FALSE)

  ## Single power law fit on the logarithms.
  L <- d[[length_col]]; W <- d[[weight_col]]
  m  <- stats::lm(log(W) ~ log(L))
  a  <- exp(stats::coef(m)[1])
  b  <- stats::coef(m)[2]
  r2 <- summary(m)$r.squared
  ## 95 percent confidence interval of the exponent b
  b_ci <- stats::confint(m, "log(L)", level = level)
  lab <- sprintf("a = %.4g\nb = %.2f [%.2f - %.2f]\nR2 = %.3f",
                 a, b, b_ci[1], b_ci[2], r2)

  ## Fitted power curve + confidence band: predicted on the log
  ## scale (where the CI is valid) then back-transformed to natural scale.
  Lseq  <- seq(min(L), max(L), length.out = 200)
  pred  <- stats::predict(m, newdata = data.frame(L = Lseq),
                          interval = "confidence", level = level)
  curve <- data.frame(.L = Lseq,
                      .W  = exp(pred[, "fit"]),
                      .lo = exp(pred[, "lwr"]),
                      .hi = exp(pred[, "upr"]))

  ci_pct <- round(100 * level)

  ## --- natural scale panel (power law) ---
  p_pow <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[length_col]],
                                           y = .data[[weight_col]])) +
    ggplot2::geom_point(alpha = 0.25, size = 0.7, colour = "#2166ac") +
    ggplot2::geom_ribbon(data = curve,
                         ggplot2::aes(x = .L, ymin = .lo, ymax = .hi),
                         fill = "#c0392b", alpha = 0.2, inherit.aes = FALSE) +
    ggplot2::geom_line(data = curve, ggplot2::aes(.L, .W),
                       colour = "#c0392b", linewidth = 0.9,
                       inherit.aes = FALSE) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = lab,
                      hjust = -0.15, vjust = 1.2, size = 3.3) +
    ggplot2::labs(title = if (scale == "both") "Natural scale" else title,
                  subtitle = sprintf("Curve W = a L^b (CI%d %% band)", ci_pct),
                  x = length_col, y = weight_col) +
    .theme_fish()

  ## --- log-log panel (linear relationship) ---
  p_log <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[length_col]],
                                           y = .data[[weight_col]])) +
    ggplot2::geom_point(alpha = 0.25, size = 0.7, colour = "#2166ac") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, colour = "black",
                         fill = "grey60", linewidth = 0.7, se = TRUE,
                         level = level) +
    ## Finite coordinates (derived from data): under scale_*_log10(),
    ## a -Inf position triggers log10(-Inf) -> "NaNs produced".
    ggplot2::annotate("text",
                      x = min(d[[length_col]], na.rm = TRUE),
                      y = max(d[[weight_col]], na.rm = TRUE),
                      label = lab, hjust = -0.05, vjust = 1.2, size = 3.3) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::labs(title = if (scale == "both") "Logarithmic scales" else title,
                  subtitle = sprintf("Logarithmic scales (CI%d %% band)",
                                     ci_pct),
                  x = length_col, y = weight_col) +
    .theme_fish()

  if (scale == "power") return(p_pow)
  if (scale == "log")   return(p_log)

  ## scale == "both": side-by-side assembly (patchwork if available).
  if (requireNamespace("patchwork", quietly = TRUE))
    return(patchwork::wrap_plots(p_pow, p_log, ncol = 2) +
             patchwork::plot_annotation(title = title))
  ## fallback: list of the two panels if patchwork is absent.
  list(power = p_pow, log = p_log)

}


#==============================================================
# 3. Von Bertalanffy growth curve
#==============================================================

#' Von Bertalanffy growth curve (VBGF)
#'
#' Plots \eqn{L(t) = L_\infty (1 - e^{-K(t - t_0)})} from growth
#' parameters, relying on \code{TropFishR::VBGF()}. Several parameter
#' sets can be compared (e.g. ELEFAN GA vs SA).
#'
#' @param Linf,K,t0 Growth parameters. \code{Linf} and \code{K} can
#'   be vectors (compares several curves); \code{t0} default 0.
#' @param labels Curve labels (if parameter vectors).
#' @param tmax Maximum age to plot (default \eqn{3/\min(K)}).
#' @param length_unit Length unit for the axis.
#' @param title Plot title.
#'
#' @return A \code{ggplot} object.
#' @references von Bertalanffy (1938); TropFishR::VBGF.
#' @examples
#' \dontrun{
#' plot_vbgf(Linf = c(42, 45), K = c(0.45, 0.38),
#'           labels = c("ELEFAN_SA", "ELEFAN_GA"))
#' }
#' @export
plot_vbgf <- function(Linf, K, t0 = 0,
                      labels = NULL,
                      tmax = NULL,
                      length_unit = "cm",
                      title = "Von Bertalanffy growth") {

  if (!requireNamespace("TropFishR", quietly = TRUE))
    stop("The 'TropFishR' package must be installed.", call. = FALSE)

  n <- max(length(Linf), length(K), length(t0))
  Linf <- rep_len(Linf, n); K <- rep_len(K, n); t0 <- rep_len(t0, n)
  if (is.null(labels)) labels <- if (n == 1) "VBGF" else paste0("Set ", seq_len(n))
  labels <- rep_len(labels, n)
  if (is.null(tmax)) tmax <- 3 / min(K)

  tgrid <- seq(0, tmax, length.out = 200)
  curves <- do.call(rbind, lapply(seq_len(n), function(i) {
    L <- TropFishR::VBGF(list(Linf = Linf[i], K = K[i], t0 = t0[i]), t = tgrid)
    data.frame(age = tgrid, L = L, jeu = labels[i], stringsAsFactors = FALSE)
  }))

  lab <- do.call(rbind, lapply(seq_len(n), function(i)
    data.frame(jeu = labels[i],
               txt = sprintf("Linf = %.1f %s ; K = %.3f/year",
                             Linf[i], length_unit, K[i]))))

  p <- ggplot2::ggplot(curves, ggplot2::aes(age, L, colour = jeu)) +
    ggplot2::geom_line(linewidth = 0.9)

  if (n == 1) {
    p <- p + ggplot2::annotate("text", x = 0, y = Inf, label = lab$txt[1],
                               hjust = -0.05, vjust = 1.4, size = 3.2) +
      ggplot2::guides(colour = "none")
  } else {
    p <- p + ggplot2::scale_colour_brewer(palette = "Dark2", name = NULL)
  }

  p +
    ggplot2::labs(title = title,
                  subtitle = "L(t) = Linf (1 - e^(-K(t - t0)))",
                  x = "Relative age (years)",
                  y = paste0("Length (", length_unit, ")")) +
    .theme_fish()

}


#==============================================================
# 4. Maturity ogive (L50)
#==============================================================

#' Maturity ogive and size at 50 percent (L50)
#'
#' Fits a logistic ogive of the proportion mature as a function of
#' size (binomial \code{stats::glm}) and extracts \eqn{L_{50}} and
#' \eqn{L_{95}}. Plots the fitted curve and the observed proportions per
#' length class.
#'
#' @param data data.frame of individual measurements.
#' @param length_col Length column.
#' @param maturity_col Maturity column (codes or logical).
#' @param mature_codes Values of \code{maturity_col} corresponding to
#'   mature individuals. If \code{NULL} and the column is logical, \code{TRUE}.
#' @param bin_width Class width for the observed points (default:
#'   range / 15).
#' @param level Confidence level of the bands and of the CI on L50 (default 0.95).
#' @param title Plot title.
#'
#' @return A \code{ggplot} object. \eqn{L_{50}} (with its confidence
#'   interval) is shown as an annotation; a confidence band surrounds
#'   the fitted ogive.
#' @references King (2007), Fisheries Biology.
#' @examples
#' \dontrun{
#' plot_maturity_ogive(cymbium_bio, "LCQ", "maturity",
#'                     mature_codes = c("M", "P"))
#' }
#' @export
plot_maturity_ogive <- function(data,
                                length_col,
                                maturity_col,
                                mature_codes = NULL,
                                bin_width = NULL,
                                level = 0.95,
                                title = "Maturity ogive") {

  .need_col(data, length_col, "length_col")
  .need_col(data, maturity_col, "maturity_col")

  d <- data[is.finite(data[[length_col]]) & !is.na(data[[maturity_col]]),
            , drop = FALSE]
  if (nrow(d) < 10) stop("Too few observations to fit the ogive.",
                         call. = FALSE)

  mat <- d[[maturity_col]]
  d$.mature <- if (is.logical(mat)) as.integer(mat)
    else {
      if (is.null(mature_codes))
        stop("'mature_codes' is required for a non-logical column.",
             call. = FALSE)
      as.integer(as.character(mat) %in% as.character(mature_codes))
    }

  L <- d[[length_col]]
  fit <- stats::glm(d$.mature ~ L, family = stats::binomial())
  co <- stats::coef(fit)
  L50 <- -co[1] / co[2]
  L95 <- (log(0.95 / 0.05) - co[1]) / co[2]

  ## CI on L50 by delta method: L50 = -b0/b1, gradient g = (-1/b1, b0/b1^2)
  V  <- stats::vcov(fit)
  g  <- c(-1 / co[2], co[1] / co[2]^2)
  se_L50 <- sqrt(as.numeric(t(g) %*% V %*% g))
  zc <- stats::qnorm(1 - (1 - level) / 2)
  L50_lo <- L50 - zc * se_L50
  L50_hi <- L50 + zc * se_L50

  ## Confidence band of the predicted proportion (computed on the
  ## logit link scale then back-transformed, which guarantees p in [0, 1]).
  grid <- data.frame(L = seq(min(L), max(L), length.out = 200))
  pr   <- stats::predict(fit, newdata = grid, type = "link", se.fit = TRUE)
  ilink <- stats::family(fit)$linkinv
  grid$p  <- ilink(pr$fit)
  grid$lo <- ilink(pr$fit - zc * pr$se.fit)
  grid$hi <- ilink(pr$fit + zc * pr$se.fit)

  if (is.null(bin_width)) bin_width <- diff(range(L)) / 15
  d$.bin <- cut(L, breaks = seq(min(L), max(L) + bin_width, by = bin_width),
                include.lowest = TRUE)
  obs <- stats::aggregate(cbind(.mature, n = 1) ~ .bin, data = d, FUN = sum)
  bins <- as.numeric(sub("[^0-9.]*([0-9.]+).*", "\\1", obs$.bin)) + bin_width / 2
  obs$mid <- bins
  obs$prop <- obs$.mature / obs$n

  ci_pct <- round(100 * level)

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = grid, ggplot2::aes(L, ymin = lo, ymax = hi),
                         fill = "#c0392b", alpha = 0.18) +
    ggplot2::geom_point(data = obs,
                        ggplot2::aes(mid, prop, size = n),
                        colour = "#2166ac", alpha = 0.6) +
    ggplot2::geom_line(data = grid, ggplot2::aes(L, p),
                       colour = "#c0392b", linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted",
                        colour = "grey45") +
    ggplot2::geom_vline(xintercept = L50, linetype = "dashed",
                        colour = "#c0392b") +
    ggplot2::annotate("rect", xmin = L50_lo, xmax = L50_hi,
                      ymin = -Inf, ymax = Inf, fill = "#c0392b", alpha = 0.08) +
    ggplot2::annotate("text", x = L50, y = 0.05,
                      label = sprintf("L50 = %.1f [%.1f - %.1f]",
                                      L50, L50_lo, L50_hi),
                      hjust = -0.05, size = 3.2, colour = "#c0392b") +
    ggplot2::scale_size_continuous(name = "n", range = c(1, 5)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = title,
                  subtitle = sprintf("L50 = %.1f (CI%d %% : %.1f - %.1f) ; L95 = %.1f",
                                     L50, ci_pct, L50_lo, L50_hi, L95),
                  x = length_col, y = "Proportion mature") +
    .theme_fish()

  p

}


#==============================================================
# 5. Catch curve (total mortality Z)
#==============================================================

#' Linear catch curve (total mortality Z)
#'
#' Represents the length-converted catch curve: logarithm of
#' the number per class (corrected for residence time) as a function of
#' relative age, with the descending part used to estimate \eqn{Z}. This
#' function relies on an object produced by \code{\link{estimate_catchcurve}}.
#'
#' @param x \code{FishStockZ} object (see \code{\link{estimate_catchcurve}}).
#' @param level Confidence level of the regression band (default 0.95).
#' @param title Plot title.
#'
#' @return A \code{ggplot} object. The slope of the descending part estimates
#'   \eqn{Z}; a confidence band surrounds the regression line and a
#'   CI on \eqn{Z} is annotated when it can be recomputed.
#' @seealso \code{\link{estimate_catchcurve}}
#' @export
plot_catch_curve <- function(x,
                             level = 0.95,
                             title = "Catch curve (total mortality Z)") {

  if (!inherits(x, "FishStockZ"))
    stop("'x' must be a FishStockZ object (estimate_catchcurve()).",
         call. = FALSE)

  cc <- x$catch_curve %||% x$cc %||% x$model
  ## TropFishR::catchCurve structure: $tplot / $lnC_plot or $t_midL / $lnC
  tt  <- cc$t_midL %||% cc$tplot %||% cc$ln_x
  lnc <- cc$lnC_dt %||% cc$lnC %||% cc$ln_y
  reg <- cc$reg_int

  if (is.null(tt) || is.null(lnc))
    stop("Unrecognized catch curve structure in the object.",
         call. = FALSE)

  df <- data.frame(t = as.numeric(tt), lnC = as.numeric(lnc))
  df <- df[is.finite(df$t) & is.finite(df$lnC), ]
  df$used <- FALSE
  if (!is.null(reg) && length(reg) == 2)
    df$used <- seq_len(nrow(df)) >= reg[1] & seq_len(nrow(df)) <= reg[2]

  Zval <- x$Z %||% x$Z$Z

  ## CI on Z (= -slope) recomputed from the regression points.
  z_lab <- if (!is.null(Zval)) sprintf("Z = %.3f /year", Zval) else NULL
  used_df <- df[df$used, , drop = FALSE]
  if (nrow(used_df) >= 3) {
    fitz <- stats::lm(lnC ~ t, data = used_df)
    ci   <- stats::confint(fitz, "t", level = level)
    ci_pct <- round(100 * level)
    z_lab <- sprintf("Z = %.3f /year (CI%d %% : %.3f - %.3f)",
                     -stats::coef(fitz)[2], ci_pct, -ci[2], -ci[1])
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(t, lnC)) +
    ggplot2::geom_point(ggplot2::aes(colour = used), size = 2) +
    ggplot2::geom_smooth(data = used_df,
                         method = "lm", formula = y ~ x,
                         colour = "#c0392b", fill = "#c0392b", alpha = 0.18,
                         se = TRUE, level = level, linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = c("TRUE" = "#c0392b",
                                            "FALSE" = "grey65"),
                                 labels = c("TRUE" = "Regression",
                                            "FALSE" = "Excluded"),
                                 name = NULL) +
    ggplot2::labs(title = title, subtitle = z_lab,
                  x = "Relative age / class", y = "ln(N / dt)") +
    .theme_fish()

  p

}


#==============================================================
# 6. Sex-ratio
#==============================================================

#' Overall sex-ratio and by length class
#'
#' Proportion of females (and males) by length class, with the
#' overall proportion as reference. Useful for detecting size
#' dimorphism or spatial segregation.
#'
#' @param data data.frame of individual measurements.
#' @param sex_col Sex column ("M"/"F" values).
#' @param length_col Length column (optional: if absent, only the
#'   overall sex-ratio is plotted).
#' @param bin_width Length class width.
#' @param level Confidence level of the binomial intervals (default 0.95).
#' @param title Plot title.
#'
#' @return A \code{ggplot} object. Each proportion is shown with its
#'   binomial confidence interval (Wilson method).
#' @examples
#' \dontrun{
#' plot_sex_ratio(cymbium_bio, "sex", "LCQ")
#' }
#' @export
plot_sex_ratio <- function(data,
                           sex_col,
                           length_col = NULL,
                           bin_width = NULL,
                           level = 0.95,
                           title = "Sex-ratio") {

  .need_col(data, sex_col, "sex_col")
  d <- data[toupper(as.character(data[[sex_col]])) %in% c("M", "F"), ,
            drop = FALSE]
  d$.sexe <- toupper(as.character(d[[sex_col]]))
  if (nrow(d) == 0) stop("No sexed individual (M/F).", call. = FALSE)

  prop_f_global <- mean(d$.sexe == "F")
  ci_pct <- round(100 * level)

  if (is.null(length_col) || !length_col %in% names(d)) {
    nf <- sum(d$.sexe == "F"); ntot <- nrow(d)
    ci <- .wilson_ci(nf, ntot, level)
    tab <- as.data.frame(table(sex = d$.sexe))
    ## binomial CI on the proportion of females, positioned on the F bar
    fbar <- data.frame(sex = "F", lo = ci[1] * ntot, hi = ci[2] * ntot)
    p <- ggplot2::ggplot(tab, ggplot2::aes(sex, Freq, fill = sex)) +
      ggplot2::geom_col(width = 0.6) +
      ggplot2::geom_errorbar(data = fbar,
                             ggplot2::aes(x = sex, ymin = lo, ymax = hi),
                             inherit.aes = FALSE, width = 0.2) +
      ggplot2::scale_fill_manual(values = .pal_sexe, guide = "none") +
      ggplot2::labs(title = title,
                    subtitle = sprintf("Proportion of females = %.1f %% (CI%d %% : %.1f - %.1f)",
                                       100 * prop_f_global, ci_pct,
                                       100 * ci[1], 100 * ci[2]),
                    x = NULL, y = "Count") +
      .theme_fish()
    return(p)
  }

  d <- d[is.finite(d[[length_col]]), , drop = FALSE]
  L <- d[[length_col]]
  if (is.null(bin_width)) bin_width <- diff(range(L)) / 15
  brks <- seq(floor(min(L)), max(L) + bin_width, by = bin_width)
  mids <- brks[-length(brks)] + bin_width / 2
  d$.bin <- cut(L, breaks = brks, include.lowest = TRUE)
  agg <- stats::aggregate(cbind(f = .sexe == "F", n = 1) ~ .bin, data = d,
                          FUN = sum)
  ## the midpoint is read from the class level (aggregate drops
  ## empty classes: 'mids' cannot be reused by position).
  agg$mid <- mids[match(agg$.bin, levels(d$.bin))]
  agg$prop_f <- agg$f / agg$n
  ## Wilson binomial interval per class
  ci_bin <- t(mapply(function(x, nn) .wilson_ci(x, nn, level),
                     agg$f, agg$n))
  agg$lo <- ci_bin[, 1]; agg$hi <- ci_bin[, 2]

  p <- ggplot2::ggplot(agg, ggplot2::aes(mid, prop_f)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted",
                        colour = "grey45") +
    ggplot2::geom_hline(yintercept = prop_f_global, linetype = "dashed",
                        colour = "#c2185b") +
    ggplot2::geom_col(ggplot2::aes(fill = prop_f), width = bin_width * 0.9) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi),
                           width = bin_width * 0.4, colour = "grey25",
                           alpha = 0.8) +
    ggplot2::geom_point(ggplot2::aes(size = n), colour = "grey25",
                        alpha = 0.5) +
    ggplot2::scale_fill_gradient2(low = "#1565c0", mid = "grey85",
                                  high = "#c2185b", midpoint = 0.5,
                                  name = "Prop. F", limits = c(0, 1)) +
    ggplot2::scale_size_continuous(name = "n", range = c(0.5, 4)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = title,
                  subtitle = sprintf("Proportion of females by class (CI%d %% binomial) ; overall = %.1f %%",
                                     ci_pct, 100 * prop_f_global),
                  x = length_col, y = "Proportion of females") +
    .theme_fish()

  p

}


#==============================================================
# 7. Condition index
#==============================================================

#' Condition index (Fulton K or Le Cren Kn)
#'
#' Computes and plots an individual condition index. The Fulton factor
#' \eqn{K = 100000 \cdot W / L^3} is comparable within the same length
#' unit; Le Cren's relative factor \eqn{K_n = W / (a L^b)} is
#' dimensionless and comparable across stocks. The distribution can be
#' shown by month (seasonal variation in condition).
#'
#' @param data data.frame of individual measurements.
#' @param length_col,weight_col Length and weight columns.
#' @param index Index type: \code{"lecren"} (default) or \code{"fulton"}.
#' @param date_col Date column for monthly variation (optional).
#' @param title Plot title.
#'
#' @return A \code{ggplot} object.
#' @references Le Cren (1951); Fulton (1904).
#' @examples
#' \dontrun{
#' plot_condition(cymbium_bio, "LCQ", "weight", date_col = "Date")
#' }
#' @export
plot_condition <- function(data,
                           length_col,
                           weight_col,
                           index = c("lecren", "fulton"),
                           date_col = NULL,
                           title = NULL) {

  index <- match.arg(index)
  .need_col(data, length_col, "length_col")
  .need_col(data, weight_col, "weight_col")

  d <- data[is.finite(data[[length_col]]) & is.finite(data[[weight_col]]) &
            data[[length_col]] > 0 & data[[weight_col]] > 0, , drop = FALSE]
  if (nrow(d) < 3) stop("Too few valid points.", call. = FALSE)

  L <- d[[length_col]]; W <- d[[weight_col]]

  if (index == "fulton") {
    d$.cond <- 1e5 * W / L^3
    ylab <- "Fulton factor K"
    ref  <- NA_real_
  } else {
    m <- stats::lm(log(W) ~ log(L))
    a <- exp(stats::coef(m)[1]); b <- stats::coef(m)[2]
    d$.cond <- W / (a * L^b)
    ylab <- "Le Cren factor Kn"
    ref  <- 1
  }
  if (is.null(title))
    title <- paste0("Condition index (",
                    if (index == "fulton") "Fulton K" else "Le Cren Kn", ")")

  if (!is.null(date_col) && date_col %in% names(d)) {
    d$.mois <- as.integer(format(as.Date(d[[date_col]]), "%m"))
    d <- d[!is.na(d$.mois), ]
    agg <- stats::aggregate(.cond ~ .mois, data = d, FUN = function(x) {
      c(m = mean(x), se = stats::sd(x) / sqrt(length(x)))
    })
    md <- data.frame(mois = agg$.mois, moy = agg$.cond[, "m"],
                     se = agg$.cond[, "se"])
    md$lo <- md$moy - 1.96 * md$se; md$hi <- md$moy + 1.96 * md$se

    p <- ggplot2::ggplot(md, ggplot2::aes(mois, moy))
    if (is.finite(ref))
      p <- p + ggplot2::geom_hline(yintercept = ref, linetype = "dotted",
                                   colour = "grey45")
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                           fill = "#2166ac", alpha = 0.2) +
      ggplot2::geom_line(colour = "#2166ac", linewidth = 0.8) +
      ggplot2::geom_point(colour = "#2166ac", size = 1.8) +
      ggplot2::scale_x_continuous(breaks = 1:12, labels = .mois_abbr) +
      ggplot2::labs(title = title,
                    subtitle = "Monthly variation (mean +/- 95% CI)",
                    x = "Month", y = ylab) +
      .theme_fish()
    return(p)
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(.cond))
  if (is.finite(ref))
    p <- p + ggplot2::geom_vline(xintercept = ref, linetype = "dotted",
                                 colour = "grey45")
  p +
    ggplot2::geom_histogram(bins = 40, fill = "#2166ac", colour = "white",
                            linewidth = 0.1) +
    ggplot2::labs(title = title, x = ylab, y = "Count") +
    .theme_fish()

}


#==============================================================
# 8. Gonadosomatic ratio (GSR / GSI)
#==============================================================

#' Gonadosomatic ratio (GSR / gonadosomatic index)
#'
#' Plots the monthly evolution of the gonadosomatic ratio
#' \eqn{GSR = 100 \times P_{gonad} / P_{somatic}}, monthly mean
#' surrounded by its confidence interval. The GSR peak marks the
#' period of gonadal maturation / spawning and complements the biological
#' rest analysis (see \code{\link{spawning_season}}). An artifact filter
#' discards GSR values deemed physiologically improbable.
#'
#' @param data data.frame of individual measurements.
#' @param gonad_col Gonad weight column (same unit as
#'   \code{weight_col}).
#' @param weight_col Somatic / total weight column of the individual.
#' @param date_col Date column (for the monthly evolution). If
#'   \code{NULL}, only the overall distribution is plotted.
#' @param sex_col Sex column (optional). If provided with
#'   \code{sex_filter}, restricts the computation (e.g. to females).
#' @param sex_filter Value(s) of \code{sex_col} to keep (e.g.
#'   \code{"F"}). Ignored if \code{sex_col} is \code{NULL}.
#' @param max_gsi Upper bound anti-artifact: GSR values above (in \%)
#'   are discarded (default 20, suited to small pelagics / crustaceans;
#'   increase for species with high gonadal investment).
#' @param level Confidence level of the monthly bands (default 0.95).
#' @param title Plot title.
#'
#' @return A \code{ggplot} object. The GSR peak month is annotated.
#' @references King (2007), Fisheries Biology; West (1990),
#'   Aust. J. Mar. Freshwater Res. 41, fecundity estimation methods
#'   and gonadic indices.
#' @seealso \code{\link{spawning_season}}, \code{\link{optimal_closure}}
#' @examples
#' \dontrun{
#' plot_gsi(penaeus_bio, "gonad_weight", "weight", date_col = "Date",
#'          sex_col = "sex", sex_filter = "F")
#' }
#' @export
plot_gsi <- function(data,
                     gonad_col,
                     weight_col,
                     date_col = NULL,
                     sex_col = NULL,
                     sex_filter = NULL,
                     max_gsi = 20,
                     level = 0.95,
                     title = "Gonadosomatic ratio (GSR)") {

  .need_col(data, gonad_col, "gonad_col")
  .need_col(data, weight_col, "weight_col")

  d <- data
  ## optional restriction to a single sex (e.g. females)
  if (!is.null(sex_col) && sex_col %in% names(d) && !is.null(sex_filter)) {
    keep <- toupper(as.character(d[[sex_col]])) %in% toupper(sex_filter)
    d <- d[keep, , drop = FALSE]
  }

  g <- suppressWarnings(as.numeric(d[[gonad_col]]))
  w <- suppressWarnings(as.numeric(d[[weight_col]]))
  ok <- is.finite(g) & is.finite(w) & w > 0 & g > 0
  d <- d[ok, , drop = FALSE]
  d$.gsi <- 100 * g[ok] / w[ok]
  ## anti-artifact filter (data entry outliers)
  d <- d[d$.gsi <= max_gsi, , drop = FALSE]
  if (nrow(d) < 3)
    stop("Too few usable GSR values after filtering.",
         call. = FALSE)

  zc <- stats::qnorm(1 - (1 - level) / 2)
  ci_pct <- round(100 * level)

  ## --- without date: overall distribution ---
  if (is.null(date_col) || !date_col %in% names(d)) {
    m  <- mean(d$.gsi)
    se <- stats::sd(d$.gsi) / sqrt(nrow(d))
    return(
      ggplot2::ggplot(d, ggplot2::aes(.gsi)) +
        ggplot2::geom_histogram(bins = 40, fill = "#7b3294",
                                colour = "white", linewidth = 0.1) +
        ggplot2::geom_vline(xintercept = m, linetype = "dashed",
                            colour = "#c0392b") +
        ggplot2::labs(title = title,
                      subtitle = sprintf("Mean GSR = %.2f %% (CI%d %% : %.2f - %.2f) ; n = %d",
                                         m, ci_pct, m - zc * se, m + zc * se,
                                         nrow(d)),
                      x = "GSR (%)", y = "Count") +
        .theme_fish()
    )
  }

  ## --- with date: monthly mean evolution +/- CI ---
  dt <- as.Date(d[[date_col]])
  d  <- d[!is.na(dt), , drop = FALSE]
  d$.mois <- as.integer(format(as.Date(d[[date_col]]), "%m"))

  agg <- lapply(sort(unique(d$.mois)), function(mo) {
    x <- d$.gsi[d$.mois == mo]; n <- length(x)
    m <- mean(x); se <- if (n > 1) stats::sd(x) / sqrt(n) else NA_real_
    data.frame(mois = mo, n = n, gsi = m,
               lo = m - zc * se, hi = m + zc * se)
  })
  agg <- do.call(rbind, agg)
  agg$mois_lab <- factor(.mois_abbr[agg$mois], levels = .mois_abbr)

  pic <- agg$mois[which.max(agg$gsi)]

  ggplot2::ggplot(agg, ggplot2::aes(mois_lab, gsi, group = 1)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "#7b3294", alpha = 0.2) +
    ggplot2::geom_line(colour = "#7b3294", linewidth = 0.8) +
    ggplot2::geom_point(colour = "#7b3294", size = 2) +
    ggplot2::geom_vline(xintercept = which(.mois_abbr == .mois_abbr[pic]),
                        linetype = "dashed", colour = "#c0392b") +
    ggplot2::annotate("text", x = which(.mois_abbr == .mois_abbr[pic]),
                      y = max(agg$hi, na.rm = TRUE),
                      label = sprintf("peak: %s", .mois_abbr[pic]),
                      hjust = -0.1, size = 3.2, colour = "#c0392b") +
    ggplot2::labs(title = title,
                  subtitle = sprintf("GSR = 100 gonad_weight/weight ; monthly mean (CI%d %%) ; n = %d",
                                     ci_pct, sum(agg$n)),
                  x = "Month", y = "Mean GSR (%)") +
    .theme_fish()

}


#==============================================================
# 9. Catch / CPUE time series
#==============================================================

#' Annual catch or CPUE series by group
#'
#' Plots time series (catches or CPUE) from a table in
#' wide format (one year column + one column per group).
#'
#' @param data data.frame in wide format.
#' @param year_col Year column.
#' @param value_cols Value columns (one per group). If \code{NULL},
#'   all numeric columns other than \code{year_col}.
#' @param ylab Label of the y-axis.
#' @param title Plot title.
#'
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' plot_catch_series(total_catches, "year", ylab = "Catches (t)")
#' }
#' @export
plot_catch_series <- function(data,
                              year_col,
                              value_cols = NULL,
                              ylab = "Value",
                              title = "Annual series") {

  .need_col(data, year_col, "year_col")
  if (is.null(value_cols))
    value_cols <- setdiff(names(data)[vapply(data, is.numeric, logical(1))],
                          year_col)
  if (length(value_cols) == 0)
    stop("No numeric value column.", call. = FALSE)

  long <- do.call(rbind, lapply(value_cols, function(cl) {
    data.frame(annee = data[[year_col]], groupe = cl,
               valeur = data[[cl]], stringsAsFactors = FALSE)
  }))
  long <- long[is.finite(long$valeur), ]

  ggplot2::ggplot(long, ggplot2::aes(annee, valeur, colour = groupe)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_brewer(palette = "Dark2", name = NULL) +
    ggplot2::labs(title = title, x = NULL, y = ylab) +
    .theme_fish()

}


#==============================================================
# S3 plot() methods for objects without visualization
#==============================================================

#' M estimate comparison plot
#'
#' Represents natural mortality estimates by method (points),
#' sorted, with consensus values (geometric mean, median) as reference.
#'
#' @param x \code{FishStockM} object (see \code{\link{estimate_M_all}}).
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @exportS3Method plot FishStockM
plot.FishStockM <- function(x, ...) {
  tab <- x$table
  tab <- tab[is.finite(tab$M) & tab$M > 0, , drop = FALSE]
  tab <- tab[order(tab$M), ]
  tab$method <- factor(tab$method, levels = tab$method)

  cons <- x$consensus
  cons_df <- data.frame(nom = names(cons), val = as.numeric(cons))

  ## Cross-method dispersion band: central 95 percent interval
  ## of the estimates (2.5-97.5 % quantiles). This is not a parametric CI
  ## (each empirical method gives a single point) but a measure of
  ## the uncertainty related to method choice.
  q <- stats::quantile(tab$M, c(0.025, 0.975), names = FALSE)
  band <- data.frame(lo = q[1], hi = q[2])

  ggplot2::ggplot(tab, ggplot2::aes(M, method, colour = source)) +
    ggplot2::geom_rect(data = band,
                       ggplot2::aes(xmin = lo, xmax = hi,
                                    ymin = -Inf, ymax = Inf),
                       inherit.aes = FALSE, fill = "grey80", alpha = 0.4) +
    ggplot2::geom_vline(data = cons_df,
                        ggplot2::aes(xintercept = val, linetype = nom),
                        colour = "grey30") +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_brewer(palette = "Set1", name = "Source") +
    ggplot2::scale_linetype_manual(values = c(geomean = "dashed",
                                              median = "dotted"),
                                   name = "Consensus") +
    ggplot2::labs(title = "Natural mortality M: comparison of methods",
                  subtitle = sprintf("Band = cross-method range 95 %% [%.3f - %.3f]",
                                     q[1], q[2]),
                  x = "M (per year)", y = NULL) +
    .theme_fish()
}


#' Fitted length-weight relationship plot
#'
#' @param x \code{FishStockLW} object (see \code{\link{fit_length_weight}}).
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @exportS3Method plot FishStockLW
plot.FishStockLW <- function(x, ...) {
  d <- x$model$model %||% x$data
  a <- x$a; b <- x$b; r2 <- x$r2
  ## reconstructs the L and W columns from the model if available
  if (!is.null(d) && ncol(d) >= 2) {
    ## the model stores log(W) ~ log(L): back to the natural scale
    Wc <- exp(d[[1]]); Lc <- exp(d[[2]])
    df <- data.frame(L = Lc, W = Wc)
  } else {
    df <- NULL
  }
  b_ci <- x$b_ci
  lab <- if (!is.null(b_ci) && length(b_ci) == 2)
    sprintf("a = %.4g\nb = %.2f [%.2f - %.2f]\nR2 = %.3f",
            a, b, b_ci[1], b_ci[2], r2)
  else sprintf("a = %.4g\nb = %.2f\nR2 = %.3f", a, b, r2)

  Lseq <- if (!is.null(df)) seq(min(df$L), max(df$L), length.out = 100)
          else seq(1, 100, length.out = 100)
  curve <- data.frame(L = Lseq, W = a * Lseq^b)

  ## 95 percent confidence band reconstructed from the stored lm model
  ## (predicted on the log scale then back-transformed).
  if (inherits(x$model, "lm")) {
    pr <- tryCatch(
      stats::predict(x$model, newdata = data.frame(L = Lseq),
                     interval = "confidence", level = 0.95),
      error = function(e) NULL)
    if (!is.null(pr)) {
      curve$lo <- exp(pr[, "lwr"]); curve$hi <- exp(pr[, "upr"])
    }
  }

  p <- ggplot2::ggplot()
  if (!is.null(df))
    p <- p + ggplot2::geom_point(data = df, ggplot2::aes(L, W),
                                 alpha = 0.2, size = 0.7, colour = "#2166ac")
  if (!is.null(curve$lo))
    p <- p + ggplot2::geom_ribbon(data = curve,
                                  ggplot2::aes(x = L, ymin = lo, ymax = hi),
                                  fill = "#c0392b", alpha = 0.2,
                                  inherit.aes = FALSE)
  p +
    ggplot2::geom_line(data = curve, ggplot2::aes(L, W),
                       colour = "#c0392b", linewidth = 0.9) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = lab,
                      hjust = -0.15, vjust = 1.2, size = 3.3) +
    ggplot2::labs(title = "Fitted length-weight relationship (W = a L^b)",
                  subtitle = "95 % confidence band",
                  x = "Length", y = "weight") +
    .theme_fish()
}


#' Mortality breakdown plot (Z, M, F)
#'
#' Stacked bars F / M composing Z, with the exploitation rate E annotated.
#'
#' @param x \code{FishStockMortality} object (see \code{\link{run_mortality}}).
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @exportS3Method plot FishStockMortality
plot.FishStockMortality <- function(x, ...) {
  Z <- x$Z$Z %||% x$Z
  ## consensus M: first the high-level field, otherwise the consensus of
  ## the FishStockM object (geometric mean), otherwise the raw value.
  M <- x$M_consensus %||%
       (if (is.list(x$M)) x$M$M_consensus %||%
          (if (!is.null(x$M$consensus)) x$M$consensus[["geomean"]] else NULL)
        else NULL) %||% x$M
  Z <- as.numeric(Z)[1]; M <- as.numeric(M)[1]
  Fv <- max(Z - M, 0)
  E  <- if (Z > 0) Fv / Z else NA_real_

  df <- data.frame(composante = factor(c("F (fishing)", "M (natural)"),
                                       levels = c("F (fishing)", "M (natural)")),
                   valeur = c(Fv, M))

  ggplot2::ggplot(df, ggplot2::aes(x = "Z", y = valeur, fill = composante)) +
    ggplot2::geom_col(width = 0.5) +
    ggplot2::scale_fill_manual(values = c("F (fishing)" = "#c0392b",
                                          "M (natural)" = "#2980b9"),
                               name = NULL) +
    ggplot2::labs(title = "Total mortality breakdown",
                  subtitle = sprintf("Z = %.3f ; M = %.3f ; F = %.3f ; E = %.2f",
                                     Z, M, Fv, E),
                  x = NULL, y = "Instantaneous rate (per year)") +
    .theme_fish()
}


#' Total allowable catch (TAC) plot
#'
#' Compares current catch, raw TAC and final TAC (after safety margin).
#'
#' @param x \code{FishStockTAC} object (see \code{\link{compute_tac}}).
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @exportS3Method plot FishStockTAC
plot.FishStockTAC <- function(x, ...) {
  catch_cur <- x$inputs$catch_current %||% NA_real_
  df <- data.frame(
    etape = c("Current catch", "Raw TAC", "Final TAC"),
    valeur = c(catch_cur, x$TAC_raw, x$TAC))
  df <- df[is.finite(df$valeur), ]
  df$etape <- factor(df$etape, levels = df$etape)

  ggplot2::ggplot(df, ggplot2::aes(etape, valeur, fill = etape)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = round(valeur, 1)),
                       vjust = -0.4, size = 3.3) +
    ggplot2::scale_fill_manual(values = c("Current catch" = "grey60",
                                          "Raw TAC" = "#f39c12",
                                          "Final TAC" = "#27ae60"),
                               guide = "none") +
    ggplot2::labs(title = "Total allowable catch",
                  subtitle = sprintf("Method: %s ; margin = %.0f %%",
                                     x$method, 100 * x$buffer),
                  x = NULL, y = "Tonnes") +
    .theme_fish()
}


#' Quota allocation by segment plot
#'
#' Bar chart of quotas allocated by segment, sorted.
#'
#' @param x \code{FishStockQuota} object (see \code{\link{allocate_quota}}).
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @exportS3Method plot FishStockQuota
plot.FishStockQuota <- function(x, ...) {
  q <- x$quotas
  q <- q[order(-q$quota_t), ]
  q$segment <- factor(q$segment, levels = q$segment)

  ggplot2::ggplot(q, ggplot2::aes(quota_t, segment)) +
    ggplot2::geom_col(fill = "#2980b9", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f (%.0f%%)",
                                                    quota_t, 100 * prorata)),
                       hjust = -0.05, size = 3) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(title = "Quota allocation by segment",
                  subtitle = sprintf("Total TAC = %.1f t", x$TAC),
                  x = "Quota (t)", y = NULL) +
    .theme_fish()
}


#' Recommended management sizes plot
#'
#' Positions candidate sizes relative to biological references
#' (maturity L50, Lopt, Lc_opt), with the verdict per size shown in colour.
#'
#' @param x \code{FishStockSizeAdvice} object (see \code{\link{recommend_sizes}}).
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @exportS3Method plot FishStockSizeAdvice
plot.FishStockSizeAdvice <- function(x, ...) {
  adv <- x$advice
  ref <- x$reference
  refs <- data.frame(
    nom = c("L50 maturity", "Lopt", "Lc_opt"),
    val = c(ref$L50_maturity, ref$Lopt, ref$Lc_opt))
  refs <- refs[is.finite(refs$val), ]

  ggplot2::ggplot(adv, ggplot2::aes(Lc, y = 1)) +
    ggplot2::geom_vline(data = refs,
                        ggplot2::aes(xintercept = val, colour = nom),
                        linetype = "dashed", linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(fill = verdict), shape = 21,
                        size = 5, colour = "grey20") +
    ggplot2::geom_text(ggplot2::aes(label = Lc), vjust = -1.4, size = 3) +
    ggplot2::scale_colour_brewer(palette = "Dark2", name = "References") +
    ggplot2::scale_y_continuous(limits = c(0.5, 1.5), breaks = NULL) +
    ggplot2::labs(title = "Recommended management sizes",
                  subtitle = "Positioning of candidate sizes (Lc)",
                  x = "First capture size Lc", y = NULL,
                  fill = "Verdict") +
    .theme_fish()
}


## fallback operator (as in ypr.R)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a