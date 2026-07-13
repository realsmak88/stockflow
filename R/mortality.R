###############################################################
#
# stockflow
#
# Module : mortality.R
#
# Estimation of mortality: Z (total), M (natural),
# F (fishing) and E (exploitation rate) and F/M.
#
# Principle: robust wrappers around the OFFICIAL functions
#   - TropFishR::catchCurve()   -> Z
#   - TropFishR::M_empirical()  -> M (multiple methods)
#   - fishmethods::M.empirical() -> M (complementary methods, optional)
# No scientific algorithm is rewritten.
#
###############################################################

#==============================================================
# Dependencies
#==============================================================

.check_tropfish_mort <- function() {

  if (!requireNamespace("TropFishR", quietly = TRUE)) {

    stop(
      "The 'TropFishR' package must be installed.",
      call. = FALSE
    )

  }

  invisible(TRUE)

}


.has_fishmethods <- function() {

  requireNamespace("fishmethods", quietly = TRUE)

}


#==============================================================
# Utility: argument filtering according to the signature
# (multi-version compatibility of packages)
#==============================================================

.match_formals <- function(fun, args) {

  keep <- names(args) %in% names(formals(fun))

  args[keep]

}


#==============================================================
# Utility: robust extraction of Linf and K
#
# Accepts:
#   - a 'FishStockELEFAN' / 'FishStockGrowth' object (via
#     extract_growth_parameters, defined in growth.R);
#   - a list containing $par$Linf and $par$K;
#   - explicit numeric values Linf / K.
#==============================================================

.growth_pars <- function(growth_model = NULL,
                         Linf = NULL,
                         K = NULL) {

  ## 1. Explicit values take priority

  if (!is.null(Linf) && !is.null(K)) {

    return(list(Linf = as.numeric(Linf),
                K    = as.numeric(K)))

  }


  if (is.null(growth_model))
    stop("Provide either 'growth_model', or 'Linf' and 'K'.")


  ## 2. stockflow growth S3 object

  if (inherits(growth_model,
               c("FishStockELEFAN", "FishStockGrowth"))) {

    p <- tryCatch(
      extract_growth_parameters(growth_model),
      error = function(e)
        stop("Unable to extract growth parameters: ",
             conditionMessage(e), call. = FALSE)
    )

    return(list(Linf = as.numeric(p$Linf),
                K    = as.numeric(p$K),
                t_anchor = suppressWarnings(as.numeric(p$t_anchor))))

  }


  ## 3. Generic list with $par

  if (is.list(growth_model) &&
      !is.null(growth_model$par)) {

    return(list(Linf = as.numeric(growth_model$par[["Linf"]]),
                K    = as.numeric(growth_model$par[["K"]])))

  }


  ## 4. Generic list with $result$par

  if (is.list(growth_model) &&
      !is.null(growth_model$result$par)) {

    return(list(Linf = as.numeric(growth_model$result$par[["Linf"]]),
                K    = as.numeric(growth_model$result$par[["K"]])))

  }


  stop("Unable to extract Linf and K from the provided growth model.")

}


#==============================================================
# Total mortality Z: linearized catch curve
#==============================================================

#' Estimation of total mortality Z (catch curve)
#'
#' Robust wrapper around \code{TropFishR::catchCurve()} to estimate the
#' instantaneous total mortality \eqn{Z} from an LFQ object and von Bertalanffy
#' growth parameters.
#'
#' @param lfq Object of class \code{lfq} (TropFishR).
#' @param growth_model stockflow growth object
#'   (\code{FishStockELEFAN} / \code{FishStockGrowth}) or a list containing
#'   \code{$par$Linf} and \code{$par$K}. Ignored if \code{Linf} and \code{K}
#'   are provided directly.
#' @param Linf,K Explicit growth parameters (optional).
#' @param catch_columns Indices of the catch columns to use. Defaults to
#'   all columns.
#' @param reg_int Vector of two integers defining the regression interval
#'   (fully recruited points); required for reproducible, non-interactive
#'   execution. If \code{NULL}, \code{catchCurve} will attempt a selection.
#' @param calc_ogive Logical. Compute the selection ogive (probability of
#'   capture)? Passed to \code{TropFishR::catchCurve()}.
#' @param plot Logical. Plot the catch curve.
#' @param M Natural mortality (optional) passed to the automatic
#'   regression interval selection (\code{\link{auto_reg_int}})
#'   as a safeguard when \code{reg_int} is not provided.
#' @param ... Additional arguments passed to
#'   \code{TropFishR::catchCurve()}.
#'
#' @return S3 object \code{FishStockZ}: a list containing \code{Z}, the
#'   confidence interval \code{Z_ci}, the complete \code{catchCurve} object and
#'   the growth parameters used.
#'
#' @references
#'   Sparre, P., Venema, S.C. (1998). Introduction to tropical fish stock
#'   assessment. FAO Fisheries Technical Paper 306/1.
#'   Mildenberger, T.K., Taylor, M.H., Wolff, M. (2017). TropFishR. MEE 8(11).
#'
#' @examples
#' \dontrun{
#'   z <- estimate_catchcurve(lfq, ga, reg_int = c(9, 21))
#'   z$Z
#' }
#' @export
estimate_catchcurve <- function(lfq,
                                growth_model = NULL,
                                Linf = NULL,
                                K = NULL,
                                catch_columns = NULL,
                                reg_int = NULL,
                                calc_ogive = TRUE,
                                plot = FALSE,
                                M = NULL,
                                ...) {

  .check_tropfish_mort()

  if (missing(lfq) || !inherits(lfq, "lfq"))
    stop("'lfq' must be an object of class 'lfq'.")


  gp <- .growth_pars(growth_model, Linf, K)


  if (is.null(catch_columns)) {

    catch_columns <- seq_len(ncol(lfq$catch))

  }


  ## reg_int not provided -> automatic selection (avoids TropFishR's
  ## interactive call, which fails in a non-interactive context).
  auto_used <- FALSE

  if (is.null(reg_int)) {

    reg_int <- auto_reg_int(lfq, Linf = gp$Linf, K = gp$K,
                            catch_columns = catch_columns, M = M)

    auto_used <- TRUE

    if (is.null(reg_int))
      stop("Automatic selection of reg_int impossible: no interval ",
           "gives a finite and positive Z. Inspect the catch curve ",
           "and provide reg_int explicitly.", call. = FALSE)

  }


  ## IMPORTANT: TropFishR::catchCurve() does NOT have 'Linf' / 'K' arguments.
  ## It looks for them INSIDE the object passed as 'param'. Passing them as
  ## named arguments has no effect (they are filtered out by the formals), and
  ## the function then fails with:
  ##   "You need to assign values to Linf and K for the catch curve ..."
  ## We therefore attach them to the lfq object.

  lfq$Linf <- gp$Linf
  lfq$K    <- gp$K

  if (!is.null(gp$t_anchor) && is.finite(gp$t_anchor))
    lfq$t_anchor <- gp$t_anchor


  args <- list(
    param         = lfq,
    catch_columns = catch_columns,
    reg_int       = reg_int,
    calc_ogive    = calc_ogive,
    plot          = plot,
    ...
  )

  ## The name of the 1st argument has changed across versions (param / x).
  args_p <- .match_formals(TropFishR::catchCurve, args)

  if (!("param" %in% names(args_p)) &&
      "x" %in% names(formals(TropFishR::catchCurve))) {

    args_p <- c(list(x = lfq), args_p)
    args_p <- .match_formals(TropFishR::catchCurve, args_p)

  }


  cc <- tryCatch(

    do.call(TropFishR::catchCurve, args_p),

    error = function(e)
      stop("catchCurve: ", conditionMessage(e), call. = FALSE)

  )


  ## Robust extraction of Z and its CI

  Z <- .extract_Z(cc)

  Z_ci <-
    if (!is.null(cc$confidenceInt))
      as.numeric(cc$confidenceInt)
  else
    c(NA_real_, NA_real_)


  structure(

    list(
      Z            = Z,
      Z_ci         = Z_ci,
      reg_int      = reg_int,
      reg_int_auto = auto_used,
      catchCurve   = cc,
      growth       = gp
    ),

    class = "FishStockZ"

  )

}


## Historical alias
#' @rdname estimate_catchcurve
#' @export
estimate_Z <- estimate_catchcurve


.extract_Z <- function(cc) {

  ## depending on version: cc$Z may be named or not
  if (!is.null(cc$Z)) {

    z <- cc$Z

    if (!is.null(names(z)) &&
        any(grepl("Z", names(z), ignore.case = TRUE))) {

      return(as.numeric(z[grep("Z", names(z), ignore.case = TRUE)][1]))

    }

    return(as.numeric(z[1]))

  }

  NA_real_

}


#' @export
print.FishStockZ <- function(x, ...) {

  cat("\n")
  cat("--------------------------------------\n")
  cat("Total mortality Z (catch curve)\n")
  cat("--------------------------------------\n\n")

  cat(sprintf("Z        : %.4f\n", x$Z))

  if (all(is.finite(x$Z_ci)))
    cat(sprintf("95%% CI   : [%.4f ; %.4f]\n",
                x$Z_ci[1], x$Z_ci[2]))

  cat(sprintf("Linf     : %.3f\n", x$growth$Linf))
  cat(sprintf("K        : %.4f\n", x$growth$K))

  invisible(x)

}


#==============================================================
# Automatic selection of the regression interval (reg_int)
#==============================================================

#' Automatic selection of the catch curve regression interval
#'
#' \code{TropFishR::catchCurve()} normally requires manually designating with
#' the mouse the two bounds of the descending, fully recruited part of the
#' linearized catch curve. This function automates that choice: it explores a
#' grid of candidate intervals, fits the curve with the official function for
#' each, and keeps the one giving the best compromise between number of
#' points and biological plausibility.
#'
#' An interval is kept only if \eqn{Z} is finite and strictly positive.
#' If \code{M} is provided, intervals leading to \eqn{Z \le M} (hence a
#' negative fishing mortality, biologically impossible) are discarded.
#'
#' @param lfq \code{lfq} object (TropFishR).
#' @param Linf,K Growth parameters.
#' @param catch_columns Catch columns to use.
#' @param M Natural mortality (optional): serves as a safeguard, a
#'   \eqn{Z \le M} being rejected.
#' @param min_points Minimum number of points in the regression.
#'
#' @return Vector of two integers usable as \code{reg_int}, or
#'   \code{NULL} if no plausible interval was found.
#'
#' @section Warning:
#'   This selection automates a visual judgment; it does not replace it.
#'   Check the produced catch curve (\code{plot = TRUE}) before
#'   releasing a \eqn{Z}.
#'
#' @examples
#' \dontrun{
#'   ri <- auto_reg_int(lfq, Linf = 45, K = 1.2, M = 1.5)
#'   z  <- estimate_catchcurve(lfq, Linf = 45, K = 1.2, reg_int = ri)
#' }
#' @export
auto_reg_int <- function(lfq,
                         Linf,
                         K,
                         catch_columns = NULL,
                         M = NULL,
                         min_points = 4) {

  .check_tropfish_mort()

  if (is.null(catch_columns))
    catch_columns <- seq_len(ncol(lfq$catch))

  Ct <- rowSums(as.matrix(lfq$catch)[, catch_columns, drop = FALSE],
                na.rm = TRUE)
  L  <- lfq$midLengths

  ## usable points: below Linf and with catches
  usable <- which(L < 0.98 * Linf & Ct > 0)

  if (length(usable) < min_points + 1)
    return(NULL)

  peak <- usable[which.max(Ct[usable])]     # mode = full recruitment
  last <- max(usable)

  ## lfq object carrying Linf/K (catchCurve reads them from the object)
  lfq2 <- lfq
  lfq2$Linf <- Linf
  lfq2$K    <- K

  ## grid of candidates: start just after the mode, end near the last
  ## usable point (the tail is often noisy due to low counts)
  starts <- seq(peak, min(peak + 3, last - min_points))
  ends   <- seq(max(last - 4, min(starts) + min_points), last)

  best <- NULL

  for (s in starts) {
    for (e in ends) {

      if (e - s + 1 < min_points) next

      cc <- tryCatch(
        suppressWarnings(suppressMessages(
          TropFishR::catchCurve(lfq2, catch_columns = catch_columns,
                                reg_int = c(s, e), calc_ogive = FALSE,
                                plot = FALSE)
        )),
        error = function(err) NULL
      )

      if (is.null(cc)) next

      Z <- .extract_Z(cc)

      if (!is.finite(Z) || Z <= 0) next
      if (!is.null(M) && is.finite(M) && Z <= M) next   # negative F: rejected

      npts  <- e - s + 1
      score <- npts        # at equal plausibility, we prefer more points

      if (is.null(best) || score > best$score)
        best <- list(reg_int = c(s, e), Z = Z, score = score)

    }
  }

  if (is.null(best)) return(NULL)

  best$reg_int

}


#==============================================================
# Total mortality Z via FSA (age-structured data)
#==============================================================

#' Estimation of Z by Chapman-Robson and catch curve (FSA)
#'
#' Wrappers around \code{FSA::chapmanRobson()} and \code{FSA::catchCurve()}
#' to estimate \eqn{Z} from an age composition (fully recruited ages).
#' Complementary to \code{estimate_catchcurve()} (length-based approach,
#' TropFishR).
#'
#' @param age Vector of ages (age classes).
#' @param catch Vector of catches (or frequencies) per age.
#' @param ages2use Ages to use (descending, fully recruited part).
#'   If \code{NULL}, all ages from the modal age onward are used.
#' @return S3 object \code{FishStockZ} containing the
#'   Chapman-Robson and catch curve estimates (\code{$methods}).
#' @references
#'   Chapman, D.G., Robson, D.S. (1960). Biometrics 16.
#'   Ogle, D.H. (2016). Introductory Fisheries Analyses with R. FSA.
#' @examples
#' \dontrun{
#'   z <- estimate_Z_fsa(age = 1:10, catch = c(5,40,60,55,40,25,15,8,4,2))
#'   z$Z
#' }
#' @export
estimate_Z_fsa <- function(age, catch, ages2use = NULL) {

  if (!requireNamespace("FSA", quietly = TRUE))
    stop("The 'FSA' package must be installed.", call. = FALSE)

  if (length(age) != length(catch))
    stop("'age' and 'catch' must have the same length.")

  df <- data.frame(age = age, catch = catch)

  if (is.null(ages2use)) {

    peak <- df$age[which.max(df$catch)]
    ages2use <- df$age[df$age >= peak]

  }


  ## Chapman-Robson
  cr <- tryCatch({

    obj <- FSA::chapmanRobson(catch ~ age, data = df, ages2use = ages2use)
    s   <- summary(obj)
    list(obj = obj,
         Z = as.numeric(s[grep("^Z$", rownames(s)), "Estimate"][1]))

  }, error = function(e)
    list(obj = NULL, Z = NA_real_))


  ## Catch curve (linearized regression)
  cc <- tryCatch({

    obj <- FSA::catchCurve(catch ~ age, data = df, ages2use = ages2use)
    s   <- summary(obj)
    list(obj = obj,
         Z = as.numeric(s[grep("^Z$", rownames(s)), "Estimate"][1]))

  }, error = function(e)
    list(obj = NULL, Z = NA_real_))


  Z <- if (is.finite(cr$Z)) cr$Z else cc$Z


  structure(

    list(
      Z       = Z,
      Z_ci    = c(NA_real_, NA_real_),
      methods = list(ChapmanRobson = cr$Z, CatchCurve = cc$Z),
      objects = list(ChapmanRobson = cr$obj, CatchCurve = cc$obj),
      growth  = NULL
    ),

    class = "FishStockZ"

  )

}


#==============================================================
# Natural mortality M: multiple empirical methods
#==============================================================

## Default longevity (tmax) from K: proxy Amax ~ 3/K
## (age at which ~95% of Linf is reached). Serves as INPUT to
## methods based on tmax; is not an estimator of M.

.default_tmax <- function(K) {

  if (is.null(K) || !is.finite(K) || K <= 0)
    return(NULL)

  3 / K

}


#' Estimation of natural mortality M by multiple methods
#'
#' Computes \eqn{M} using all applicable empirical methods via the
#' official functions \code{TropFishR::M_empirical()} and, if available,
#' \code{fishmethods::M.empirical()}, plus the Jensen (1996) estimator
#' based on K. Returns a comparative table and a consensus.
#'
#' @param Linf,K von Bertalanffy growth parameters.
#' @param temp Mean water temperature (degrees C) for Pauly's method.
#' @param tmax Maximum age (longevity). If \code{NULL}, a proxy \eqn{3/K} is
#'   used for the methods that require it (Hoenig, Then_tmax,
#'   AlversonCarney).
#' @param tm50 Age at 50 percent maturity (for Rikhter-Efanov / Roff /
#'   Jensen based on maturity). Optional.
#' @param t0 Theoretical age at zero length (von Bertalanffy). Optional;
#'   used by some \code{FSA::metaM()} methods.
#' @param Winf Asymptotic weight (weight-based methods). Optional.
#' @param Lmean Mean (or representative) length of the stock. Optional;
#'   required by size-dependent methods (Gislason, Charnov).
#' @param methods Vector of methods or \code{"auto"} (all applicable methods
#'   according to the available inputs).
#' @param use_fsa Logical. Include the ~30 estimators from
#'   \code{FSA::metaM()} (see \code{FSA::Mmethods()}) in addition to those
#'   from \code{TropFishR::M_empirical()}.
#' @param consensus Summary function(s): \code{"geomean"} and/or
#'   \code{"median"}.
#'
#' @return S3 object \code{FishStockM}: a list with \code{table} (data.frame
#'   method / M / source), \code{consensus} (named numeric) and the inputs.
#'
#' @references
#'   Pauly, D. (1980). ICES J. Mar. Sci. 39(2).
#'   Then, A.Y. et al. (2015). ICES J. Mar. Sci. 72(1).
#'   Hoenig, J.M. (1983). Fish. Bull. 82.
#'   Jensen, A.L. (1996). Can. J. Fish. Aquat. Sci. 53.
#'   Gislason, H. et al. (2010). Fish and Fisheries 11.
#'   Lorenzen, K. (1996). J. Fish Biol. 49.
#'
#' @examples
#' \dontrun{
#'   m <- estimate_M_all(Linf = 86.4, K = 0.18, temp = 25)
#'   m$table
#'   m$consensus
#' }
#' @export
estimate_M_all <- function(Linf,
                           K,
                           temp = 25,
                           tmax = NULL,
                           tm50 = NULL,
                           t0 = NULL,
                           Winf = NULL,
                           Lmean = NULL,
                           methods = "auto",
                           use_fsa = TRUE,
                           consensus = c("geomean", "median")) {

  .check_tropfish_mort()

  if (missing(Linf) || missing(K))
    stop("'Linf' and 'K' are required.")


  if (is.null(tmax))
    tmax <- .default_tmax(K)


  ## ---- catalog of TropFishR::M_empirical methods --------------
  ## each entry: required arguments present?

  inputs <- list(
    Linf = Linf, Winf = Winf, K_l = K, K_w = K,
    temp = temp, tmax = tmax, tm50 = tm50
  )

  tf_methods <- c(
    "Pauly_Linf",     # Linf, K_l, temp
    "Pauly_Winf",     # Winf, K_w, temp
    "Hoenig",         # tmax
    "Then_growth",    # Linf, K_l
    "Then_tmax",      # tmax
    "AlversonCarney", # K_l, tmax
    "RikhterEfanov",  # tm50
    "Roff",           # K_l, tm50
    "Gislason"        # Linf, K_l (+ Bl)
  )

  needs <- list(
    Pauly_Linf     = c("Linf", "K_l", "temp"),
    Pauly_Winf     = c("Winf", "K_w", "temp"),
    Hoenig         = c("tmax"),
    Then_growth    = c("Linf", "K_l"),
    Then_tmax      = c("tmax"),
    AlversonCarney = c("K_l", "tmax"),
    RikhterEfanov  = c("tm50"),
    Roff           = c("K_l", "tm50"),
    Gislason       = c("Linf", "K_l")
  )

  if (!identical(methods, "auto"))
    tf_methods <- intersect(tf_methods, methods)


  rows <- list()

  for (meth in tf_methods) {

    req <- needs[[meth]]

    have <- all(vapply(req,
                       function(a) !is.null(inputs[[a]]) &&
                                   is.finite(inputs[[a]]),
                       logical(1)))

    if (!have) next


    val <- tryCatch({

      args <- .match_formals(
        TropFishR::M_empirical,
        c(inputs, list(method = meth))
      )

      out <- do.call(TropFishR::M_empirical, args)

      ## M_empirical returns a named matrix/value
      as.numeric(out)[1]

    },
    error = function(e) NA_real_)


    if (is.finite(val))
      rows[[length(rows) + 1]] <-
        data.frame(method = meth, M = val,
                   source = "TropFishR",
                   stringsAsFactors = FALSE)

  }


  ## ---- Jensen (1996) based on K (published identity) ----------------

  if (identical(methods, "auto") || "Jensen_K" %in% methods)
    rows[[length(rows) + 1]] <-
      data.frame(method = "Jensen_K", M = 1.5 * K,
                 source = "Jensen1996", stringsAsFactors = FALSE)

  if ((identical(methods, "auto") || "Jensen_tm50" %in% methods) &&
      !is.null(tm50) && is.finite(tm50))
    rows[[length(rows) + 1]] <-
      data.frame(method = "Jensen_tm50", M = 1.65 / tm50,
                 source = "Jensen1996", stringsAsFactors = FALSE)


  ## ---- complementary methods via fishmethods (optional) -------

  if (.has_fishmethods()) {

    fm <- tryCatch({

      args <- .match_formals(
        fishmethods::M.empirical,
        list(Linf = Linf, Kl = K, `T` = temp,
             tmax = tmax, tm = tm50, method = 1)
      )

      do.call(fishmethods::M.empirical, args)

    },
    error = function(e) NULL)

    if (!is.null(fm)) {

      fm <- as.data.frame(fm)

      rows[[length(rows) + 1]] <-
        data.frame(
          method = paste0("fishmethods_", rownames(fm)),
          M      = as.numeric(fm[[1]]),
          source = "fishmethods",
          stringsAsFactors = FALSE
        )

    }

  }


  ## ---- FSA empirical methods (FSA::metaM on FSA::Mmethods) --
  ## Covers ~30 estimators (Hoenig*, Pauly*, Jensen*, Gislason,
  ## AlversonCarney, Charnov, ZhangMegrey, RikhterEfanov, ChenWatanabe,
  ## PetersonWroblewski, etc.). Each method is called via the OFFICIAL
  ## function; only those whose inputs are available succeed.

  if (isTRUE(use_fsa) &&
      requireNamespace("FSA", quietly = TRUE)) {

    fsa_methods <-
      if (identical(methods, "auto"))
        FSA::Mmethods("all")
      else
        intersect(FSA::Mmethods("all"), methods)

    ## Arguments recognized by FSA::metaM() (see ?FSA::metaM).
    ## Each method uses only a subset of them; superfluous
    ## inputs are simply ignored by metaM.
    fsa_args <- list(
      Linf = Linf, Winf = Winf, K = K, t0 = t0, b = 3,
      tmax = tmax, t50 = tm50, Temp = temp, L = Lmean
    )
    fsa_args <- fsa_args[!vapply(fsa_args, is.null, logical(1))]

    for (meth in fsa_methods) {

      ## metaM() is called method by method: this way a method whose
      ## inputs are missing (e.g. QuinnDeriso) does not stop the
      ## others. The output is a data.frame (M column) depending on the version.
      val <- tryCatch({

        out <- do.call(FSA::metaM,
                       c(list(meth, verbose = FALSE), fsa_args))
        if (is.data.frame(out)) out$M[1] else as.numeric(out)[1]

      },
        error   = function(e) NA_real_,
        warning = function(w) NA_real_
      )

      if (is.finite(val) && val > 0)
        rows[[length(rows) + 1]] <-
          data.frame(method = paste0("FSA_", meth), M = val,
                     source = "FSA", stringsAsFactors = FALSE)

    }

  }


  if (length(rows) == 0)
    stop("No M method applicable with the provided inputs.")


  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL


  ## ---- consensus --------------------------------------------------

  cons <- numeric(0)
  Mv <- tab$M[is.finite(tab$M) & tab$M > 0]

  if ("geomean" %in% consensus)
    cons["geomean"] <- exp(mean(log(Mv)))

  if ("median" %in% consensus)
    cons["median"] <- stats::median(Mv)


  structure(

    list(
      table     = tab,
      consensus = cons,
      inputs    = list(Linf = Linf, K = K, temp = temp,
                       tmax = tmax, tm50 = tm50, Winf = Winf)
    ),

    class = "FishStockM"

  )

}


## Single-method wrapper (backward compatibility)
#' Estimation of M by a single method
#' @inheritParams estimate_M_all
#' @param method Name of a method (default \code{"Then_growth"}, recommended
#'   by Then et al. 2015 in the absence of age data).
#' @return Numeric value of M.
#' @export
estimate_M <- function(Linf,
                       K,
                       method = "Then_growth",
                       temp = 25,
                       tmax = NULL,
                       tm50 = NULL,
                       Winf = NULL) {

  res <- estimate_M_all(
    Linf = Linf, K = K, temp = temp,
    tmax = tmax, tm50 = tm50, Winf = Winf,
    methods = method
  )

  res$table$M[match(method, res$table$method)]

}


#' @export
print.FishStockM <- function(x, ...) {

  cat("\n")
  cat("--------------------------------------\n")
  cat("Natural mortality M (multiple methods)\n")
  cat("--------------------------------------\n\n")

  print(x$table, row.names = FALSE)

  cat("\nConsensus:\n")

  for (nm in names(x$consensus))
    cat(sprintf("  %-8s : %.4f\n", nm, x$consensus[[nm]]))

  invisible(x)

}


#==============================================================
# F, E, F/M and exploitation diagnosis
#==============================================================

#' Fishing mortality F, exploitation rate E and F/M ratio
#'
#' @param Z Total mortality.
#' @param M Natural mortality (single value, e.g. consensus).
#' @param E_target Exploitation reference (default 0.5; Gulland 1971:
#'   \eqn{E_{opt}\approx0.5}).
#' @return data.frame: Z, M, F, E, FM (=F/M) and exploitation status.
#' @references Gulland, J.A. (1971). The fish resources of the ocean.
#' @export
estimate_FE <- function(Z, M, E_target = 0.5) {

  Fmort <- Z - M
  E     <- Fmort / Z
  FM    <- Fmort / M

  statut <- if (is.finite(E)) {

    if (E > E_target + 0.05) "Overexploitation (E > 0.5)"
    else if (E < E_target - 0.05) "Underexploitation (E < 0.5)"
    else "Exploitation close to optimum (E ~ 0.5)"

  } else NA_character_


  data.frame(
    Z = Z, M = M, F = Fmort, E = E, FM = FM,
    statut = statut,
    stringsAsFactors = FALSE
  )

}


## Historical alias
#' @rdname estimate_FE
#' @export
calculate_F_E <- function(Z, M) estimate_FE(Z, M)[, c("Z", "M", "F", "E")]


#==============================================================
# Comparison of M methods (table + associated F/E)
#==============================================================

#' Comparison of M methods with associated F, E and F/M
#'
#' @param Z Total mortality (from \code{estimate_catchcurve}).
#' @param M_obj \code{FishStockM} object (from \code{estimate_M_all}).
#' @param E_target Exploitation reference.
#' @return Sorted data.frame, one row per M method.
#' @export
compare_M_methods <- function(Z, M_obj, E_target = 0.5) {

  if (!inherits(M_obj, "FishStockM"))
    stop("'M_obj' must be a 'FishStockM' object.")

  tab <- M_obj$table

  fe <- do.call(rbind, lapply(tab$M, function(m)
    estimate_FE(Z, m, E_target)))

  out <- cbind(tab, fe[, c("F", "E", "FM", "statut")])

  out <- out[order(out$M), ]
  rownames(out) <- NULL

  out

}


#==============================================================
# Orchestrator: complete mortality workflow
#==============================================================

#' Complete mortality workflow (Z, M, F, E)
#'
#' Chains the estimation of Z (catch curve), the computation of M by multiple
#' methods, and the derivation of F, E and F/M for each method, with a
#' consensus.
#'
#' @param lfq \code{lfq} object.
#' @param growth_model stockflow growth object (or list with parameters).
#' @param temp Temperature (Pauly's method).
#' @param tmax,tm50,t0,Winf,Lmean Optional inputs passed to
#'   \code{\link{estimate_M_all}} (longevity, age at maturity, t0, asymptotic
#'   weight, mean length).
#' @param reg_int Regression interval for the catch curve
#'   (reproducibility).
#' @param catch_columns Catch columns.
#' @param E_target Exploitation reference (default 0.5).
#' @param plot Plot the catch curve.
#' @param ... Passed to \code{estimate_catchcurve()}.
#'
#' @return S3 object \code{FishStockMortality}.
#' @examples
#' \dontrun{
#'   mort <- run_mortality(lfq, ga, temp = 25, reg_int = c(9, 21))
#'   summary(mort)
#' }
#' @export
run_mortality <- function(lfq,
                          growth_model,
                          temp = 25,
                          tmax = NULL,
                          tm50 = NULL,
                          t0 = NULL,
                          Winf = NULL,
                          Lmean = NULL,
                          reg_int = NULL,
                          catch_columns = NULL,
                          E_target = 0.5,
                          plot = FALSE,
                          ...) {

  ## M is estimated BEFORE Z: it serves as a safeguard for the automatic
  ## selection of reg_int (an interval giving Z <= M would imply F < 0).
  gp0 <- .growth_pars(growth_model)

  message("Estimating M (multiple methods)...")

  m_obj <- estimate_M_all(
    Linf = gp0$Linf, K = gp0$K, temp = temp,
    tmax = tmax, tm50 = tm50, t0 = t0,
    Winf = Winf, Lmean = Lmean
  )

  M_star <- if ("geomean" %in% names(m_obj$consensus))
    m_obj$consensus[["geomean"]] else stats::median(m_obj$table$M)


  message("Estimating Z (catch curve)...")

  z_obj <- estimate_catchcurve(
    lfq = lfq,
    growth_model = growth_model,
    catch_columns = catch_columns,
    reg_int = reg_int,
    plot = plot,
    M = M_star,
    ...
  )


  gp <- z_obj$growth

  if (is.finite(z_obj$Z) && is.finite(M_star) && z_obj$Z <= M_star)
    warning("Z (", round(z_obj$Z, 3), ") is less than or equal to M (",
            round(M_star, 3), "): fishing mortality would be negative. ",
            "The regression interval (reg_int) is very likely poorly ",
            "placed. Inspect the catch curve.", call. = FALSE)


  message("Computing F, E and F/M...")

  comparison <- compare_M_methods(z_obj$Z, m_obj, E_target)

  ## M_star has already been computed above (it serves as a safeguard for reg_int)
  summary_fe <- estimate_FE(z_obj$Z, M_star, E_target)


  structure(

    list(
      Z          = z_obj,
      M          = m_obj,
      comparison = comparison,
      M_consensus = M_star,
      summary    = summary_fe,
      growth     = gp
    ),

    class = "FishStockMortality"

  )

}


#' @export
print.FishStockMortality <- function(x, ...) {

  cat("\n=====================================\n")
  cat(" stockflow : Mortality\n")
  cat("=====================================\n\n")

  cat(sprintf("Z         : %.4f\n", x$Z$Z))
  cat(sprintf("M (cons.) : %.4f\n", x$M_consensus))
  cat(sprintf("F         : %.4f\n", x$summary$F))
  cat(sprintf("E = F/Z   : %.4f\n", x$summary$E))
  cat(sprintf("F/M       : %.4f\n", x$summary$FM))
  cat(sprintf("Status    : %s\n", x$summary$statut))

  invisible(x)

}


#' @export
summary.FishStockMortality <- function(object, ...) {

  print(object)

  cat("\n-- Detail by M method --\n\n")

  print(object$comparison, row.names = FALSE)

  invisible(object)

}
