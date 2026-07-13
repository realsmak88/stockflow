###############################################################
#
# stockflow
#
# Module : growth.R
#
# Estimation of growth parameters
# Based exclusively on official functions
# of the TropFishR package
#
###############################################################

#==============================================================
# Dependency check
#==============================================================

check_tropfish <- function() {
  
  if (!requireNamespace("TropFishR", quietly = TRUE)) {
    
    stop(
      "The 'TropFishR' package must be installed.",
      call. = FALSE
    )
    
  }
  
  invisible(TRUE)
  
}

#==============================================================
# LFQ object verification
#==============================================================

#' Verify the validity of an LFQ object
#'
#' Checks that an object is indeed of class \code{lfq} (TropFishR) and that it
#' contains the required fields (\code{dates}, \code{midLengths}, \code{catch}).
#'
#' @param lfq Object to verify.
#' @return \code{TRUE} invisibly if the object is valid ; raises an
#'   error otherwise.
#' @seealso \code{\link{prepare_tropfish}}, \code{\link{summary_lfq}}
#' @examples
#' \dontrun{
#'   check_lfq(lfq)
#' }
#' @export
check_lfq <- function(lfq) {
  
  check_tropfish()
  
  if (missing(lfq))
    stop("Missing LFQ object.")
  
  if (!inherits(lfq, "lfq"))
    stop("The object must be of class 'lfq'.")
  
  required <- c(
    "dates",
    "midLengths",
    "catch"
  )
  
  miss <- setdiff(required, names(lfq))
  
  if (length(miss) > 0) {
    
    stop(
      "Invalid LFQ object. Missing fields : ",
      paste(miss, collapse = ", ")
    )
    
  }
  
  invisible(TRUE)
  
}

#==============================================================
# LFQ summary
#==============================================================

#' Summary of an LFQ object
#'
#' Displays the main characteristics of an \code{lfq} object : number of
#' sampling dates, number of length classes, total count and
#' length range.
#'
#' @param lfq Object of class \code{lfq} (TropFishR).
#' @return The \code{lfq} object, invisibly (called for its display
#'   effect).
#' @seealso \code{\link{check_lfq}}, \code{\link{plot_lfq}}
#' @examples
#' \dontrun{
#'   summary_lfq(lfq)
#' }
#' @export
summary_lfq <- function(lfq) {
  
  check_lfq(lfq)
  
  cat("\n")
  cat("====================================\n")
  cat(" stockflow - LFQ summary\n")
  cat("====================================\n\n")
  
  cat("Number of dates :", length(lfq$dates), "\n")
  cat("Number of classes :", length(lfq$midLengths), "\n")
  cat("Total count :", sum(lfq$catch), "\n")
  cat("Min length :", min(lfq$midLengths), "\n")
  cat("Max length :", max(lfq$midLengths), "\n")
  
  invisible(lfq)
  
}

#==============================================================
# Plot LFQ
#==============================================================

#' Plot an LFQ object
#'
#' Plots the length frequencies, with optional restructuring via
#' \code{TropFishR::lfqRestructure()} (highlighting of cohorts).
#'
#' @param lfq Object of class \code{lfq} (TropFishR).
#' @param restructure Logical. Restructure the object before plotting ?
#' @param MA Width of the moving average used by the restructuring.
#' @param ... Additional arguments passed to the plotting method.
#' @return The plotted object, invisibly.
#' @references Mildenberger, T.K., Taylor, M.H., Wolff, M. (2017). TropFishR.
#'   \emph{Methods in Ecology and Evolution} 8(11).
#' @seealso \code{\link{summary_lfq}}
#' @examples
#' \dontrun{
#'   plot_lfq(lfq, restructure = TRUE, MA = 5)
#' }
#' @export
plot_lfq <- function(lfq,
                     restructure = FALSE,
                     MA = 5,
                     ...) {
  
  check_lfq(lfq)
  
  obj <- lfq
  
  if (restructure) {
    
    obj <- TropFishR::lfqRestructure(
      obj,
      MA = MA
    )
    
  }
  
  plot(obj, ...)
  
  invisible(obj)
  
}

#==============================================================
# FishStockGrowth class
#==============================================================

new_growth <- function(results,
                       lfq,
                       metadata = list()) {
  
  structure(
    
    list(
      
      results = results,
      
      lfq = lfq,
      
      metadata = metadata,
      
      date = Sys.time()
      
    ),
    
    class = "FishStockGrowth"
    
  )
  
}

#==============================================================
# print()
#==============================================================

#' @export
print.FishStockGrowth <- function(x, ...) {
  
  cat("\n")
  cat("====================================\n")
  cat(" FishStockGrowth object\n")
  cat("====================================\n\n")
  
  cat("Created on :",
      format(x$date),
      "\n\n")
  
  if (!is.null(x$results$best_model)) {
    
    cat("Best model :",
        x$results$best_model,
        "\n")
    
  }
  
  invisible(x)
  
}

#==============================================================
# summary()
#==============================================================

#' @export
summary.FishStockGrowth <- function(object,
                                    ...) {
  
  print(object)
  
  if (!is.null(object$results$comparison)) {
    
    print(object$results$comparison)
    
  }
  
  invisible(object)
  
}

#==============================================================
# Plot FishStockGrowth object
#==============================================================

#' @export
plot.FishStockGrowth <- function(x,
                                 type = c(
                                   "lfq",
                                   "growth"
                                 ),
                                 ...) {
  
  type <- match.arg(type)
  
  if (type == "lfq") {
    
    plot_lfq(x$lfq)
    
  }
  
  invisible(x)
  
}

#==============================================================
# ELEFAN verification
#==============================================================

check_elefan <- function(model) {
  
  if (is.null(model))
    stop("NULL ELEFAN object.")
  
  if (is.null(model$par))
    stop("ELEFAN parameters missing.")
  
  if (!all(c("Linf", "K") %in% names(model$par))) {
    
    stop(
      "Linf or K missing."
    )
    
  }
  
  invisible(TRUE)
  
}

#==============================================================
# Linf extraction
#==============================================================

get_Linf <- function(model) {
  
  check_elefan(model)
  
  model$par$Linf
  
}

#==============================================================
# K extraction
#==============================================================

get_K <- function(model) {
  
  check_elefan(model)
  
  model$par$K
  
}

#==============================================================
# t_anchor extraction
#==============================================================

get_tanchor <- function(model) {
  
  check_elefan(model)
  
  if ("t_anchor" %in% names(model$par)) {
    
    return(model$par$t_anchor)
    
  }
  
  NA
  
}

#==============================================================
# LFQ object preparation
#==============================================================

prepare_growth_lfq <- function(lfq,
                               MA = 5) {
  
  check_lfq(lfq)
  
  if (!"rcounts" %in% names(lfq)) {
    
    lfq <- TropFishR::lfqRestructure(
      
      lfq,
      
      MA = MA
      
    )
    
  }
  
  lfq
  
}

###############################################################
#
# Block 2A
#
# Powell-Wetherall
# ELEFAN
# Official TropFishR wrappers
#
###############################################################

#==============================================================
# Execution time measurement
#==============================================================

.time_it <- function(expr){

  start <- Sys.time()

  ## Some TropFishR model-fitting functions (ELEFAN, powell_wetherall, ...)
  ## draw a diagnostic plot even when plotting is not requested. On a small
  ## graphics device (e.g. the RStudio plot pane during devtools::test())
  ## this raises "figure margins too large" and aborts the fit. Redirect any
  ## such forced plotting to a null PDF device with ample margins so the fit
  ## always succeeds, regardless of the active device.
  grDevices::pdf(file = NULL, width = 7, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)

  value <- force(expr)

  elapsed <- Sys.time() - start

  list(
    value=value,
    elapsed=elapsed
  )

}

#==============================================================
# Powell result verification
#==============================================================

check_powell <- function(model){

  if(is.null(model))
    stop(
      "Powell-Wetherall returned NULL. ",
      "TropFishR::powell_wetherall() requires an INTERACTIVE selection of ",
      "the regression interval when 'reg_int' is not provided : in a ",
      "non-interactive context (knitr, Rscript, R CMD check) it fails ",
      "silently. Explicitly provide reg_int, for example : ",
      "run_powell(lfq, reg_int = c(10, 25)).",
      call. = FALSE
    )

  invisible(TRUE)

}

#==============================================================
# Powell-Wetherall
#==============================================================

#' Linf diagnostic via the Powell-Wetherall method
#'
#' Wrapper around \code{TropFishR::powell_wetherall()}.
#'
#' @section Important:
#'   \code{TropFishR::powell_wetherall()} requires an \strong{interactive}
#'   selection of the regression interval when \code{reg_int} is not provided. In
#'   a non-interactive context (knitr, Rscript, \code{R CMD check}), it returns
#'   \code{NULL}. Therefore explicitly pass \code{reg_int}, for example
#'   \code{run_powell(lfq, reg_int = c(10, 25))}.
#'
#' @param lfq LFQ object (see \code{prepare_tropfish()}).
#' @param catch_columns catch columns to use (default : all).
#' @param ... arguments passed to \code{TropFishR::powell_wetherall()},
#'   notably \code{reg_int} (bounds of the regression interval).
#' @return S3 object \code{FishStockPowell}.
#' @seealso \code{\link{run_elefan_ga}}, \code{\link{run_growth_analysis}}
#' @examples
#' \dontrun{
#'   pw <- run_powell(lfq, reg_int = c(10, 25))
#'   summary(pw)
#' }
#' @export
run_powell <- function(lfq,
                       catch_columns=NULL,
                       ...){
  
  check_lfq(lfq)
  
  if(is.null(catch_columns))
    catch_columns <- seq_len(ncol(lfq$catch))
  
  res <- tryCatch(
    
    .time_it(
      
      TropFishR::powell_wetherall(
        
        lfq,
        
        catch_columns=catch_columns,
        
        ...
        
      )
      
    ),
    
    error=function(e){
      
      stop(
        "Powell-Wetherall error : ",
        conditionMessage(e),
        call.=FALSE
      )
      
    }
    
  )
  
  check_powell(res$value)
  
  structure(
    
    list(
      
      method="Powell-Wetherall",
      
      result=res$value,
      
      elapsed=res$elapsed
      
    ),
    
    class="FishStockPowell"
    
  )
  
}

#==============================================================
# Powell summary
#==============================================================

#' Summary of a FishStockPowell object
#' @param object FishStockPowell object
#' @param ... ignored
#' @export
summary.FishStockPowell <- function(object,...){
  
  cat("\n")
  cat("--------------------------------------\n")
  cat("Powell-Wetherall\n")
  cat("--------------------------------------\n")
  
  print(object$result)
  
  cat("\nTime :", object$elapsed,"\n")
  
  invisible(object)
  
}

#==============================================================
# K-Scan (fixed Linf from Powell-Wetherall)
#==============================================================

#' K-Scan : estimation of K with fixed \eqn{L_{inf}} (Powell-Wetherall)
#'
#' Second step of the TropFishR growth protocol. Once
#' \eqn{L_{inf}} is estimated by \code{\link{run_powell}}, it is fixed and a
#' range of \eqn{K} is scanned by ELEFAN : the value of \eqn{K} that maximizes the
#' restructuring score \eqn{R_n} provides a first estimate, serving
#' as an informed starting point for \code{\link{run_elefan_sa}} /
#' \code{\link{run_elefan_ga}}. This is a special case of the response
#' surface search, restricted to a single dimension (\eqn{K}).
#'
#' @param lfq Object of class \code{lfq} (TropFishR).
#' @param Linf Either a numeric value of \eqn{L_{inf}}, or a
#'   \code{FishStockPowell} object returned by \code{\link{run_powell}} (in which case
#'   \code{Linf_est} is extracted from it).
#' @param K_range Vector of \eqn{K} values to scan
#'   (default \code{seq(0.1, 2, 0.1)}).
#' @param MA Width of the moving average for restructuring (default 5).
#' @param ... Additional arguments passed to \code{TropFishR::ELEFAN()}.
#'
#' @return Object of class \code{FishStockKScan} : list containing
#'   \code{result} (raw ELEFAN output, including optimal \code{par$K} and
#'   \code{score_mat}), \code{Linf}, \code{K} (optimal K), \code{scan}
#'   (data.frame \code{K}/\code{score}) and \code{elapsed}.
#' @references Pauly, D. (1980). ELEFAN I. \emph{ICLARM Fishbyte}.
#'   Wetherall, J.A. (1986). A new method for estimating growth and mortality
#'   parameters from length-frequency data. \emph{ICLARM Fishbyte}.
#' @seealso \code{\link{run_powell}}, \code{\link{plot_response_surface}},
#'   \code{\link{run_elefan_sa}}, \code{\link{run_growth_analysis}}
#' @examples
#' \dontrun{
#'   pw <- run_powell(lfq)
#'   ks <- run_kscan(lfq, Linf = pw, K_range = seq(0.1, 1.5, 0.05))
#'   ks$K
#'   plot(ks)
#' }
#' @export
run_kscan <- function(lfq,
                      Linf,
                      K_range = seq(0.1, 2, 0.1),
                      MA = 5,
                      ...) {

  check_lfq(lfq)

  ## Accept a Powell object or a numeric value
  if (inherits(Linf, "FishStockPowell")) {
    Linf_val <- Linf$result$Linf_est
  } else if (is.numeric(Linf) && length(Linf) == 1L) {
    Linf_val <- Linf
  } else {
    stop("'Linf' must be a numeric value or a FishStockPowell object.",
         call. = FALSE)
  }

  if (!is.finite(Linf_val) || Linf_val <= 0)
    stop("Invalid Linf (", Linf_val, ").", call. = FALSE)

  if (length(K_range) < 2)
    stop("K_range must contain at least two values to scan.", call. = FALSE)

  res <- tryCatch(
    .time_it(
      TropFishR::ELEFAN(
        lfq,
        Linf_fix = Linf_val,
        K_range  = K_range,
        MA       = MA,
        plot     = FALSE,
        contour  = FALSE,
        hide.progressbar = TRUE,
        ...
      )
    ),
    error = function(e)
      stop("K-Scan error : ", conditionMessage(e), call. = FALSE)
  )

  out <- res$value

  ## Extract the score profile (K -> Rn)
  sm    <- as.matrix(out$score_mat)
  scan  <- data.frame(
    K     = as.numeric(rownames(sm)),
    score = as.numeric(sm[, 1])
  )

  structure(
    list(
      method  = "K-Scan (Linf fixed)",
      result  = out,
      Linf    = Linf_val,
      K       = out$par$K,
      scan    = scan,
      elapsed = res$elapsed
    ),
    class = "FishStockKScan"
  )
}

#' Summary of a K-Scan
#' @param object Object \code{FishStockKScan}.
#' @param ... Ignored.
#' @export
summary.FishStockKScan <- function(object, ...) {
  cat("\n--------------------------------------\n")
  cat("K-Scan (Linf fixed)\n")
  cat("--------------------------------------\n")
  cat("Fixed Linf :", round(object$Linf, 2), "\n")
  cat("Optimal K :", round(object$K, 3), "\n")
  best <- object$scan[which.max(object$scan$score), ]
  cat("Max Rn score :", round(best$score, 3), "at K =", round(best$K, 3), "\n")
  cat("Time :", object$elapsed, "\n")
  invisible(object)
}

#==============================================================
# Classic ELEFAN
#==============================================================

#' Classic ELEFAN (grid search)
#'
#' Wrapper around \code{TropFishR::ELEFAN()} : estimates \eqn{L_{inf}} and
#' \eqn{K} by grid search from length frequencies.
#'
#' @param lfq Object of class \code{lfq} (TropFishR).
#' @param Linf_range Vector of two values bounding the search for
#'   \eqn{L_{inf}} (optional).
#' @param K_range Vector of two values bounding the search for \eqn{K}
#'   (optional).
#' @param ... Additional arguments passed to \code{TropFishR::ELEFAN()}.
#' @return An object of class \code{FishStockELEFAN} : list containing
#'   \code{method}, \code{result} (raw TropFishR output) and \code{elapsed}.
#' @references Pauly, D. (1980). ELEFAN I. \emph{ICLARM Fishbyte}.
#' @seealso \code{\link{run_elefan_ga}}, \code{\link{run_elefan_sa}},
#'   \code{\link{run_growth_analysis}}
#' @examples
#' \dontrun{
#'   fit <- run_elefan(lfq, Linf_range = c(40, 55), K_range = c(0.4, 2))
#' }
#' @export
run_elefan <- function(lfq,
                       Linf_range=NULL,
                       K_range=NULL,
                       ...){
  
  
  check_lfq(lfq)
  
  
  check_growth_ranges <- function(Linf_range = NULL,
                                  K_range = NULL){
    
    
    warnings <- character()
    
    
    #---------------------------------
    # Linf check
    #---------------------------------
    
    if(!is.null(Linf_range)){
      
      if(length(Linf_range)!=2)
        stop("Linf_range must contain two values.")
      
      
      if(diff(Linf_range) <= 0)
        stop("Invalid Linf_range.")
      
      
      if(diff(Linf_range) > 200){
        
        warnings <- c(
          warnings,
          "The Linf range is very wide."
        )
        
      }
      
    }
    
    
    #---------------------------------
    # K check
    #---------------------------------
    
    if(!is.null(K_range)){
      
      if(length(K_range)!=2)
        stop("K_range must contain two values.")
      
      
      if(diff(K_range) <= 0)
        stop("Invalid K_range.")
      
      
      if(K_range[2] > 2){
        
        warnings <- c(
          warnings,
          paste0(
            "K_max = ",
            K_range[2],
            
            " is very high for an ELEFAN analysis."
          )
        )
        
      }
      
      
      if(K_range[2] > 5){
        
        warnings <- c(
          warnings,
          "The response surface risks being dominated by unrealistic K values."
        )
        
      }
      
    }
    
    
    if(length(warnings)>0){
      
      warning(
        paste(
          warnings,
          collapse="\n"
        )
      )
      
    }
    
    
    invisible(TRUE)
    
  }
  
  
  res <- tryCatch(
    
    .time_it(
      
      TropFishR::ELEFAN(
        
        lfq,
        
        Linf_range=Linf_range,
        
        K_range=K_range,
        
        ...
        
      )
      
    ),
    
    error=function(e){
      
      stop(
        "ELEFAN error : ",
        conditionMessage(e),
        call.=FALSE
      )
      
    }
    
  )
  
  
  structure(
    
    list(
      
      method="ELEFAN",
      
      result=res$value,
      
      elapsed=res$elapsed
      
    ),
    
    class="FishStockELEFAN"
    
  )
  
}

#==============================================================
# ELEFAN parameter verification
#==============================================================

#' Extract the growth parameters of an ELEFAN model
#'
#' Robustly retrieves \eqn{L_{inf}}, \eqn{K}, \code{t_anchor}, the performance
#' index \eqn{\phi'} as well as the fit scores
#' (\code{Rn_max}, \code{ASP}, \code{fESP}) and the computation time.
#'
#' \eqn{\phi'} is recomputed if missing : \eqn{\phi' = \log_{10}(K) +
#' 2\log_{10}(L_{inf})}.
#'
#' @param model Object \code{FishStockELEFAN} (see \code{\link{run_elefan_ga}}).
#' @return A named list : \code{Linf}, \code{K}, \code{t_anchor},
#'   \code{phiL}, \code{Rn_max}, \code{ASP}, \code{fESP}, \code{elapsed}.
#' @references Pauly, D., Munro, J.L. (1984). Once more on the comparison of
#'   growth in fish and invertebrates. \emph{ICLARM Fishbyte} 2(1).
#' @seealso \code{\link{growth_table}}, \code{\link{compare_growth_models}}
#' @examples
#' \dontrun{
#'   p <- extract_growth_parameters(fit)
#'   p$Linf; p$K; p$phiL
#' }
#' @export
extract_growth_parameters <- function(model){
  
  if(is.null(model$result$par))
    stop("Growth parameters missing.")
  
  
  p <- model$result$par
  
  result <- model$result
  
  
  out <- list(
    
    Linf = p[["Linf"]],
    
    K = p[["K"]],
    
    t_anchor =
      if("t_anchor" %in% names(p))
        p[["t_anchor"]]
    else
      NA_real_,
    
    
    phiL =
      if("phiL" %in% names(p))
        p[["phiL"]]
    else
      log10(p[["K"]]) +
      2*log10(p[["Linf"]])
    
  )
  
  
  ## TropFishR ELEFAN indicators
  
  out$Rn_max <-
    if("Rn_max" %in% names(result))
      result$Rn_max
  else
    NA_real_
  
  
  out$ASP <-
    if("ASP" %in% names(result))
      result$ASP
  else
    NA_real_
  
  
  out$fESP <-
    if("fESP" %in% names(result))
      result$fESP
  else
    NA_real_
  
  
  ## Time
  
  out$elapsed <-
    
    if(!is.null(model$elapsed))
      
      as.numeric(
        model$elapsed,
        units="secs"
      )
  
  else
    
    NA_real_
  
  
  out
  
}

#==============================================================
# ELEFAN summary
#==============================================================

#' @export
summary.FishStockELEFAN <- function(object,...){
  
  par <- extract_growth_parameters(object)
  
  cat("\n")
  
  cat("--------------------------------------\n")
  
  cat(object$method,"\n")
  
  cat("--------------------------------------\n\n")
  
  cat(sprintf("Linf      : %.3f\n", par$Linf))
  
  cat(sprintf("K         : %.4f\n", par$K))
  
  cat(sprintf("t_anchor  : %.3f\n", par$t_anchor))
  
  cat(sprintf("phiL      : %.3f\n", par$phiL))
  
  cat("\n")
  
  cat("Time :", round(as.numeric(object$elapsed, units="secs"),2),"sec\n")
  
  invisible(object)
  
}

#==============================================================
# Utility function
#==============================================================

#' Table of parameters of a growth model
#'
#' Formats the estimated parameters of an ELEFAN model as a single-row
#' \code{data.frame}, convenient for assembling reports.
#'
#' @param model Object \code{FishStockELEFAN}.
#' @return A single-row \code{data.frame} (Linf, K, t_anchor, phiL, Rn_max...).
#' @seealso \code{\link{extract_growth_parameters}}
#' @examples
#' \dontrun{
#'   growth_table(fit)
#' }
#' @export
growth_table <- function(model){
  
  p <- extract_growth_parameters(model)
  
  data.frame(
    
    method = model$method,
    
    Linf = p$Linf,
    
    K = p$K,
    
    t_anchor = p$t_anchor,
    
    phiL = p$phiL,
    
    Rn_max = p$Rn_max,
    
    ASP = p$ASP,
    
    fESP = p$fESP,
    
    Rn_max = p$Rn_max,
    
    elapsed = p$elapsed,
    
    stringsAsFactors = FALSE
    
  )
  
}

###############################################################
#
# Block 2B
#
# ELEFAN_GA
# ELEFAN_SA
# Model comparison
#
###############################################################

#==============================================================
# ELEFAN GA
#==============================================================

#' ELEFAN by genetic algorithm (ELEFAN_GA)
#'
#' Wrapper around \code{TropFishR::ELEFAN_GA()} : optimization of growth
#' parameters by genetic algorithm. The control arguments are filtered according
#' to the signature of the installed version of TropFishR, which ensures
#' compatibility across versions.
#'
#' @param lfq Object of class \code{lfq} (TropFishR).
#' @param control List of arguments passed to \code{ELEFAN_GA()} : typically
#'   \code{low_par}, \code{up_par}, \code{popSize}, \code{maxiter},
#'   \code{seasonalised}.
#' @param ... Additional arguments passed to \code{ELEFAN_GA()}.
#' @return An object of class \code{FishStockELEFAN}.
#' @references Taylor, M.H., Mildenberger, T.K. (2017). Extending
#'   electronic length frequency analysis in R. \emph{Fisheries Management and
#'   Ecology} 24(4).
#' @seealso \code{\link{run_elefan_sa}}, \code{\link{compare_growth_models}}
#' @examples
#' \dontrun{
#'   fit <- run_elefan_ga(lfq, control = list(
#'     low_par = list(Linf = 40, K = 0.4, t_anchor = 0),
#'     up_par  = list(Linf = 55, K = 2.0, t_anchor = 1),
#'     popSize = 40, maxiter = 30))
#' }
#' @export
run_elefan_ga <- function(lfq,
                          control=list(),
                          ...){
  
  check_lfq(lfq)
  
  #==============================================================
  # Biological check of growth search bounds
  #==============================================================
  
  check_growth_ranges <- function(Linf_range = NULL,
                                  K_range = NULL){
    
    
    warnings <- character()
    
    
    #---------------------------------
    # Linf check
    #---------------------------------
    
    if(!is.null(Linf_range)){
      
      if(length(Linf_range)!=2)
        stop("Linf_range must contain two values.")
      
      
      if(diff(Linf_range) <= 0)
        stop("Invalid Linf_range.")
      
      
      if(diff(Linf_range) > 200){
        
        warnings <- c(
          warnings,
          "The Linf range is very wide."
        )
        
      }
      
    }
    
    
    #---------------------------------
    # K check
    #---------------------------------
    
    if(!is.null(K_range)){
      
      if(length(K_range)!=2)
        stop("K_range must contain two values.")
      
      
      if(diff(K_range) <= 0)
        stop("Invalid K_range.")
      
      
      if(K_range[2] > 2){
        
        warnings <- c(
          warnings,
          paste0(
            "K_max = ",
            K_range[2],
            
            " is very high for an ELEFAN analysis."
          )
        )
        
      }
      
      
      if(K_range[2] > 5){
        
        warnings <- c(
          warnings,
          "The response surface risks being dominated by unrealistic K values."
        )
        
      }
      
    }
    
    
    if(length(warnings)>0){
      
      warning(
        paste(
          warnings,
          collapse="\n"
        )
      )
      
    }
    
    
    invisible(TRUE)
    
  }
  
  f <- get("ELEFAN_GA",
           envir=asNamespace("TropFishR"))
  
  args_fun <- names(formals(f))
  
  call_args <- c(
    list(lfq=lfq),
    control,
    list(...)
  )
  
  call_args <- call_args[
    names(call_args) %in% args_fun
  ]
  
  res <- tryCatch(
    
    .time_it(
      
      do.call(
        f,
        call_args
      )
      
    ),
    
    error=function(e){
      
      stop(
        "ELEFAN_GA : ",
        conditionMessage(e),
        call.=FALSE
      )
      
    }
    
  )
  
  structure(
    
    list(
      
      method="ELEFAN_GA",
      
      result=res$value,
      
      elapsed=res$elapsed
      
    ),
    
    class="FishStockELEFAN"
    
  )
  
}

#==============================================================
# ELEFAN SA
#==============================================================

#' ELEFAN by simulated annealing (ELEFAN_SA)
#'
#' Wrapper around \code{TropFishR::ELEFAN_SA()} : optimization of growth
#' parameters by simulated annealing.
#'
#' @param lfq Object of class \code{lfq} (TropFishR).
#' @param control List of arguments passed to \code{ELEFAN_SA()} : typically
#'   \code{low_par}, \code{up_par}, \code{SA_time}, \code{SA_temp},
#'   \code{seasonalised}.
#' @param ... Additional arguments passed to \code{ELEFAN_SA()}.
#' @return An object of class \code{FishStockELEFAN}.
#' @references Taylor, M.H., Mildenberger, T.K. (2017). \emph{Fisheries
#'   Management and Ecology} 24(4).
#' @seealso \code{\link{run_elefan_ga}}, \code{\link{compare_growth_models}}
#' @examples
#' \dontrun{
#'   fit <- run_elefan_sa(lfq, control = list(
#'     low_par = list(Linf = 40, K = 0.4, t_anchor = 0),
#'     up_par  = list(Linf = 55, K = 2.0, t_anchor = 1)))
#' }
#' @export
run_elefan_sa <- function(lfq,
                          control=list(),
                          ...){
  
  check_lfq(lfq)
  
  f <- get("ELEFAN_SA",
           envir=asNamespace("TropFishR"))
  
  args_fun <- names(formals(f))
  
  call_args <- c(
    list(lfq=lfq),
    control,
    list(...)
  )
  
  call_args <- call_args[
    names(call_args) %in% args_fun
  ]
  
  res <- tryCatch(
    
    .time_it(
      
      do.call(
        f,
        call_args
      )
      
    ),
    
    error=function(e){
      
      stop(
        "ELEFAN_SA : ",
        conditionMessage(e),
        call.=FALSE
      )
      
    }
    
  )
  
  structure(
    
    list(
      
      method="ELEFAN_SA",
      
      result=res$value,
      
      elapsed=res$elapsed
      
    ),
    
    class="FishStockELEFAN"
    
  )
  
}

#==============================================================
# ELEFAN score
#==============================================================

#' Fit score of an ELEFAN model
#'
#' Returns the fit score available in the TropFishR output, testing
#' successively \code{Rn_max}, \code{score} then \code{Rn} (compatibility
#' across versions).
#'
#' @param model Object \code{FishStockELEFAN}.
#' @return The numeric score, or \code{NA} if none is available.
#' @seealso \code{\link{best_growth_model}}
#' @examples
#' \dontrun{
#'   growth_score(fit)
#' }
#' @export
growth_score <- function(model){
  
  r <- model$result
  
  if("Rn_max" %in% names(r))
    return(r$Rn_max)
  
  if("score" %in% names(r))
    return(r$score)
  
  if("Rn" %in% names(r))
    return(r$Rn)
  
  NA
  
}

#==============================================================
# Comparative table
#==============================================================

#' Compare several growth models
#'
#' Assembles a comparative table of fitted ELEFAN models (classic, GA, SA),
#' with the estimated parameters, fit scores, computation time and an
#' indicator flagging estimates hitting the search bounds
#' (\code{boundary_flag}) - a classic sign of a poorly calibrated search space.
#'
#' @param ... Objects \code{FishStockELEFAN} to compare.
#' @param Linf_range Search bounds for \eqn{L_{inf}}, used to
#'   detect solutions on the bounds (optional).
#' @param K_range Search bounds for \eqn{K} (optional).
#' @return A \code{data.frame} : one row per method, with \code{method},
#'   \code{Linf}, \code{K}, \code{t_anchor}, \code{phiL}, \code{Rn_max},
#'   \code{ASP}, \code{fESP}, \code{elapsed}, \code{boundary_flag} and
#'   \code{rank_score}.
#' @seealso \code{\link{best_growth_model}}, \code{\link{run_growth_analysis}}
#' @examples
#' \dontrun{
#'   cmp <- compare_growth_models(ga, sa, Linf_range = c(40, 55))
#' }
#' @export
compare_growth_models <- function(...,
                                  Linf_range = NULL,
                                  K_range = NULL){
  
  models <- list(...)
  
  
  out <- lapply(models, function(m){
    
    p <- extract_growth_parameters(m)
    
    
    safe <- function(x){
      
      if(is.null(x) || length(x)==0){
        return(NA_real_)
      }
      
      as.numeric(x[1])
      
    }
    
    
    Linf <- safe(p$Linf)
    K <- safe(p$K)
    
    
    ##################################################
    # Bounds detection
    ##################################################
    
    boundary_flag <- FALSE
    
    
    if(!is.null(Linf_range)){
      
      if(Linf <= Linf_range[1] ||
         Linf >= Linf_range[2]){
        
        boundary_flag <- TRUE
        
      }
    }
    
    
    if(!is.null(K_range)){
      
      if(K <= K_range[1] ||
         K >= K_range[2]){
        
        boundary_flag <- TRUE
        
      }
    }
    
    
    data.frame(
      
      method = ifelse(
        is.null(m$method),
        NA,
        m$method
      ),
      
      Linf = Linf,
      
      K = K,
      
      t_anchor = safe(p$t_anchor),
      
      phiL = safe(p$phiL),
      
      Rn_max = safe(p$Rn_max),
      
      ASP = safe(p$ASP),
      
      fESP = safe(p$fESP),
      
      elapsed = safe(p$elapsed),
      
      boundary_flag = boundary_flag,
      
      stringsAsFactors = FALSE
      
    )
    
  })
  
  
  out <- do.call(rbind,out)
  
  
  rownames(out) <- NULL
  
  
  ##################################################
  # Ranking
  ##################################################
  
  out$rank_score <- rank(
    -out$Rn_max,
    ties.method="min"
  )
  
  
  out
  
}

#==============================================================
# Best model selection
#==============================================================

#' Select the model with the best score
#'
#' Simple selection : returns the row of the comparative table whose
#' fit score is the highest. The \code{Rn_max} column is used with
#' priority, falling back on \code{score} (compatibility across
#' TropFishR versions).
#'
#' @param table \code{data.frame} produced by \code{\link{compare_growth_models}}.
#' @return The \code{data.frame} reduced to the row of the best model.
#' @seealso \code{\link{best_growth_model}} for a multi-criteria selection.
#' @examples
#' \dontrun{
#'   choose_best_growth(cmp)
#' }
#' @export
choose_best_growth <- function(table){

  ## Robust scoring criterion : Rn_max (TropFishR) with priority,
  ## falling back on 'score' if present (compatibility with older versions).

  score_col <-
    if("Rn_max" %in% names(table))
      "Rn_max"
  else if("score" %in% names(table))
    "score"
  else
    NA_character_


  if(is.na(score_col)){

    warning("No score column (Rn_max/score) : returning the first model.")

    return(table[1, ])

  }


  ok <- !is.na(table[[score_col]])

  if(sum(ok)==0){

    return(table[1, ])

  }


  table[

    which.max(table[[score_col]]),

  ]

}

#==============================================================
# Print table
#==============================================================

#' Display the comparative table of growth models
#'
#' Formatted printout of the ELEFAN method comparison table.
#'
#' @param tab \code{data.frame} produced by \code{\link{compare_growth_models}}.
#' @return The table, invisibly (called for its display effect).
#' @seealso \code{\link{compare_growth_models}}
#' @examples
#' \dontrun{
#'   print_growth_comparison(cmp)
#' }
#' @export
print_growth_comparison <- function(tab){
  
  cat("\n")
  
  cat("---------------------------------\n")
  
  cat("Comparison of growth models\n")
  
  cat("---------------------------------\n\n")
  
  print(tab,row.names=FALSE)
  
  invisible(tab)
  
}

#' Multi-criteria selection of the best growth model
#'
#' Ranks the compared models according to the chosen criterion and returns the best one.
#' The \code{"balanced"} criterion combines the fit score (70 percent) and the
#' stability of parameters around the mean of the methods (30 percent), so as
#' not to retain a model with a high score but atypical parameters.
#'
#' @param comparison \code{data.frame} produced by
#'   \code{\link{compare_growth_models}}.
#' @param criterion Selection criterion : \code{"Rn_max"} (default),
#'   \code{"score"}, \code{"balanced"}, \code{"fastest"}, \code{"Linf"} or
#'   \code{"K"}.
#' @param verbose Logical. Display the selected model in the console.
#' @return An (invisible) list : \code{best} (row of the selected model),
#'   \code{ranking} (ordered table) and \code{criterion}.
#' @seealso \code{\link{compare_growth_models}}, \code{\link{choose_best_growth}}
#' @examples
#' \dontrun{
#'   sel <- best_growth_model(cmp, criterion = "balanced")
#'   sel$best
#' }
#' @export
best_growth_model <- function(comparison,
                              criterion = c(
                                "Rn_max",
                                "score",
                                "balanced",
                                "fastest",
                                "Linf",
                                "K"
                              ),
                              verbose = TRUE){
  
  criterion <- match.arg(criterion)
  
  
  df <- comparison
  
  
  if(criterion=="Rn_max"){
    
    df <- df[
      order(df$Rn_max,
            decreasing=TRUE),
    ]
    
  }
  
  
  if(criterion=="score"){
    
    df <- df[
      order(df$score,
            decreasing=TRUE),
    ]
    
  }
  
  
  if(criterion=="fastest"){
    
    df <- df[
      order(df$elapsed,
            decreasing=FALSE),
    ]
    
  }
  
  
  if(criterion=="Linf"){
    
    df <- df[
      order(df$Linf,
            decreasing=TRUE),
    ]
    
  }
  
  
  if(criterion=="K"){
    
    df <- df[
      order(df$K,
            decreasing=FALSE),
    ]
    
  }
  
  
  if(criterion=="balanced"){
    
    ## Normalization
    
    scale01 <- function(x){
      
      (x-min(x,na.rm=TRUE)) /
        (max(x,na.rm=TRUE)-min(x,na.rm=TRUE))
      
    }
    
    
    df$balanced_score <-
      
      0.70*scale01(df$Rn_max) +
      
      0.15*(1-scale01(abs(df$Linf -
                            mean(df$Linf)))) +
      
      0.15*(1-scale01(abs(df$K -
                            mean(df$K))))
    
    
    df <- df[
      order(df$balanced_score,
            decreasing=TRUE),
    ]
    
  }
  
  
  rownames(df) <- NULL
  
  
  best <- df[1,]
  
  
  if(verbose){
    
    cat("\n")
    cat("--------------------------------------\n")
    cat("Best growth model\n")
    cat("--------------------------------------\n\n")
    
    cat("Criterion :",criterion,"\n\n")
    
    cat("Selected method :",
        best$method,"\n\n")
    
    cat(sprintf(
      "Linf      : %.3f\n",
      best$Linf
    ))
    
    cat(sprintf(
      "K         : %.4f\n",
      best$K
    ))
    
    cat(sprintf(
      "phiL      : %.3f\n",
      best$phiL
    ))
    
    cat(sprintf(
      "Rn_max    : %.4f\n",
      best$Rn_max
    ))
    
    cat(sprintf(
      "Time      : %.2f sec\n",
      best$elapsed
    ))
    
  }
  
  
  invisible(
    list(
      best=best,
      ranking=df,
      criterion=criterion
    )
  )

}

#==============================================================
# von Bertalanffy curve (predictor)
#==============================================================

#' von Bertalanffy growth curve
#'
#' Generates the predicted length as a function of age via the official
#' function \code{TropFishR::VBGF()} (no algorithm rewritten).
#'
#' @param Linf,K Growth parameters.
#' @param t_anchor Time anchor point (TropFishR) ; ignored if \code{t0}
#'   is provided.
#' @param t0 Theoretical age at zero length (optional, alternative to
#'   \code{t_anchor}).
#' @param ages Vector of ages (years) for the prediction.
#' @return data.frame with columns \code{age} and \code{length}.
#' @examples
#' \dontrun{
#'   crb <- vbgf_curve(Linf = 86.4, K = 0.18, t_anchor = 0.83, ages = seq(0, 15, .1))
#' }
#' @export
vbgf_curve <- function(Linf,
                       K,
                       t_anchor = 0,
                       t0 = NULL,
                       ages = seq(0, 15, by = 0.1)) {

  check_tropfish()

  par <- list(Linf = Linf, K = K)

  if (!is.null(t0))
    par$t0 <- t0
  else
    par$t_anchor <- t_anchor


  L <- tryCatch(

    do.call(TropFishR::VBGF, list(param = par, t = ages)),

    error = function(e)
      stop("VBGF : ", conditionMessage(e), call. = FALSE)

  )


  data.frame(age = ages, length = as.numeric(L))

}

#==============================================================
# Orchestrator : complete growth analysis
#==============================================================

#' Complete growth analysis (multi-method + automatic selection)
#'
#' Runs the requested ELEFAN methods (classic, GA, SA), builds the
#' multi-criteria comparative table and automatically selects the best
#' model. Powell-Wetherall can be included as a diagnostic for initial
#' values (not compared, as it has a different structure).
#'
#' @param lfq Object \code{lfq}.
#' @param methods Subset of
#'   \code{c("elefan", "elefan_ga", "elefan_sa")}.
#' @param powell Logical. Include Powell-Wetherall as a diagnostic.
#' @param Linf_range,K_range Search ranges passed to the ELEFAN methods
#'   and to the bounds check.
#' @param criterion Selection criterion (see \code{best_growth_model}).
#' @param control List of control arguments passed to each method.
#' @param verbose Display progress and the best model.
#' @param ... Additional arguments passed to the ELEFAN functions.
#'
#' @return S3 object \code{FishStockGrowthAnalysis} : list \code{models},
#'   \code{comparison} (data.frame), \code{best}, \code{criterion} and
#'   \code{powell} (optional).
#' @examples
#' \dontrun{
#'   ga_an <- run_growth_analysis(lfq,
#'                                Linf_range = c(70, 100),
#'                                K_range = c(0.1, 0.5),
#'                                criterion = "Rn_max")
#'   ga_an$comparison
#'   ga_an$best
#' }
#' @export
run_growth_analysis <- function(lfq,
                                methods = c("elefan_ga", "elefan_sa", "elefan"),
                                powell = TRUE,
                                Linf_range = NULL,
                                K_range = NULL,
                                criterion = "Rn_max",
                                control = list(),
                                verbose = TRUE,
                                ...) {

  check_lfq(lfq)

  methods <- match.arg(
    methods,
    choices = c("elefan", "elefan_ga", "elefan_sa"),
    several.ok = TRUE
  )


  runners <- list(
    elefan    = run_elefan,
    elefan_ga = run_elefan_ga,
    elefan_sa = run_elefan_sa
  )


  models <- list()

  for (m in methods) {

    if (verbose) message("Fitting : ", m, " ...")

    ## Argument routing : run_elefan() takes Linf_range/K_range as
    ## direct arguments ; run_elefan_ga/sa() expect them in 'control'.
    fit <- tryCatch(

      if (m == "elefan") {
        do.call(runners[[m]],
                c(list(lfq, Linf_range = Linf_range, K_range = K_range),
                  control, list(...)))
      } else {
        ctrl <- control
        if (!is.null(Linf_range)) ctrl$Linf_range <- Linf_range
        if (!is.null(K_range))    ctrl$K_range    <- K_range
        runners[[m]](lfq, control = ctrl, ...)
      },

      error = function(e) {

        warning(m, " failed : ", conditionMessage(e))

        NULL

      }

    )

    if (!is.null(fit))
      models[[m]] <- fit

  }


  if (length(models) == 0)
    stop("No growth method succeeded.")


  ## Powell-Wetherall diagnostic (optional, not compared)

  pw    <- NULL
  kscan <- NULL

  if (isTRUE(powell)) {

    pw <- tryCatch(
      run_powell(lfq),
      error = function(e) {
        warning("Powell-Wetherall failed : ", conditionMessage(e))
        NULL
      }
    )

    ## K-Scan : optimal K at fixed Linf(Powell) (diagnostic, initial value)
    if (!is.null(pw)) {
      k_lo <- if (!is.null(K_range)) K_range[1] else 0.1
      k_hi <- if (!is.null(K_range)) K_range[2] else 2
      kscan <- tryCatch(
        run_kscan(lfq, Linf = pw,
                  K_range = seq(k_lo, k_hi,
                                length.out = 15)),
        error = function(e) {
          warning("K-Scan failed : ", conditionMessage(e))
          NULL
        }
      )
      if (verbose && !is.null(kscan))
        message("K-Scan (Linf=", round(kscan$Linf, 2),
                " fixed) -> K = ", round(kscan$K, 3))
    }

  }


  ## Comparison + selection

  comparison <- do.call(
    compare_growth_models,
    c(models, list(Linf_range = Linf_range, K_range = K_range))
  )

  sel <- best_growth_model(
    comparison,
    criterion = criterion,
    verbose = verbose
  )


  structure(

    list(
      models     = models,
      comparison = comparison,
      best       = sel$best,
      ranking    = sel$ranking,
      criterion  = criterion,
      powell     = pw,
      kscan      = kscan
    ),

    class = "FishStockGrowthAnalysis"

  )

}


#' @export
print.FishStockGrowthAnalysis <- function(x, ...) {

  cat("\n=====================================\n")
  cat(" stockflow : Growth analysis\n")
  cat("=====================================\n\n")

  cat("Fitted methods :",
      paste(names(x$models), collapse = ", "), "\n\n")

  print(x$comparison, row.names = FALSE)

  cat("\nBest model (", x$criterion, ") : ",
      x$best$method, "\n", sep = "")

  invisible(x)

}

#==============================================================
# ELEFAN response surface (confidence contour)
#==============================================================

#' ELEFAN response surface with confidence contour
#'
#' Plots the restructuring score matrix \eqn{R_n} as a function of
#' \eqn{L_{inf}} and \eqn{K} -- the \dQuote{banana-shaped} heat map of the
#' TropFishR protocol. The iso-contours delimit the confidence region
#' around the optimum ; the point of maximum score is marked. For
#' one-dimensional input (\code{FishStockKScan}), the function plots the
#' score profile \eqn{R_n(K)} at fixed \eqn{L_{inf}}.
#'
#' @param model Object \code{FishStockELEFAN} (grid search via
#'   \code{\link{run_elefan}}) or \code{FishStockKScan}
#'   (\code{\link{run_kscan}}).
#' @param conf Fraction of the maximum score delimiting the confidence
#'   region shown (default 0.95, i.e. 95\% of the optimal \eqn{R_n}).
#' @param title Chart title (optional).
#'
#' @return A \code{ggplot} object.
#' @seealso \code{\link{run_elefan}}, \code{\link{run_kscan}},
#'   \code{\link{plot_growth_modes}}
#' @examples
#' \dontrun{
#'   fit <- run_elefan(lfq, Linf_range = c(38, 46), K_range = c(0.2, 1))
#'   plot_response_surface(fit)
#' }
#' @export
plot_response_surface <- function(model, conf = 0.95, title = NULL) {

  ## --- one-dimensional case : K-Scan ---
  if (inherits(model, "FishStockKScan")) {
    d <- model$scan
    kbest <- model$K
    p <- ggplot2::ggplot(d, ggplot2::aes(.data$K, .data$score)) +
      ggplot2::geom_line(colour = "#2166ac", linewidth = 0.9) +
      ggplot2::geom_point(size = 1.4, colour = "#2166ac") +
      ggplot2::geom_vline(xintercept = kbest, linetype = 2,
                          colour = "#c0392b") +
      ggplot2::annotate("text", x = kbest,
                        y = max(d$score, na.rm = TRUE),
                        label = sprintf("K = %.3f", kbest),
                        hjust = -0.1, vjust = 1, size = 3.3,
                        colour = "#c0392b") +
      ggplot2::labs(
        title = title %||% "K-Scan : score profile (Linf fixed)",
        subtitle = sprintf("Linf = %.2f fixed", model$Linf),
        x = "K (an^-1)", y = "Restructuring score Rn") +
      .theme_fish()
    return(p)
  }

  if (!inherits(model, "FishStockELEFAN"))
    stop("'model' must be a FishStockELEFAN or FishStockKScan object.",
         call. = FALSE)

  sm <- model$result$score_mat
  if (is.null(sm))
    stop("No score matrix (score_mat) : rerun run_elefan() with ",
         "Linf_range and K_range ranges.", call. = FALSE)

  sm <- as.matrix(sm)
  if (ncol(sm) < 2 || nrow(sm) < 2)
    stop("The response surface requires a 2D grid (Linf x K). Use ",
         "plot_response_surface() on a K-Scan for a 1D profile.",
         call. = FALSE)

  ## Matrix -> long data.frame (rownames = K, colnames = Linf)
  Kv    <- as.numeric(rownames(sm))
  Linfv <- as.numeric(colnames(sm))
  grid  <- expand.grid(K = Kv, Linf = Linfv)
  grid$score <- as.numeric(sm)

  ## Optimum : the solution retained by ELEFAN (authoritative), not the
  ## which.max of the raw grid (an isolated peak may not be the refined
  ## solution).
  par <- model$result$par
  best <- if (!is.null(par$Linf) && !is.null(par$K)) {
    data.frame(Linf = par$Linf, K = par$K,
               score = model$result$Rn_max %||% max(grid$score, na.rm = TRUE))
  } else {
    grid[which.max(grid$score), ]
  }
  thr  <- conf * max(grid$score, na.rm = TRUE)

  ## interpolate=FALSE : each cell displays the actually computed score
  ## (interpolation smooths and may mask an isolated peak = the optimum).
  p <- ggplot2::ggplot(grid, ggplot2::aes(.data$Linf, .data$K)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$score),
                         interpolate = FALSE) +
    ggplot2::geom_contour(ggplot2::aes(z = .data$score),
                          colour = "white", linewidth = 0.3,
                          alpha = 0.6) +
    ggplot2::geom_contour(ggplot2::aes(z = .data$score),
                          breaks = thr, colour = "#c0392b",
                          linewidth = 0.9) +
    ggplot2::geom_point(data = best, size = 3, shape = 21,
                        fill = "#c0392b", colour = "white", stroke = 0.8) +
    ggplot2::annotate("text", x = best$Linf, y = best$K,
                      label = sprintf("Linf=%.1f\nK=%.2f", best$Linf, best$K),
                      hjust = -0.12, vjust = 0.5, size = 3.1,
                      colour = "#c0392b") +
    ggplot2::scale_fill_viridis_c(option = "D", name = "Rn") +
    ggplot2::labs(
      title = title %||% "ELEFAN response surface",
      subtitle = sprintf("Confidence contour : %.0f %% of the maximum Rn ; optimum marked",
                         conf * 100),
      x = "Linf", y = "K (an^-1)") +
    .theme_fish()

  p
}

#==============================================================
# Growth curve across modes (restructured LFQ)
#==============================================================

#' von Bertalanffy growth curve overlaid on modes (restructured LFQ)
#'
#' Plots the histogram of \strong{restructured} length frequencies
#' (positive / negative scores from \code{TropFishR::lfqRestructure()}) and
#' overlays one or more von Bertalanffy growth curves -- the
#' final visual check of the protocol : a good solution passes through
#' the modes (positive restructured bins) across the samples.
#'
#' @param lfq Object of class \code{lfq} (raw frequencies).
#' @param model Fitted growth object : \code{FishStockELEFAN},
#'   \code{FishStockKScan}, \code{FishStockGrowthAnalysis} (the best model
#'   is used), or a list \code{list(Linf=, K=, t_anchor=)}.
#' @param MA Width of the moving average for restructuring (default 5).
#'   Must match the one used during fitting.
#' @param n_cohorts Number of cohorts (shifted curves) to plot (default 6).
#' @param title Chart title (optional).
#'
#' @return A \code{ggplot} object.
#' @seealso \code{\link[TropFishR]{lfqRestructure}},
#'   \code{\link{plot_response_surface}}, \code{\link{run_growth_analysis}}
#' @examples
#' \dontrun{
#'   fit <- run_elefan_ga(lfq, Linf_range = c(38, 46), K_range = c(0.2, 1))
#'   plot_growth_modes(lfq, fit)
#' }
#' @export
plot_growth_modes <- function(lfq, model, MA = 5, n_cohorts = 6,
                              title = NULL) {

  check_lfq(lfq)

  ## --- extract Linf, K, t_anchor depending on model type ---
  gp <- .extract_growth_par(model)

  ## --- restructuring the frequencies ---
  lfqr <- tryCatch(
    TropFishR::lfqRestructure(lfq, MA = MA),
    error = function(e)
      stop("LFQ restructuring : ", conditionMessage(e), call. = FALSE)
  )

  ## --- restructured histogram in long format ---
  rc    <- lfqr$rcounts                # matrix midLengths x dates
  mids  <- lfqr$midLengths
  dates <- lfqr$dates
  df <- expand.grid(mid = mids, date = seq_along(dates))
  df$rc   <- as.numeric(rc)
  df$date <- dates[df$date]
  df$sign <- ifelse(df$rc >= 0, "positive", "negative")

  ## bar width = class step
  bw <- if (length(mids) > 1) stats::median(diff(mids)) else 1
  ## time scale : convert dates to decimal years
  yr <- as.numeric(format(dates, "%Y")) +
        (as.numeric(format(dates, "%j")) - 1) / 365
  date_num <- yr[match(df$date, dates)]
  ## visual amplitude of the bars (fraction of a year)
  span   <- if (length(unique(yr)) > 1) diff(range(yr)) else 1
  scale_bar <- 0.35 * span / max(abs(df$rc), na.rm = TRUE)

  df$xleft  <- date_num
  df$xright <- date_num + df$rc * scale_bar

  ## --- shifted VBGF curves (cohorts) ---
  t0_anchor <- gp$t_anchor %||% 0
  age_seq   <- seq(0, 15, by = 0.05)
  ymin <- min(yr); ymax <- max(yr) + 1
  ## one cohort roughly every year, shifted in time ; the model's t_anchor
  ## shifts the growth phase within the year.
  starts <- seq(floor(ymin) - n_cohorts, ceiling(ymax), by = 1)
  curves <- do.call(rbind, lapply(starts, function(s) {
    len <- gp$Linf * (1 - exp(-gp$K * (age_seq - t0_anchor)))
    len[len < 0] <- 0
    data.frame(cohort = s, t = s + age_seq, length = len)
  }))
  curves <- curves[curves$t >= ymin - 0.1 & curves$t <= ymax + 0.1 &
                   curves$length <= max(mids) * 1.02, ]

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = df,
      ggplot2::aes(xmin = pmin(.data$xleft, .data$xright),
                   xmax = pmax(.data$xleft, .data$xright),
                   ymin = .data$mid - bw / 2, ymax = .data$mid + bw / 2,
                   fill = .data$sign),
      colour = NA) +
    ggplot2::geom_line(
      data = curves,
      ggplot2::aes(x = .data$t, y = .data$length, group = .data$cohort),
      colour = "#c0392b", linewidth = 0.7, alpha = 0.9) +
    ggplot2::scale_fill_manual(
      values = c(positif = "#2166ac", negatif = "grey80"),
      name = "Score") +
    ggplot2::labs(
      title = title %||% "Growth curve across modes",
      subtitle = sprintf("VBGF (Linf=%.1f, K=%.2f) on restructured frequencies (MA=%d)",
                         gp$Linf, gp$K, MA),
      x = "year", y = "Length") +
    .theme_fish()

  p
}

## Internal helper : extract (Linf, K, t_anchor) from different objects
.extract_growth_par <- function(model) {
  if (is.list(model) && !is.null(model$Linf) && !is.null(model$K) &&
      is.null(attr(model, "class"))) {
    return(list(Linf = model$Linf, K = model$K,
                t_anchor = model$t_anchor %||% 0))
  }
  if (inherits(model, "FishStockGrowthAnalysis")) {
    ## $best (a row of the comparison) already carries Linf/K/t_anchor of
    ## the best model retained ; we read them directly (the keys of $models
    ## are lowercase 'elefan_ga', the method column uppercase).
    b <- model$best
    Linf <- b$Linf; K <- b$K
    ta   <- if (!is.null(b$t_anchor) && is.finite(b$t_anchor)) b$t_anchor else 0
    if (is.null(Linf) || is.null(K) || !is.finite(Linf) || !is.finite(K))
      stop("Growth parameters missing from the best model.",
           call. = FALSE)
    return(list(Linf = Linf, K = K, t_anchor = ta))
  }
  if (inherits(model, "FishStockKScan")) {
    ## object with flat fields (Linf, K) ; t_anchor under $result$par
    return(list(Linf = model$Linf, K = model$K,
                t_anchor = model$result$par$t_anchor %||% 0))
  }
  if (inherits(model, "FishStockELEFAN")) {
    res <- model$result
    return(list(Linf = res$par$Linf, K = res$par$K,
                t_anchor = res$par$t_anchor %||% 0))
  }
  stop("Model type not recognized for growth parameter extraction.",
       call. = FALSE)
}

#' Plots a K-Scan (Rn score profile as a function of K)
#' @param x Object \code{FishStockKScan}.
#' @param ... Arguments passed to \code{\link{plot_response_surface}}.
#' @return A \code{ggplot} object (invisible).
#' @export
plot.FishStockKScan <- function(x, ...) {
  print(plot_response_surface(x, ...))
  invisible(x)
}

#' Plots the response surface of a grid-search ELEFAN fit
#' @param x Object \code{FishStockELEFAN}.
#' @param ... Arguments passed to \code{\link{plot_response_surface}}.
#' @return A \code{ggplot} object (invisible).
#' @export
plot.FishStockELEFAN <- function(x, ...) {
  print(plot_response_surface(x, ...))
  invisible(x)
}

#==============================================================
# VBGF growth from age-length data (FSA)
#==============================================================

#' von Bertalanffy fit from age-length data (FSA)
#'
#' Complement to the ELEFAN methods (length-frequency) : fits the
#' von Bertalanffy model when individual ages are available, via the
#' official functions \code{FSA::vbStarts()} and \code{FSA::vbFuns()} and
#' \code{stats::nls()}. Bootstrap confidence intervals are computed
#' with \code{FSA::bootCase()} / \code{car::Boot()} if available.
#'
#' @param data data.frame containing age and length.
#' @param age_col,length_col Names of the age and length columns.
#' @param type von Bertalanffy parameterization (\code{"typical"} by default ;
#'   see \code{FSA::vbFuns()}).
#' @param boot Logical. Compute bootstrap CIs.
#' @param nboot Number of bootstrap resamples.
#'
#' @return S3 object \code{FishStockVBGF} (compatible with the parameter
#'   extractors of the mortality module via \code{$par}).
#'
#' @references
#'   von Bertalanffy, L. (1938). Human Biology 10.
#'   Ogle, D.H. (2016). Introductory Fisheries Analyses with R.
#'
#' @examples
#' \dontrun{
#'   vb <- fit_vbgf_age(agedata, "age", "TL", boot = TRUE)
#'   vb$par
#' }
#' @export
fit_vbgf_age <- function(data,
                         age_col,
                         length_col,
                         type = "typical",
                         boot = FALSE,
                         nboot = 200) {

  if (!requireNamespace("FSA", quietly = TRUE))
    stop("The 'FSA' package must be installed.", call. = FALSE)

  if (!all(c(age_col, length_col) %in% names(data)))
    stop("Age/length columns not found in 'data'.")


  d <- data.frame(
    age = suppressWarnings(as.numeric(data[[age_col]])),
    len = suppressWarnings(as.numeric(data[[length_col]]))
  )

  d <- d[is.finite(d$age) & is.finite(d$len), ]

  if (nrow(d) < 10)
    stop("Insufficient sample size (<10) to fit the VBGF.")


  ## initial values and official FSA function
  sv <- tryCatch(
    FSA::vbStarts(len ~ age, data = d, type = type),
    error = function(e)
      stop("vbStarts : ", conditionMessage(e), call. = FALSE)
  )

  vb <- FSA::vbFuns(type)


  fit <- tryCatch(

    stats::nls(len ~ vb(age, Linf, K, t0),
               data = d, start = sv,
               control = stats::nls.control(maxiter = 200)),

    error = function(e)
      stop("nls (VBGF) : ", conditionMessage(e), call. = FALSE)

  )


  par <- as.list(stats::coef(fit))

  par$phiL <- log10(par$K) + 2 * log10(par$Linf)


  ## bootstrap confidence intervals (optional)
  ci <- NULL

  if (isTRUE(boot)) {

    ci <- tryCatch({

      bc <-
        if (exists("bootCase", where = asNamespace("FSA")))
          FSA::bootCase(fit, B = nboot)
        else if (requireNamespace("car", quietly = TRUE))
          car::Boot(fit, R = nboot)
        else
          NULL

      if (!is.null(bc))
        apply(as.matrix(bc), 2,
              stats::quantile, probs = c(0.025, 0.975),
              na.rm = TRUE)
      else
        NULL

    }, error = function(e) NULL)

  }


  structure(

    list(
      method = paste0("VBGF_nls_", type),
      par    = par,
      ci     = ci,
      model  = fit,
      n      = nrow(d)
    ),

    class = "FishStockVBGF"

  )

}


#' @export
print.FishStockVBGF <- function(x, ...) {

  cat("\n")
  cat("--------------------------------------\n")
  cat("von Bertalanffy (age-length, FSA/nls)\n")
  cat("--------------------------------------\n\n")

  cat(sprintf("Linf      : %.3f\n", x$par$Linf))
  cat(sprintf("K         : %.4f\n", x$par$K))
  cat(sprintf("t0        : %.3f\n", x$par$t0))
  cat(sprintf("phiL      : %.3f\n", x$par$phiL))
  cat(sprintf("n         : %d\n", x$n))

  if (!is.null(x$ci)) {
    cat("\n95%% bootstrap CI :\n")
    print(round(x$ci, 4))
  }

  invisible(x)

}