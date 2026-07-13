###############################################################
#
# stockflow
#
# Module : lbb.R
#
# LBB - Length-Based Bayesian Biomass estimation
# (Froese, Winker, Coro et al. 2018)
#
# Principle: robust wrapper around the OFFICIAL 'LBB' package.
# The Bayesian algorithm (JAGS) is never rewritten; only the
# data preparation, calling and indicator extraction steps
# are handled here. The length reference points
# (Lopt, Lc_opt) use the published formulas of Froese et al.
#
###############################################################


#==============================================================
# Preparing data in LBB format
#==============================================================

#' Preparation of length data in LBB format
#'
#' Builds a data.frame in the format expected by the \code{LBB} package
#' (columns \code{Stock}, \code{Year}, \code{Length}, \code{CatchNo}) from
#' individual lengths or a frequency table.
#'
#' @param x data.frame of individual lengths OR of frequencies.
#' @param length_col Name of the length column.
#' @param year_col Name of the year column.
#' @param freq_col Name of the count/frequency column. If \code{NULL},
#'   each row counts as one individual (individual data).
#' @param stock Stock name (identifier).
#' @param bin_width Width of the length classes (binning). If
#'   \code{NULL}, no binning.
#'
#' @return data.frame \code{Stock, Year, Length, CatchNo}.
#' @references Froese, R. et al. (2018). ICES J. Mar. Sci. 75(6).
#' @export
prepare_lbb_data <- function(x,
                             length_col,
                             year_col,
                             freq_col = NULL,
                             stock = "stock",
                             bin_width = NULL) {

  if (!all(c(length_col, year_col) %in% names(x)))
    stop("Length/year columns not found.")


  d <- data.frame(
    Year   = x[[year_col]],
    Length = suppressWarnings(as.numeric(x[[length_col]])),
    CatchNo = if (is.null(freq_col)) 1
              else suppressWarnings(as.numeric(x[[freq_col]]))
  )

  d <- d[is.finite(d$Length) & d$Length > 0 & !is.na(d$Year), ]


  ## binning into length classes if requested
  if (!is.null(bin_width)) {

    d$Length <- floor(d$Length / bin_width) * bin_width + bin_width / 2

  }


  agg <- stats::aggregate(CatchNo ~ Year + Length, data = d, FUN = sum)

  out <- data.frame(
    Stock   = stock,
    Year    = agg$Year,
    Length  = agg$Length,
    CatchNo = agg$CatchNo
  )

  out <- out[order(out$Year, out$Length), ]
  rownames(out) <- NULL

  out

}



#==============================================================
# Backend "Froese script" (LBB_ggplot.R, R2jags)
#
# R. Froese's official script (LBB_ggplot.R, ggplot2 adaptation by
# K. Ba) is procedural: it reads an ID file (priors) + a
# lengths file (Year, Length, CatchNo), runs JAGS and writes
# LBB_Results/<Stock>_Tableau_Parametres.csv. These functions prepare
# the files, run the script without modifying it scientifically,
# and retrieve the results table.
#==============================================================

#' Writing the ID file (priors) for Froese's LBB script
#'
#' @param stock Stock identifier.
#' @param lf_file Name of the lengths file (columns Year, Length, CatchNo).
#' @param species Species name.
#' @param mm Logical. Lengths in mm (\code{TRUE}) or cm (\code{FALSE}).
#' @param MK_prior,Linf_prior,Lm50,Lc_user,Lstart Optional priors
#'   (\code{MK.user}, \code{Linf.user}, \code{Lm50}, \code{Lc.user},
#'   \code{Lstart.user}).
#' @param gaussian_sel Logical. Gaussian selectivity (gillnets)?
#' @param merge_lf Logical. Merge years with low sample size?
#' @param pile Correction for piling effect (0, 1 or 999).
#' @param start_year,end_year Optional time bounds.
#' @param file Output path of the ID file (CSV).
#' @return The ID data.frame (invisible), written to \code{file}.
#' @export
write_lbb_id <- function(stock,
                         lf_file,
                         species = stock,
                         mm = TRUE,
                         MK_prior = 1.5,
                         Linf_prior = NA,
                         Lm50 = NA,
                         Lc_user = NA,
                         Lstart = NA,
                         gaussian_sel = FALSE,
                         merge_lf = FALSE,
                         pile = 0,
                         start_year = NA,
                         end_year = NA,
                         file = "LBB_ID.csv") {

  id <- data.frame(
    Stock        = stock,
    Species      = species,
    File         = lf_file,
    mm.user      = mm,
    GausSel      = gaussian_sel,
    MergeLF      = merge_lf,
    Pile         = pile,
    Gears.user   = NA,
    StartYear    = start_year,
    EndYear      = end_year,
    Years.user   = NA,
    Linf.user    = Linf_prior,
    Lcut.user    = NA,
    Year.select  = NA,
    Lc.user      = Lc_user,
    MK.user      = MK_prior,
    Lstart.user  = Lstart,
    Lm50         = Lm50,
    Comment      = "generated by stockflow::write_lbb_id()",
    stringsAsFactors = FALSE
  )

  utils::write.csv(id, file, row.names = FALSE)
  invisible(id)

}


#' Writing the input files for the LBB script (data + ID)
#'
#' Writes to \code{workdir} the lengths file (Year, Length, CatchNo) and
#' the ID file (priors), ready to feed Froese's \code{LBB_ggplot.R} script.
#'
#' @param data data.frame in the \code{prepare_lbb_data()} format (Stock, Year,
#'   Length, CatchNo).
#' @param stock Stock identifier (default: first value of
#'   \code{data$Stock}).
#' @param workdir Working directory.
#' @param ... Priors passed to \code{write_lbb_id()} (mm, MK_prior, ...).
#'
#' @section Units (classic pitfall):
#'   Froese's script imposes an asymmetric convention:
#'   \itemize{
#'     \item the lengths in the data file are \strong{always in mm};
#'     \item the priors in the ID file (\code{Linf.user}, \code{Lc.user},
#'           \code{Lm50}) are in the working unit, indicated by
#'           \code{mm.user}.
#'   }
#'   With \code{mm = FALSE} (working in cm), the script itself divides the
#'   lengths by 10. So convert your lengths to mm before calling, and
#'   express the priors in cm.
#'
#' @return List: \code{lf_file}, \code{id_file}, \code{workdir}, \code{stock}.
#' @export
write_lbb_inputs <- function(data,
                             stock = NULL,
                             workdir = tempdir(),
                             ...) {

  req <- c("Year", "Length", "CatchNo")
  if (!all(req %in% names(data)))
    stop("'data' must contain: ", paste(req, collapse = ", "),
         " (see prepare_lbb_data()).")

  if (is.null(stock))
    stock <- if ("Stock" %in% names(data)) as.character(data$Stock[1]) else "stock"

  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

  safe_stock <- gsub("[^A-Za-z0-9_.-]", "_", stock)
  lf_name    <- paste0(safe_stock, "_LF.csv")

  ## Froese's script matches the lengths file and the ID file by
  ## the 'Stock' column. Without it (or with a different value), it stops
  ## with: "Stock ID in ID file does not correspond to the Stock ID in DAT ...".
  lf <- data.frame(
    Stock   = stock,
    Year    = data$Year,
    Length  = data$Length,
    CatchNo = data$CatchNo,
    stringsAsFactors = FALSE
  )

  utils::write.csv(lf, file.path(workdir, lf_name), row.names = FALSE)

  id_file <- file.path(workdir, paste0(safe_stock, "_ID.csv"))
  write_lbb_id(stock = stock, lf_file = lf_name, file = id_file, ...)

  list(lf_file = lf_name,
       id_file = basename(id_file),
       workdir = workdir,
       stock = stock)

}


## Non-destructive patch of the script (Stock / ID.File / rm() / output_dir)
.patch_lbb_script <- function(script_path, stock, id_file, out_path) {

  lines <- readLines(script_path, warn = FALSE)

  lines <- vapply(lines, function(l) {

    if (grepl("^\\s*rm\\(list\\s*=\\s*ls", l))
      return(paste0("# (disabled by stockflow) ", l))

    if (grepl("^\\s*Stock\\s*<-", l))
      return(sprintf('Stock <- "%s"', stock))

    if (grepl("^\\s*ID\\.File\\s*<-", l))
      return(sprintf('ID.File <- "%s"', id_file))

    if (grepl('^\\s*output_dir\\s*<-', l))
      return(sprintf('output_dir <- "%s"', out_path))

    l

  }, character(1), USE.NAMES = FALSE)

  patched <- tempfile(fileext = ".R")
  writeLines(lines, patched)
  patched

}


#' Running Froese's LBB script (R2jags backend)
#'
#' Prepares the input files, runs the official
#' \code{LBB_ggplot.R} script (provided in \code{inst/scripts/}) without
#' modifying its scientific core, then retrieves the annual results table.
#'
#' @param data data.frame \code{prepare_lbb_data()} (Stock, Year, Length,
#'   CatchNo).
#' @param stock Stock identifier.
#' @param workdir Working directory (outputs, including ggplot2 figures,
#'   are written there under \code{LBB_Results/}).
#' @param script Path to the LBB script. Defaults to the one provided by
#'   the package.
#' @param ... Priors passed to \code{write_lbb_inputs()} / \code{write_lbb_id()}
#'   (mm, MK_prior, Linf_prior, Lm50, ...).
#'
#' @return S3 object \code{FishStockLBB} (\code{fit} = raw Ldat data.frame,
#'   \code{summary} = harmonized indicators, \code{ref_levels},
#'   \code{output_dir}).
#'
#' @section Prerequisites:
#'   Packages \code{R2jags}, \code{Hmisc}, \code{ggplot2}, \code{crayon} and
#'   the \strong{JAGS} software must be installed.
#'
#' @references Froese, R. et al. (2018). ICES J. Mar. Sci. 75(6): 2004-2015.
#' @examples
#' \dontrun{
#'   dd  <- prepare_lbb_data(freq, "LCT", "annee", "n", stock = "Penaeus")
#'   fit <- run_lbb_froese(dd, mm = TRUE, MK_prior = 1.5, Linf_prior = 45)
#'   fit$summary
#' }
#' @export
run_lbb_froese <- function(data,
                           stock = NULL,
                           workdir = file.path(tempdir(), "LBB_run"),
                           script = NULL,
                           ...) {

  if (is.null(script))
    script <- system.file("scripts", "LBB_ggplot.R", package = "stockflow")

  if (script == "" || !file.exists(script))
    stop("Script 'LBB_ggplot.R' not found. Please provide 'script'.", call. = FALSE)

  if (!requireNamespace("R2jags", quietly = TRUE))
    stop("The 'R2jags' package (and JAGS) must be installed.", call. = FALSE)


  inp <- write_lbb_inputs(data, stock = stock, workdir = workdir, ...)

  out_dir <- file.path(workdir, "LBB_Results")

  patched <- .patch_lbb_script(script, inp$stock, inp$id_file, "LBB_Results")


  ## execution in workdir (script's relative paths) via a dedicated env
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(workdir)

  env <- new.env()

  tryCatch(
    sys.source(patched, envir = env),
    error = function(e)
      stop("Running the LBB script: ", conditionMessage(e), call. = FALSE)
  )


  safe_stock <- gsub("[^A-Za-z0-9_.-]", "_", inp$stock)
  res_csv <- file.path("LBB_Results",
                       paste0(inp$stock, "_Tableau_Parametres.csv"))

  if (!file.exists(res_csv))
    res_csv <- file.path("LBB_Results",
                         paste0(safe_stock, "_Tableau_Parametres.csv"))

  if (!file.exists(res_csv))
    stop("LBB results table not found in ", out_dir, ".")

  Ldat <- utils::read.csv(res_csv, stringsAsFactors = FALSE)

  smy <- extract_lbb(Ldat)


  structure(
    list(fit = Ldat,
         summary = smy,
         ref_levels = .lbb_ref_levels(smy),
         output_dir = normalizePath(out_dir, mustWork = FALSE)),
    class = "FishStockLBB"
  )

}


#==============================================================
# Extracting LBB indicators
#==============================================================

#' Extraction of annual indicators from an LBB output
#'
#' Robustly searches for the key indicators in the object returned by
#' the \code{LBB} package and assembles them into an annual data.frame.
#'
#' @param fit Object returned by the LBB function.
#' @return data.frame: Year, Linf, Lc, Lopt, MK, FK, FM, ZK, BB0, BBmsy,
#'   LmeanLopt (columns available depending on the output).
#' @export
extract_lbb <- function(fit) {

  ## LBB often stores one data.frame per year; we try several
  ## locations and harmonize the names via patterns.

  df <- .find_lbb_table(fit)

  if (is.null(df))
    stop("Unable to extract an indicator table from the LBB output.")


  pick <- function(patterns) {

    for (p in patterns) {

      hit <- grep(p, names(df), ignore.case = TRUE, value = TRUE)

      if (length(hit) > 0)
        return(df[[hit[1]]])

    }

    rep(NA_real_, nrow(df))

  }


  out <- data.frame(
    Year      = pick(c("^year$", "yr")),
    Linf      = pick(c("^Linf", "L_inf")),
    Lc        = pick(c("^Lc$", "Lc_", "L50")),
    Lopt      = pick(c("^Lopt", "r\\.?Lopt", "L_opt")),
    MK        = pick(c("^MK$", "M_K", "M\\.K", "MdivK")),
    FK        = pick(c("^FK$", "F_K", "F\\.K", "FdivK")),
    FM        = pick(c("^FM$", "F_M", "F\\.M", "FdivM")),
    ZK        = pick(c("^ZK$", "Z_K", "Z\\.K")),
    BB0       = pick(c("BB0", "B_B0", "B\\.B0", "BdivB0")),
    BBmsy     = pick(c("BBmsy", "B_Bmsy", "B\\.Bmsy", "BdivBmsy")),
    LmeanLopt = pick(c("Lmean.*Lopt", "LmeanLopt")),
    stringsAsFactors = FALSE
  )

  ## F/M derived if absent but FK and MK are present
  if (all(is.na(out$FM)) &&
      !all(is.na(out$FK)) && !all(is.na(out$MK)))
    out$FM <- out$FK / out$MK

  out[order(out$Year), ]

}


.find_lbb_table <- function(fit) {

  ## 1. the object is already a data.frame
  if (is.data.frame(fit)) return(fit)

  ## 2. look for a data.frame element containing years + indicators
  if (is.list(fit)) {

    for (nm in names(fit)) {

      el <- fit[[nm]]

      if (is.data.frame(el) &&
          any(grepl("year", names(el), ignore.case = TRUE)) &&
          any(grepl("BB0|Bmsy|FM|MK", names(el), ignore.case = TRUE)))
        return(el)

    }

    ## fallback: first data.frame found
    for (nm in names(fit))
      if (is.data.frame(fit[[nm]]))
        return(fit[[nm]])

  }

  NULL

}


.lbb_ref_levels <- function(smy) {

  if (is.null(smy) || nrow(smy) == 0)
    return(NULL)

  last <- smy[nrow(smy), ]

  data.frame(
    Year_last = last$Year,
    BB0       = last$BB0,
    BBmsy     = last$BBmsy,
    FM        = last$FM,
    statut    = .lbb_status(last$BBmsy, last$FM),
    stringsAsFactors = FALSE
  )

}


.lbb_status <- function(BBmsy, FM) {

  if (is.na(BBmsy) || is.na(FM))
    return(NA_character_)

  if (BBmsy >= 1 && FM <= 1) "Healthy (B>=Bmsy, F<=M)"
  else if (BBmsy < 1 && FM > 1) "Overexploited (B<Bmsy, F>M)"
  else "Intermediate"

}


#==============================================================
# Length reference points (published Froese formulas)
#==============================================================

#' Length reference points (Lopt, Lc_opt)
#'
#' Computes the optimal length \eqn{L_{opt}} and the optimal length at
#' first capture \eqn{L_{c,opt}} from the ratios \eqn{M/K} and
#' \eqn{F/M}, according to the published formulas of Froese et al. This is not
#' the LBB algorithm but analytical reference relationships.
#'
#' @param Linf Asymptotic length.
#' @param MK Ratio \eqn{M/K} (default 1.5).
#' @param FM Ratio \eqn{F/M} (default 1, proxy \eqn{F = M}).
#' @return List: \code{Lopt}, \code{Lopt_Linf}, \code{Lc_opt},
#'   \code{Lc_opt_Linf}.
#' @references
#'   Froese, R. et al. (2016). Minimizing the impact of fishing.
#'   Fish and Fisheries 17(3).
#'   Froese, R. et al. (2018). ICES J. Mar. Sci. 75(6).
#' @examples
#'   lbb_reference_points(Linf = 45, MK = 1.5, FM = 1)
#' @export
lbb_reference_points <- function(Linf, MK = 1.5, FM = 1) {

  ## Lopt / Linf = 3 / (3 + M/K)
  Lopt_Linf <- 3 / (3 + MK)

  ## Lc_opt / Linf = (2 + 3 * F/M) / ((1 + F/M) * (3 + M/K))
  Lc_opt_Linf <- (2 + 3 * FM) / ((1 + FM) * (3 + MK))

  list(
    Lopt        = Lopt_Linf * Linf,
    Lopt_Linf   = Lopt_Linf,
    Lc_opt      = Lc_opt_Linf * Linf,
    Lc_opt_Linf = Lc_opt_Linf
  )

}


#==============================================================
# S3 methods
#==============================================================

#' @export
print.FishStockLBB <- function(x, ...) {

  cat("\n=====================================\n")
  cat(" stockflow : LBB (Length-Based Bayesian Biomass)\n")
  cat("=====================================\n\n")

  if (!is.null(x$summary)) {

    cat("Annual indicators (excerpt):\n\n")
    print(utils::tail(x$summary, 5), row.names = FALSE)

  }

  if (!is.null(x$ref_levels)) {

    cat("\nReference levels (last year):\n")
    print(x$ref_levels, row.names = FALSE)

  }

  invisible(x)

}


#' @export
summary.FishStockLBB <- function(object, ...) {

  print(object)

  invisible(object$summary)

}


#' Plot of LBB indicators (B/B0, B/Bmsy, F/M)
#' @param x Object \code{FishStockLBB}.
#' @param ... Ignored.
#' @export
plot.FishStockLBB <- function(x, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {

    ## basic graphical fallback
    with(x$summary, {
      graphics::plot(Year, BBmsy, type = "b",
                     ylab = "B/Bmsy", xlab = "year")
      graphics::abline(h = 1, lty = 2)
    })

    return(invisible(x))

  }


  df <- x$summary

  long <- data.frame(
    Year = rep(df$Year, 3),
    ind  = rep(c("B/B0", "B/Bmsy", "F/M"), each = nrow(df)),
    val  = c(df$BB0, df$BBmsy, df$FM)
  )

  p <- ggplot2::ggplot(long, ggplot2::aes(Year, val)) +
    ggplot2::geom_line(colour = "#12507b") +
    ggplot2::geom_point(size = 1) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
    ggplot2::facet_wrap(~ ind, scales = "free_y") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "LBB: status indicators") +
    ggplot2::theme_bw()

  print(p)

  invisible(p)

}