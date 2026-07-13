###############################################################
#
# stockflow
#
# Module : report.R
#
# Automatic generation of reports (HTML / PDF / DOCX).
#
# Principle: the S3 objects produced by the modules (growth,
# mortality, allometry, LBB, LBSPR, MSE) are aggregated into an
# object 'FishStockResults', then rendered via the OFFICIAL package
# 'rmarkdown' from a parameterized template supplied in
# inst/rmarkdown/. No result is recomputed here.
#
###############################################################

#==============================================================
# Dependencies
#==============================================================

.check_rmarkdown <- function() {

  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("The 'rmarkdown' package must be installed.", call. = FALSE)

  if (!requireNamespace("knitr", quietly = TRUE))
    stop("The 'knitr' package must be installed.", call. = FALSE)

  invisible(TRUE)

}


#==============================================================
# Aggregation of results
#==============================================================

#' Aggregation of pipeline results
#'
#' Gathers the objects produced by the different modules into a single
#' object, consumable by \code{fishstock_report()}. All arguments are
#' optional: only the available sections will be rendered.
#'
#' @param growth Object \code{FishStockGrowthAnalysis} (or \code{FishStockELEFAN}).
#' @param mortality Object \code{FishStockMortality}.
#' @param allometry Object \code{FishStockLW}.
#' @param lbb Object \code{FishStockLBB}.
#' @param lbspr Object \code{FishStockLBSPR}.
#' @param mse Object \code{FishStockMSE}.
#' @param maturity data.frame or list with \code{L50}, \code{L95} (optional).
#' @param ypr Object \code{FishStockYPR} (per-recruit analysis; \code{run_ypr()}).
#' @param tac Object \code{FishStockTAC} (\code{compute_tac()}).
#' @param quota Object \code{FishStockQuota} (\code{allocate_quota()}).
#' @param closure Object \code{FishStockClosure} (\code{optimal_closure()}).
#' @param size_advice Object \code{FishStockSizeAdvice} (\code{recommend_sizes()}).
#' @param meta List of metadata: \code{species}, \code{stock}, \code{area},
#'   \code{period}, \code{n_ind}, \code{length_unit}, \code{author}, ...
#'
#' @return S3 object \code{FishStockResults}.
#' @examples
#' \dontrun{
#'   res <- collect_results(growth = growth, mortality = mort,
#'                          allometry = lw, lbspr = fit_sp, mse = mse,
#'                          meta = list(species = "Penaeus notialis",
#'                                      stock = "Example EEZ"))
#' }
#' @export
collect_results <- function(growth = NULL,
                            mortality = NULL,
                            allometry = NULL,
                            lbb = NULL,
                            lbspr = NULL,
                            mse = NULL,
                            maturity = NULL,
                            ypr = NULL,
                            tac = NULL,
                            quota = NULL,
                            closure = NULL,
                            size_advice = NULL,
                            meta = list()) {

  structure(

    list(
      growth      = growth,
      mortality   = mortality,
      allometry   = allometry,
      lbb         = lbb,
      lbspr       = lbspr,
      mse         = mse,
      maturity    = maturity,
      ypr         = ypr,
      tac         = tac,
      quota       = quota,
      closure     = closure,
      size_advice = size_advice,
      meta        = meta,
      created     = Sys.time()
    ),

    class = "FishStockResults"

  )

}


#' @export
print.FishStockResults <- function(x, ...) {

  cat("\n=====================================\n")
  cat(" stockflow : aggregated results\n")
  cat("=====================================\n\n")

  if (length(x$meta)) {
    for (nm in names(x$meta))
      cat(sprintf("%-12s : %s\n", nm, paste(x$meta[[nm]], collapse = ", ")))
    cat("\n")
  }

  mods <- c("growth", "mortality", "allometry", "maturity",
            "lbb", "lbspr", "mse")

  for (m in mods)
    cat(sprintf("%-12s : %s\n", m,
                if (is.null(x[[m]])) "-" else "available"))

  invisible(x)

}


#==============================================================
# Extraction of summary tables (used by the template)
#==============================================================

#' Summary table of stock status
#'
#' Assembles into a data.frame the available key indicators (growth,
#' mortality, LBB, LBSPR) for the report header.
#'
#' @param results Object \code{FishStockResults}.
#' @return data.frame: Indicator, Value, Source.
#' @export
summary_table <- function(results) {

  if (!inherits(results, "FishStockResults"))
    stop("'results' must be a FishStockResults object.")

  rows <- list()

  add <- function(ind, val, src) {
    if (is.null(val) || length(val) == 0 || all(is.na(val))) return(invisible(NULL))
    rows[[length(rows) + 1]] <<-
      data.frame(Indicator = ind,
                 Value = if (is.numeric(val)) format(round(val, 3)) else as.character(val),
                 Source = src, stringsAsFactors = FALSE)
  }

  ## Growth
  g <- results$growth
  if (!is.null(g)) {
    b <- if (!is.null(g$best)) g$best else NULL
    if (!is.null(b)) {
      add("Linf", b$Linf, "Growth (ELEFAN)")
      add("K",    b$K,    "Growth (ELEFAN)")
      add("phiL", b$phiL, "Growth (ELEFAN)")
    }
  }

  ## Allometry
  a <- results$allometry
  if (!is.null(a)) {
    add("a (length-weight)", a$a, "Allometry")
    add("b (length-weight)", a$b, "Allometry")
    add("Allometry type", a$allometry, "Allometry")
  }

  ## maturity
  m <- results$maturity
  if (!is.null(m)) {
    add("L50", m$L50, "maturity")
    add("L95", m$L95, "maturity")
  }

  ## Mortality
  mo <- results$mortality
  if (!is.null(mo)) {
    add("Z", mo$Z$Z, "Mortality")
    add("M (consensus)", mo$M_consensus, "Mortality")
    add("F", mo$summary$F, "Mortality")
    add("E = F/Z", mo$summary$E, "Mortality")
    add("Status (E)", mo$summary$statut, "Mortality")
  }

  ## LBB
  l <- results$lbb
  if (!is.null(l) && !is.null(l$ref_levels)) {
    add("B/B0",   l$ref_levels$BB0,   "LBB")
    add("B/Bmsy", l$ref_levels$BBmsy, "LBB")
    add("F/M",    l$ref_levels$FM,    "LBB")
    add("Status (LBB)", l$ref_levels$statut, "LBB")
  }

  ## LBSPR
  s <- results$lbspr
  if (!is.null(s) && !is.null(s$summary) && nrow(s$summary) > 0) {
    last <- s$summary[nrow(s$summary), ]
    add("SPR (last year)", last$SPR, "LBSPR")
    add("F/M (LBSPR)", last$FM, "LBSPR")
    add("Status (SPR)", last$statut, "LBSPR")
  }

  ## YPR (per-recruit reference points)
  y <- results$ypr
  if (!is.null(y) && !is.null(y$reference_points)) {
    rp <- y$reference_points[1, , drop = FALSE]
    if (!is.null(rp$F01))  add("F0.1 (YPR)",  rp$F01,  "YPR")
    if (!is.null(rp$Fmax)) add("Fmax (YPR)",  rp$Fmax, "YPR")
    if (!is.null(rp$E01))  add("E0.1 (YPR)",  rp$E01,  "YPR")
    if (!is.null(rp$Emax)) add("Emax (YPR)",  rp$Emax, "YPR")
  }

  ## TAC
  tc <- results$tac
  if (!is.null(tc) && inherits(tc, "FishStockTAC")) {
    add("TAC (t)", tc$TAC, "TAC")
    add("TAC method", tc$method, "TAC")
  }

  ## Biological rest period
  cl <- results$closure
  if (!is.null(cl) && !is.null(cl$best)) {
    add("Rest window", cl$best$months_label, "Biological rest period")
    add("Rest duration (months)", cl$best$duration, "Biological rest period")
    add("Share of protected score", cl$best$protected, "Biological rest period")
  }

  ## Management sizes
  sa <- results$size_advice
  if (!is.null(sa) && !is.null(sa$reference)) {
    add("Lopt (cm)",   sa$reference$Lopt,   "Technical measures")
    add("Lc_opt (cm)", sa$reference$Lc_opt, "Technical measures")
  }

  if (length(rows) == 0)
    return(data.frame(Indicator = character(0),
                      Value = character(0),
                      Source = character(0)))

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out

}


#==============================================================
# Report generation
#==============================================================

#' Automatic report generation (HTML / PDF / DOCX)
#'
#' Renders the parameterized template supplied with the package
#' (\code{inst/rmarkdown/rapport_stockflow.Rmd}) via
#' \code{rmarkdown::render()}. Only the sections whose results are
#' present in \code{results} are included.
#'
#' @param results Object \code{FishStockResults} (see \code{collect_results()}).
#' @param file Path of the output file (the extension is inferred from the
#'   format if absent).
#' @param format Output format: \code{"html"}, \code{"pdf"} or
#'   \code{"docx"}. PDF requires LaTeX (e.g. \code{tinytex}).
#' @param title,author Title and author of the report.
#' @param template Path of an alternative template (\code{.Rmd}).
#' @param quiet Logical. Silent during rendering.
#' @param ... Additional arguments passed to \code{rmarkdown::render()}.
#'
#' @return Path of the generated file (invisible).
#' @examples
#' \dontrun{
#'   res <- collect_results(growth = growth, mortality = mort, mse = mse,
#'                          meta = list(species = "Penaeus notialis"))
#'   fishstock_report(res, file = "rapport", format = "html")
#' }
#' @export
fishstock_report <- function(results,
                             file = "rapport_stockflow",
                             format = c("html", "pdf", "docx"),
                             title = "Stock assessment - stockflow",
                             author = "stockflow",
                             template = NULL,
                             quiet = TRUE,
                             ...) {

  .check_rmarkdown()

  if (!inherits(results, "FishStockResults"))
    stop("'results' must be a FishStockResults object ",
         "(see collect_results()).", call. = FALSE)

  format <- match.arg(format)


  if (is.null(template))
    template <- system.file("rmarkdown", "rapport_stockflow.Rmd",
                            package = "stockflow")

  if (template == "" || !file.exists(template))
    stop("Template not found. Provide 'template'.", call. = FALSE)


  out_format <- switch(
    format,
    html = "html_document",
    pdf  = "pdf_document",
    docx = "word_document"
  )

  ext <- switch(format, html = ".html", pdf = ".pdf", docx = ".docx")

  if (!grepl("\\.(html|pdf|docx)$", file))
    file <- paste0(file, ext)

  file <- normalizePath(file, mustWork = FALSE)


  ## the template is copied to a temporary folder (avoids writing
  ## into inst/ and collisions of intermediate files)
  tmp_rmd <- tempfile(fileext = ".Rmd")
  file.copy(template, tmp_rmd, overwrite = TRUE)


  out <- tryCatch(

    rmarkdown::render(
      input         = tmp_rmd,
      output_format = out_format,
      output_file   = basename(file),
      output_dir    = dirname(file),
      params        = list(results = results,
                           title   = title,
                           author  = author),
      envir         = new.env(parent = globalenv()),
      quiet         = quiet,
      ...
    ),

    error = function(e)
      stop("Report rendering: ", conditionMessage(e), call. = FALSE)

  )

  message("Report generated: ", out)

  invisible(out)

}