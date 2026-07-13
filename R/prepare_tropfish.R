#' Construction of an LFQ object compatible with TropFishR
#'
#' Restructures a length data.frame (columns \code{Length}, \code{Year},
#' \code{Month}) into an \code{lfq} object usable by \pkg{TropFishR}.
#'
#' @param stock data.frame of lengths (Length, Year, Month).
#' @param bin_size width of the length classes.
#' @param species species name.
#' @param stock_name stock identifier.
#' @param length_unit length unit (\code{"cm"} or \code{"mm"}).
#' @param Lmin lower bound of the length classes. By default
#'   (\code{NULL}), the smallest observed length is used, with no
#'   truncation. Only set this if selectivity justifies discarding
#'   small classes: a value that is too high truncates small-sized
#'   stocks (e.g. octopus LM 4-25 cm, volutes LCQ 2.6-41 cm).
#' @return \code{lfq} object (TropFishR list).
#' @export
prepare_tropfish <- function(stock,
                             bin_size = 1.0,
                             species = NA,
                             stock_name = NA,
                             length_unit = "cm",
                             Lmin = NULL){
  
  if(!requireNamespace("TropFishR", quietly = TRUE)){
    stop("The TropFishR package must be installed.")
  }
  
  
  ##################################################
  # Extraction
  ##################################################
  
  # Accepts either a direct data.frame,
  # or an object containing $lengths
  if(is.data.frame(stock) || tibble::is_tibble(stock)){
    
    dat <- stock
    
  } else if("lengths" %in% names(stock)){
    
    dat <- stock$lengths
    
  } else {
    
    stop("The provided object must be a table with Length, Year, Month or contain $lengths")
    
  }
  
  
  ##################################################
  # Column check
  ##################################################
  
  required <- c(
    "Length",
    "Year",
    "Month"
  )
  
  
  missing <- setdiff(required, names(dat))
  
  
  if(length(missing)>0){
    
    stop(
      "Missing columns: ",
      paste(missing, collapse=", ")
    )
    
  }
  
  
  ##################################################
  # Cleaning
  ##################################################
  
  # Lower bound: by default the smallest observed length (no
  # truncation). A hard-coded value would truncate small-sized
  # stocks (e.g. octopus LM 4-25 cm, volutes LCQ 2.6-41 cm).
  # If Lmin is provided, truncation is applied; otherwise nothing is
  # touched and TropFishR is left to determine the lower bound itself.
  #
  # WARNING: NEVER inject an Lmin computed here. TropFishR::lfqCreate()
  # applies 'length_unit', so an Lmin computed in the data's unit
  # (e.g. 8 mm) could end up being compared against converted bounds (4.4 cm),
  # hence the "wrong sign in 'by' argument" error.
  dat <- dat[!is.na(dat$Length), ]

  if (!is.null(Lmin))
    dat <- dat[dat$Length >= Lmin, ]
  
  
  ##################################################
  # Date creation
  ##################################################
  
  dat$SampleDate <- as.Date(
    sprintf(
      "%04d-%02d-15",
      dat$Year,
      dat$Month
    )
  )
  
  
  ##################################################
  # TropFishR LFQ
  ##################################################
  
  args <- list(
    data            = dat,
    Lname           = "Length",
    Dname           = "SampleDate",
    bin_size        = bin_size,
    species         = species,
    stock           = stock_name,
    length_unit     = length_unit,
    aggregate_dates = TRUE,
    plot            = FALSE
  )

  # Lmin is only passed if it was explicitly requested.
  if (!is.null(Lmin)) args$Lmin <- Lmin

  lfq <- do.call(TropFishR::lfqCreate, args)


  return(lfq)
}
