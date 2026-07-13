############################################################
## [ARCHIVE — inst/legacy/] Obsolete module, DO NOT put back
## into R/: it redefines check_tropfish(), run_elefan_ga(),
## plot_lfq() with old signatures, conflicting with
## R/growth.R (active versions). Kept for reference.
############################################################
##
## stockflow
##
## prepare_lfq.R
##
## Interface to TropFishR
##
############################################################

#===========================================================
# Checking TropFishR
#===========================================================

check_tropfish <- function(){
  
  if(!requireNamespace("TropFishR", quietly = TRUE))
    stop(
      "The 'TropFishR' package must be installed.",
      call. = FALSE
    )
  
  invisible(TRUE)
}

#===========================================================
# Conversion FishStockLFQ -> lfq TropFishR
#===========================================================

as_tropfish_lfq <- function(lfq,
                            species="Unknown",
                            units="cm"){
  
  check_tropfish()
  
  if(!inherits(lfq,"FishStockLFQ"))
    stop("FishStockLFQ object expected.")
  
  obj <- list(
    
    dates      = as.Date(paste0(lfq$dates,"-15")),
    
    midLengths = as.numeric(lfq$mids),
    
    catch       = as.matrix(lfq$catch),
    
    bin_size    = attr(lfq,"bin_width"),
    
    species     = species,
    
    units       = units
    
  )
  
  class(obj) <- "lfq"
  
  obj
  
}

#===========================================================
# Checking the lfq object
#===========================================================

check_lfq_object <- function(lfq){
  
  stopifnot(inherits(lfq,"lfq"))
  
  stopifnot(is.matrix(lfq$catch))
  
  stopifnot(length(lfq$dates)==nrow(lfq$catch))
  
  stopifnot(length(lfq$midLengths)==ncol(lfq$catch))
  
  invisible(TRUE)
  
}

#===========================================================
# Restructuring
#===========================================================

restructure_lfq <- function(lfq,
                            MA=5,
                            addl.sqrt=FALSE){
  
  check_tropfish()
  
  check_lfq_object(lfq)
  
  TropFishR::lfqRestructure(
    lfq,
    MA=MA,
    addl.sqrt=addl.sqrt
  )
  
}

#===========================================================
# Powell-Wetherall
#===========================================================

run_powell_wetherall <- function(lfq){
  
  check_tropfish()
  
  TropFishR::powell_wetherall(lfq)
  
}

#===========================================================
# ELEFAN (genetic algorithm)
#===========================================================

run_elefan_ga <- function(lfq,
                          Linf_range=NULL,
                          K_range=c(0.1,1)){
  
  check_tropfish()
  
  if(is.null(Linf_range)){
    
    Linf_range <-
      
      c(
        
        max(lfq$midLengths),
        
        max(lfq$midLengths)*1.5
        
      )
    
  }
  
  TropFishR::ELEFAN_GA(
    
    lfq,
    
    low_par=list(
      
      Linf=Linf_range[1],
      
      K=K_range[1]
      
    ),
    
    up_par=list(
      
      Linf=Linf_range[2],
      
      K=K_range[2]
      
    )
    
  )
  
}

#===========================================================
# Diagnostic
#===========================================================

plot_lfq <- function(lfq){
  
  graphics::image(
    
    x=1:nrow(lfq$catch),
    
    y=lfq$midLengths,
    
    z=t(lfq$catch),
    
    xlab="Time",
    
    ylab="Length"
    
  )
  
}

#===========================================================
# Main function
#===========================================================

prepare_lfq <- function(stock,
                        species="Unknown",
                        units="cm"){
  
  raw <- prepare_lengths(stock)
  
  tf <- as_tropfish_lfq(
    
    raw,
    
    species=species,
    
    units=units
    
  )
  
  check_lfq_object(tf)
  
  tf
  
}

############################################################
##
## END
##
############################################################
