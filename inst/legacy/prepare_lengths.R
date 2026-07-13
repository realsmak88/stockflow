############################################################
## [ARCHIVE — inst/legacy/] Obsolete module, DO NOT put back
## into R/: it redefines check_lfq() (and others) already
## present and more developed in R/growth.R. Kept for
## reference only.
############################################################
##
## stockflow
##
## prepare_lengths.R
##
## Preparation of length data
##
############################################################

#============================================================
# Format detection
#============================================================

detect_length_format <- function(data){
  
  nms <- tolower(names(data))
  
  if(all(c("length","count") %in% nms))
    return("aggregated")
  
  if("length" %in% nms)
    return("individual")
  
  if(any(grepl("^20", nms)))
    return("lfq")
  
  stop("Unknown data format.", call.=FALSE)
}


#============================================================
# Cleaning
#============================================================

clean_lengths <- function(data,
                          min_length=0,
                          max_length=Inf){
  
  stopifnot(is.data.frame(data))
  
  names(data) <- tolower(names(data))
  
  data <- data[!is.na(data$length), ]
  
  data <- data[data$length > min_length, ]
  
  data <- data[data$length <= max_length, ]
  
  data <- unique(data)
  
  rownames(data) <- NULL
  
  data
}


#============================================================
# Bin width
#============================================================

choose_bin_width <- function(lengths,
                             method="sturges"){
  
  n <- length(lengths)
  
  if(method=="sturges"){
    
    k <- ceiling(log2(n)+1)
    
    return(diff(range(lengths))/k)
    
  }
  
  if(method=="fd"){
    
    return(2*IQR(lengths)/(n^(1/3)))
    
  }
  
  if(method=="scott"){
    
    return(3.5*sd(lengths)/(n^(1/3)))
    
  }
  
  stop("Unknown method.")
}


#============================================================
# Class boundaries
#============================================================

create_length_breaks <- function(lengths,
                                 bin_width){
  
  seq(
    floor(min(lengths)),
    ceiling(max(lengths))+bin_width,
    by=bin_width
  )
  
}


#============================================================
# Class midpoints
#============================================================

create_mid_lengths <- function(breaks){
  
  head(breaks,-1)+diff(breaks)/2
  
}


#============================================================
# Global histogram
#============================================================

build_histogram <- function(lengths,
                            breaks){
  
  hist(
    lengths,
    breaks=breaks,
    plot=FALSE
  )
  
}


#============================================================
# Monthly aggregation
#============================================================

aggregate_lengths <- function(data,
                              breaks){
  
  if(!all(c("year","month") %in% names(data)))
    stop("year and month columns required.")
  
  data$date <- sprintf(
    "%04d-%02d",
    data$year,
    data$month
  )
  
  dates <- sort(unique(data$date))
  
  mids <- create_mid_lengths(breaks)
  
  freq <- matrix(
    0,
    nrow=length(dates),
    ncol=length(mids)
  )
  
  rownames(freq) <- dates
  colnames(freq) <- mids
  
  for(i in seq_along(dates)){
    
    tmp <- data[data$date==dates[i],]
    
    h <- hist(
      tmp$length,
      breaks=breaks,
      plot=FALSE
    )
    
    freq[i,] <- h$counts
    
  }
  
  list(
    
    dates=dates,
    
    mids=mids,
    
    catch=freq
    
  )
  
}


#============================================================
# LFQ check
#============================================================

check_lfq <- function(lfq){
  
  stopifnot(is.list(lfq))
  
  if(any(rowSums(lfq$catch)==0))
    warning("Months with no catches.")
  
  if(any(colSums(lfq$catch)==0))
    warning("Empty size classes.")
  
  invisible(TRUE)
  
}


#============================================================
# TropFishR-compatible LFQ object
#============================================================

build_lfq_object <- function(data,
                             bin_width=NULL){
  
  data <- clean_lengths(data)
  
  if(is.null(bin_width))
    
    bin_width <- choose_bin_width(
      data$length
    )
  
  breaks <- create_length_breaks(
    data$length,
    bin_width
  )
  
  lfq <- aggregate_lengths(
    data,
    breaks
  )
  
  check_lfq(lfq)
  
  class(lfq) <- c("FishStockLFQ","list")
  
  attr(lfq,"bin_width") <- bin_width
  
  attr(lfq,"breaks") <- breaks
  
  lfq
  
}


#============================================================
# Summary
#============================================================

summary.FishStockLFQ <- function(object,...){
  
  cat("\n")
  
  cat("FishStock LFQ\n")
  
  cat("----------------------\n")
  
  cat("Dates:",length(object$dates),"\n")
  
  cat("Classes:",length(object$mids),"\n")
  
  cat("Observations:",
      sum(object$catch),
      "\n")
  
  cat("Bin width:",
      attr(object,"bin_width"),
      "\n")
  
}


#============================================================
# Printing
#============================================================

print.FishStockLFQ <- function(x,...){
  
  summary(x)
  
  invisible(x)
  
}


#============================================================
# Histogram
#============================================================

plot_length_distribution <- function(data,
                                     bin_width=NULL){
  
  data <- clean_lengths(data)
  
  if(is.null(bin_width))
    
    bin_width <- choose_bin_width(
      data$length
    )
  
  hist(
    
    data$length,
    
    breaks=create_length_breaks(
      data$length,
      bin_width
    ),
    
    main="Length distribution",
    
    xlab="Length",
    
    col="grey80",
    
    border="white"
    
  )
  
}


#============================================================
# LFQ heatmap
#============================================================

plot_lfq_heatmap <- function(lfq){
  
  if(!requireNamespace("ggplot2", quietly=TRUE))
    stop("Install ggplot2.")
  
  df <- as.data.frame(as.table(lfq$catch))
  
  names(df) <- c(
    "Date",
    "Length",
    "Frequency"
  )
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x=Date,
      y=Length,
      fill=Frequency
    )
  ) +
    
    ggplot2::geom_tile() +
    
    ggplot2::scale_fill_viridis_c() +
    
    ggplot2::theme_bw() +
    
    ggplot2::theme(
      
      axis.text.x=
        ggplot2::element_text(
          angle=90,
          hjust=1
        )
      
    )
  
}


#============================================================
# Main function
#============================================================

prepare_lengths <- function(stock,
                            bin_width=NULL){
  
  if(inherits(stock,"FishStockData"))
    
    data <- stock$lengths
  
  else
    
    data <- stock
  
  format <- detect_length_format(data)
  
  message("Detected format: ",format)
  
  if(format!="individual")
    
    stop(
      "Current version: individual data only."
    )
  
  lfq <- build_lfq_object(
    
    data,
    
    bin_width
    
  )
  
  return(lfq)
  
}

############################################################
##
## END prepare_lengths.R
############################################################