############################################################
##
## stockflow
##
## validation.R
##
## Quality control module for fisheries data
##
############################################################

#============================================================
# DESCRIPTION
#============================================================

# This module checks:
#
# - presence of required variables
# - missing values
# - duplicates
# - outlier values
# - biological consistency
#
# The functions produce diagnostic objects
# subsequently used in validate_stock()

#============================================================
# Checking missing values
#============================================================

#' Check missing values
#'
#' @param data data.frame
#'
#' @return list

check_missing <- function(data){
  
  if(!is.data.frame(data)){
    
    stop(
      "L'objet fourni doit \u00eatre un data.frame",
      call.=FALSE
    )
    
  }
  
  missing_count <- sapply(
    data,
    function(x)
      sum(is.na(x))
  )
  
  missing_count <- missing_count[
    missing_count > 0
  ]
  
  result <- list(
    
    status = length(missing_count) == 0,
    
    missing = missing_count
    
  )
  
  
  return(result)
  
}

#============================================================
# Duplicate detection
#============================================================

#' Detect duplicates
#'
#' @param data data.frame
#'
#' @return list

check_duplicates <- function(data){

  duplicates <- sum(
    duplicated(data)
  )
  
  result <- list(
    
    status = duplicates == 0,
    
    n_duplicates = duplicates
    
  )
  
  return(result)
  
}

#============================================================
# Outlier detection
#============================================================

#' Detect extreme values
#'
#' Uses the IQR method
#'
#' @param x numeric vector
#'
#' @return list

detect_outliers <- function(x){
  
  if(!is.numeric(x)){
    
    stop(
      "Le vecteur doit \u00eatre num\u00e9rique",
      call.=FALSE
    )
    
  }
  
  q1 <- stats::quantile(
    x,
    0.25,
    na.rm=TRUE
  )
  
  q3 <- stats::quantile(
    x,
    0.75,
    na.rm=TRUE
  )
  
  iqr <- q3 - q1

  lower <- q1 - 1.5 * iqr
  
  upper <- q3 + 1.5 * iqr
  
  outliers <- which(
    
    x < lower |
      x > upper
    
  )

  result <- list(
    
    status = length(outliers) == 0,
    
    n_outliers = length(outliers),
    
    indices = outliers,
    
    limits = c(
      lower,
      upper
    )
    
  )
  
  return(result)
  
}

#============================================================
# Checking required columns
#============================================================

#' Check required columns
#'
#' @param data data.frame
#' @param required vector of names
#'
#' @return list

check_required_columns <- function(
    
  data,
  
  required
  
){
  
  missing <- setdiff(
    required,
    names(data)
  )
  
  result <- list(
    
    status = length(missing)==0,
    
    missing_columns = missing
    
  )
  
  return(result)
  
}

#============================================================
# Quality score
#============================================================

#' Calculate a quality score
#'
#' @param checks list of checks
#'
#' @return score between 0 and 100

calculate_quality_score <- function(checks){

    status <- unlist(
    lapply(
      checks,
      function(x)
        x$status
    )
  )
  score <- mean(
    status
  ) * 100

  return(
    round(
      score,
      1
    )
  )
  
}

#============================================================
# Summary of a check
#============================================================

#' Display a check
#'
#' @param name test name
#' @param result result

print_check <- function(
    
  name,
  
  result
  
){
  
  if(result$status){
    
    
    cat(
      "\u2713 ",
      name,
      "\n"
    )
    
    
  }else{
    
    
    cat(
      "\u26a0 ",
      name,
      "\n"
    )
    
  }
  
}

############################################################
##
## BIOLOGICAL VALIDATION
##
## Length data
##
############################################################

#============================================================
# Checking biologically possible lengths
#============================================================

#' Check lengths
#'
#' Checks:
#'
#' - missing length
#' - length <= 0
#' - missing values
#' - extreme values
#'
#' @param data length data
#' @param species_config species configuration list (from
#'   \code{\link{read_config}} or \code{\link{create_default_config}});
#'   pass \code{list()} if no configuration is available.
#' @param min_length minimum acceptable length
#' @param max_length maximum acceptable length
#'
#' @return validation object
#' @export

validate_lengths <- function(

  data,

  species_config = list(),

  min_length=0,

  max_length=NULL

){
 
  cat(
    "\n================================\n"
  )
  
  cat(
    "Length validation\n"
  )
  
  cat(
    "================================\n\n"
  )

  checks <- list()

  #------------------------------------------------------------
  # Required columns
  #------------------------------------------------------------
  
  checks$columns <-
    
    check_required_columns(
      
      data,
      
      c("length")
      
    )
  
  #------------------------------------------------------------
  # Missing values
  #------------------------------------------------------------
  
  checks$missing <-
    
    check_missing(data)

  #------------------------------------------------------------
  # Duplicates
  #------------------------------------------------------------
 
  checks$duplicates <-
    
    check_duplicates(data)

  #------------------------------------------------------------
  # Positive values
  #------------------------------------------------------------
 
  invalid_length <-
    
    which(
      
      data$length <= min_length |
        
        is.na(data$length)
      
    )
  
  checks$biological_range <- list(
    
    status =
      
      length(invalid_length)==0,
    
    n_invalid =
      
      length(invalid_length),
    
    indices =
      
      invalid_length
    
  )

  #------------------------------------------------------------
  # Biological upper limit
  #------------------------------------------------------------

  if(!is.null(max_length)){

    too_large <-
      
      which(
        
        data$length > max_length
        
      )

  }else{

    too_large <- integer(0)
    
    
  }
  
  checks$maximum_length <- list(
    
    status =
      
      length(too_large)==0,
    
    n_invalid =
      
      length(too_large),
    
    indices =
      
      too_large
    
  )
  
  #------------------------------------------------------------
  # Statistical extreme values
  #------------------------------------------------------------

  checks$outliers <-
    
    detect_outliers(
      
      data$length
      
    )

  #------------------------------------------------------------
  # Month
  #------------------------------------------------------------

  if("month" %in% names(data)){

    invalid_month <-
      
      which(
        
        data$month < 1 |
          
          data$month > 12
        
      )

    checks$months <- list(
      
      status =
        
        length(invalid_month)==0,
      
      n_invalid =
        
        length(invalid_month),
      
      indices =
        
        invalid_month
      
    )

  }
  
  #------------------------------------------------------------
  # Years
  #------------------------------------------------------------
  
  
  
  if("year" %in% names(data)){
    
    
    
    invalid_year <-
      
      which(
        
        data$year < 1900 |
          
          data$year > as.numeric(
            
            format(Sys.Date(),"%Y")
            
          )
        
      )
    
    
    
    checks$years <- list(
      
      status =
        
        length(invalid_year)==0,
      
      n_invalid =
        
        length(invalid_year),
      
      indices =
        
        invalid_year
      
    )
    
    
  }
  
  
  
  
  
  
  
  
  #------------------------------------------------------------
  # Overall score
  #------------------------------------------------------------
  
  
  
  score <-
    
    calculate_quality_score(
      
      checks
      
    )
  
  
  
  
  
  result <- list(
    
    dataset="lengths",
    
    checks=checks,
    
    quality_score=score
    
  )
  
  
  
  class(result) <-
    
    "FishStockValidation"
  
  
  
  return(result)
  
}







#============================================================
# Length validation display
#============================================================

#' Display a FishStockValidation object
#' @param x FishStockValidation object
#' @param ... ignored
#' @export
print.FishStockValidation <- function(x,...){
  
  
  
  cat(
    "\n================================\n"
  )
  
  cat(
    "stockflow Validation\n"
  )
  
  cat(
    "================================\n\n"
  )
  
  
  
  cat(
    "Dataset : ",
    x$dataset,
    "\n\n"
  )
  
  for(i in names(x$checks)){
    
    print_check(
      
      i,
      
      x$checks[[i]]
      
    )
    
  }
  
  cat(
    "\nOverall quality : ",
    x$quality_score,
    "%\n"
  )
  
  cat(
    "\n================================\n"
  )
  
}

############################################################
##
## VALIDATION OF FISHING AND REPRODUCTION DATA
##
############################################################

#============================================================
#
# MATURITY VALIDATION
#
#============================================================

#' Validate maturity data
#'
#' Checks:
#'
#' - length and maturity columns
#' - missing values
#' - maturity between 0 and 1
#' - length/maturity consistency
#'
#' @param data maturity data.frame
#'
#' @return FishStockValidation

#' Validate maturity data
#'
#' Checks the quality of the maturity data: required columns, missing
#' values, duplicates, consistency of maturity codes and lengths.
#'
#' @param data maturity \code{data.frame} (columns \code{length} and
#'   \code{maturity}).
#' @return A \code{FishStockValidation} object containing the checks
#'   performed and a quality score.
#' @seealso \code{\link{validate_lengths}}, \code{\link{validate_stock}}
#' @examples
#' \dontrun{
#'   validate_maturity(maturity_df)
#' }
#' @export
validate_maturity <- function(data){
  
  
  cat(
    "\n================================\n"
  )
  
  cat(
    "Maturity validation\n"
  )
  
  cat(
    "================================\n\n"
  )
  
  checks <- list()
  
  #------------------------------------------------------------
  # Columns
  #------------------------------------------------------------
  
  
  checks$columns <-
    
    check_required_columns(
      
      data,
      
      c(
        "length",
        "maturity"
      )
      
    )
  
  #------------------------------------------------------------
  # Missing values
  #------------------------------------------------------------
  
  checks$missing <-
    
    check_missing(data)
  
  #------------------------------------------------------------
  # Maturity values between 0 and 1
  #------------------------------------------------------------
  
  
  invalid <-
    
    which(
      
      data$maturity < 0 |
        
        data$maturity > 1 |
        
        is.na(data$maturity)
      
    )
  
  
  
  
  checks$maturity_range <- list(
    
    status =
      length(invalid)==0,
    
    n_invalid =
      length(invalid),
    
    indices =
      invalid
    
  )
  
  
  
  
  
  
  
  #------------------------------------------------------------
  # Biological consistency
  #------------------------------------------------------------
  
  
  negative_length <-
    
    which(
      
      data$length <= 0
      
    )
  
  
  
  checks$length_valid <- list(
    
    status =
      length(negative_length)==0,
    
    n_invalid =
      length(negative_length),
    
    indices =
      negative_length
    
  )
  
  
  
  
  
  
  score <-
    
    calculate_quality_score(
      checks
    )
  
  
  
  result <- list(
    
    dataset="maturity",
    
    checks=checks,
    
    quality_score=score
    
  )
  
  class(result) <-
    
    "FishStockValidation"
  
  return(result)
  
}

#============================================================
#
# CATCH VALIDATION
#
#============================================================

#' Validate catch data
#'
#' Checks:
#'
#' - year
#' - positive catch
#' - missing values
#'
#'
#' @param data catch data.frame
#'
#' @return FishStockValidation

#' Validate catch data
#'
#' Checks the quality of the catch series: required columns, missing
#' or negative values, duplicate years, outlier values.
#'
#' @param data catch \code{data.frame} (columns \code{year} and
#'   \code{catch}).
#' @return A \code{FishStockValidation} object.
#' @seealso \code{\link{validate_effort}}, \code{\link{validate_stock}}
#' @examples
#' \dontrun{
#'   validate_catches(catches_df)
#' }
#' @export
validate_catches <- function(data){
  
  cat(
    "\n================================\n"
  )
  
  cat(
    "Catch validation\n"
  )
  
  cat(
    "================================\n\n"
  )
  
  checks <- list()
  
  checks$columns <-
    
    check_required_columns(
      
      data,
      
      c(
        "year",
        "catch"
      )
      
    )
  
  checks$missing <-
    
    check_missing(data)
  
  #------------------------------------------------------------
  # Positive catch
  #------------------------------------------------------------
  
  negative <-
    
    which(
      
      data$catch < 0
      
    )
  
  checks$catch_range <- list(
    
    status =
      length(negative)==0,
    
    n_invalid =
      length(negative),
    
    indices =
      negative
    
  )
  
  #------------------------------------------------------------
  # Years
  #------------------------------------------------------------
  
  if("year" %in% names(data)){
    
    
    bad_year <-
      
      which(
        
        data$year < 1900 |
          
          data$year >
          as.numeric(format(Sys.Date(),"%Y"))
        
      )
    
    
    
    checks$years <- list(
      
      status =
        length(bad_year)==0,
      
      n_invalid =
        length(bad_year),
      
      indices =
        bad_year
      
    )
    
    
  }
  
  
  
  
  
  score <-
    
    calculate_quality_score(
      checks
    )
  
  
  
  
  result <- list(
    
    dataset="catches",
    
    checks=checks,
    
    quality_score=score
    
  )
  
  
  
  class(result) <-
    
    "FishStockValidation"
  
  
  
  return(result)
  
  
}







#============================================================
#
# EFFORT VALIDATION
#
#============================================================



#' Validate effort data
#'
#' Checks:
#'
#' - positive effort
#' - valid year
#'
#'
#' @param data effort data.frame
#'
#' @return FishStockValidation



#' Validate fishing effort data
#'
#' Checks the quality of the effort series: required columns, missing,
#' null, or negative values, duplicates, outlier values.
#'
#' @param data effort \code{data.frame} (columns \code{year} and
#'   \code{effort}).
#' @return A \code{FishStockValidation} object.
#' @seealso \code{\link{validate_catches}}, \code{\link{validate_stock}}
#' @examples
#' \dontrun{
#'   validate_effort(effort_df)
#' }
#' @export
validate_effort <- function(data){
  
  
  
  cat(
    "\n================================\n"
  )
  
  cat(
    "Effort validation\n"
  )
  
  cat(
    "================================\n\n"
  )
  
  
  
  checks <- list()
  
  
  
  
  
  checks$columns <-
    
    check_required_columns(
      
      data,
      
      c(
        "year",
        "effort"
      )
      
    )
  
  
  
  
  
  
  
  checks$missing <-
    
    check_missing(data)
  
  
  
  
  
  
  
  #------------------------------------------------------------
  # Positive effort
  #------------------------------------------------------------
  
  
  bad_effort <-
    
    which(
      
      data$effort <= 0
      
    )
  
  
  
  checks$effort_range <- list(
    
    status =
      length(bad_effort)==0,
    
    n_invalid =
      length(bad_effort),
    
    indices =
      bad_effort
    
  )
  
  
  
  
  
  
  score <-
    
    calculate_quality_score(
      checks
    )
  
  
  
  
  
  result <- list(
    
    dataset="effort",
    
    checks=checks,
    
    quality_score=score
    
  )
  
  
  
  class(result) <-
    
    "FishStockValidation"
  
  
  
  return(result)
  
  
}







############################################################
#
# END BLOCK 3/4
#
############################################################

############################################################
##
## VALIDATION ORCHESTRATOR
##
############################################################

#============================================================
# Creating the validation object
#============================================================

create_validation_object <- function(results){
  
  scores <- vapply(
    results,
    function(x) x$quality_score,
    numeric(1)
  )
  
  global_score <- round(mean(scores), 1)
  
  structure(
    list(
      date = Sys.time(),
      global_score = global_score,
      datasets = results
    ),
    class = "FishStockValidationReport"
  )
  
}


#============================================================
# Full project validation
#============================================================

#' Validate all data of a stock
#'
#' Chains together the validations of all datasets present in the
#' \code{FishStockData} object (lengths, maturity, catches, effort) and
#' produces an overall report with a quality score.
#'
#' @param stock \code{FishStockData} object (see \code{\link{import_data}}).
#' @return A \code{FishStockValidationReport} object, with \code{print} and
#'   \code{summary} methods.
#' @seealso \code{\link{validate_lengths}}, \code{\link{validate_catches}}
#' @examples
#' \dontrun{
#'   rapport <- validate_stock(stock)
#'   summary(rapport)
#' }
#' @export
validate_stock <- function(stock){
  
  if(!inherits(stock, "FishStockData")){
    stop(
      "L'objet doit \u00eatre de classe 'FishStockData'.",
      call. = FALSE
    )
  }
  
  cat("\n")
  cat("=========================================\n")
  cat(" stockflow - Full validation\n")
  cat("=========================================\n")
  
  results <- list()
  
  if(!is.null(stock$lengths)){
    results$lengths <- validate_lengths(stock$lengths)
  }
  
  if(!is.null(stock$maturity)){
    results$maturity <- validate_maturity(stock$maturity)
  }
  
  if(!is.null(stock$catches)){
    results$catches <- validate_catches(stock$catches)
  }
  
  if(!is.null(stock$effort)){
    results$effort <- validate_effort(stock$effort)
  }
  
  report <- create_validation_object(results)
  
  cat("\nValidation complete.\n")
  cat("Overall score :", report$global_score, "%\n\n")
  
  return(report)
  
}


#============================================================
# print() method
#============================================================

#' @export
print.FishStockValidationReport <- function(x, ...){
  
  cat("\n")
  cat("=========================================\n")
  cat(" stockflow - Validation report\n")
  cat("=========================================\n\n")
  
  cat("Date :", format(x$date), "\n")
  cat("Overall score :", x$global_score, "%\n\n")
  
  cat("Datasets validated\n")
  cat("-----------------------\n")
  
  for(n in names(x$datasets)){
    
    score <- x$datasets[[n]]$quality_score
    
    cat(
      sprintf("%-12s : %5.1f %%\n",
              n,
              score)
    )
    
  }
  
  cat("\n=========================================\n")
  
  invisible(x)
  
}


#============================================================
# summary() method
#============================================================

#' @export
summary.FishStockValidationReport <- function(object, ...){
  
  out <- data.frame(
    
    dataset = names(object$datasets),
    
    quality_score =
      vapply(
        object$datasets,
        function(x) x$quality_score,
        numeric(1)
      )
    
  )
  
  return(out)
  
}


############################################################
##
## END validation.R
##
############################################################