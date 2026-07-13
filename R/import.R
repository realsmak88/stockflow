############################################################
##
## stockflow
##
## import.R
##
## Module for importing and standardizing fisheries
## data
##
############################################################


#' Importing fisheries data
#'
#' This module allows reading, cleaning and standardizing
#' the data required for analyses:
#'
#' - length frequencies
#' - maturity
#' - catches
#' - effort
#' - parameters
#'
############################################################


#============================================================
# Dependencies
#============================================================

# Packages used:
#
# tidyverse
# readxl
# janitor
# yaml
# here



#============================================================
# Format detection
#============================================================


#' Detect the format of a file
#'
#' @param file file path
#'
#' @return file extension
#'
#' @keywords internal


detect_file_type <- function(file){
  
  if(!file.exists(file)){
    
    stop(
      "File not found: ",
      file,
      call.=FALSE
    )
    
  }
  
  
  ext <- tools::file_ext(file)
  
  ext <- tolower(ext)
  
  
  supported <- c(
    "csv",
    "txt",
    "xlsx",
    "xls",
    "rds"
  )
  
  
  if(!ext %in% supported){
    
    stop(
      "Unsupported format: ",
      ext,
      "\nAccepted formats: ",
      paste(supported,collapse=", "),
      call.=FALSE
    )
    
  }
  
  
  return(ext)
  
}



#============================================================
# CSV separator detection
#============================================================



detect_separator <- function(file){
  
  line <- readLines(
    file,
    n = 1,
    warn = FALSE
  )
  
  separators <- c(
    "," = stringr::str_count(line, ","),
    ";" = stringr::str_count(line, ";"),
    "\t" = stringr::str_count(line, "\t")
  )
  
  sep <- names(separators)[which.max(separators)]
  
  return(sep)
  
}



#============================================================
# Generic reading
#============================================================



read_generic <- function(file){
  
  
  type <- detect_file_type(file)
  
  
  
  data <- switch(
    
    type,
    
    
    csv = {
      
      readr::read_delim(
        file,
        delim=detect_separator(file),
        show_col_types=FALSE
      )
      
    },
    
    
    txt = {
      
      readr::read_delim(
        file,
        delim=detect_separator(file),
        show_col_types=FALSE
      )
      
    },
    
    
    xlsx = {
      
      readxl::read_excel(file)
      
    },
    
    
    xls = {
      
      readxl::read_excel(file)
      
    },
    
    
    rds = {
      
      readRDS(file)
      
    }
    
  )
  
  
  data <- as.data.frame(data)
  
  
  return(data)
  
}



#============================================================
# Standardization of column names
#============================================================



standardize_names <- function(data){
  
  
  names(data) <-
    
    janitor::make_clean_names(
      names(data)
    )
  
  
  
  dictionary <- c(
    
    
    # length
    
    "length"="length",
    
    "len"="length",
    
    "size"="length",
    
    "taille"="length",
    
    "longueur"="length",
    
    "tl"="length",
    
    
    # year
    
    "year"="year",
    
    "annee"="year",
    
    "yr"="year",
    
    
    # month
    
    "month"="month",
    
    "mois"="month",
    
    "mnth"="month",
    
    
    # catch
    
    "catch"="catch",
    
    "capture"="catch",
    
    "captures"="catch",
    
    "landings"="catch",
    
    
    # effort
    
    "effort"="effort",
    
    "fishing_effort"="effort",
    
    
    # maturity
    
    "maturity"="maturity",
    
    "mature"="maturity",
    
    "mat"="maturity"
    
  )
  
  
  
  current <- names(data)
  
  
  
  for(i in seq_along(current)){
    
    
    old <- current[i]
    
    
    if(old %in% names(dictionary)){
      
      
      names(data)[i] <-
        
        dictionary[[old]]
      
    }
    
  }
  
  
  return(data)
  
}



#============================================================
# Type cleaning
#============================================================



clean_types <- function(data){
  
  
  data <- standardize_names(data)
  
  
  
  numeric_vars <- c(
    
    "length",
    
    "catch",
    
    "effort",
    
    "maturity"
    
  )
  
  
  for(v in numeric_vars){
    
    
    if(v %in% names(data)){
      
      
      data[[v]] <-
        
        suppressWarnings(
          as.numeric(data[[v]])
        )
      
      
    }
    
  }
  
  
  
  integer_vars <- c(
    
    "year",
    
    "month"
    
  )
  
  
  for(v in integer_vars){
    
    
    if(v %in% names(data)){
      
      
      data[[v]] <-
        
        suppressWarnings(
          as.integer(data[[v]])
        )
      
    }
    
  }
  
  
  
  return(data)
  
}

############################################################
#
# SPECIALIZED READERS
#
############################################################



#============================================================
# Import of length data
#============================================================


#' Import length frequency data
#'
#' The minimal expected data are:
#'
#' - length
#'
#' Optional:
#'
#' - year
#' - month
#' - sex
#' - station
#'
#' @param file file path
#'
#' @return standardized data.frame
#'
#' @export


read_lengths <- function(file){
  
  
  message("Importing length data...")
  
  
  data <- read_generic(file)
  
  
  data <- clean_types(data)
  
  
  
  required <- c(
    "length"
  )
  
  
  
  missing <- setdiff(
    required,
    names(data)
  )
  
  
  
  if(length(missing)>0){
    
    
    stop(
      "Invalid length data.\n",
      "Missing columns: ",
      paste(
        missing,
        collapse=", "
      ),
      call.=FALSE
    )
    
  }
  
  
  
  data$type <- "length"
  
  
  
  return(data)
  
}






#============================================================
# Import of maturity data
#============================================================


#' Import sexual maturity data
#'
#' Expected columns:
#'
#' - length
#' - maturity
#'
#'
#' @param file file path
#'
#' @return data.frame


#' Import sexual maturity data
#'
#' Reads a maturity file (CSV, TSV or Excel), automatically detects the
#' separator, harmonizes column names and standardizes types.
#'
#' Minimal expected columns: \code{length} and \code{maturity}.
#'
#' @param file Path of the file to import.
#' @return A standardized \code{data.frame}, with a \code{type} column
#'   equal to \code{"maturity"}.
#' @seealso \code{\link{read_lengths}}, \code{\link{import_data}}
#' @examples
#' \dontrun{
#'   mat <- read_maturity("data/raw/maturity.csv")
#' }
#' @export
read_maturity <- function(file){
  
  
  message(
    "Importing maturity data..."
  )
  
  
  
  data <- read_generic(file)
  
  
  data <- clean_types(data)
  
  
  
  required <- c(
    
    "length",
    
    "maturity"
    
  )
  
  
  
  missing <- setdiff(
    required,
    names(data)
  )
  
  
  
  if(length(missing)>0){
    
    
    stop(
      "Invalid maturity data.\n",
      "Missing columns: ",
      paste(
        missing,
        collapse=", "
      ),
      call.=FALSE
    )
    
  }
  
  
  
  data$type <- "maturity"
  
  
  
  return(data)
  
}






#============================================================
# Import of catches
#============================================================


#' Import catch data
#'
#' Expected columns:
#'
#' - year
#' - catch
#'
#'
#' @param file file path
#'
#' @return data.frame


#' Import catch data
#'
#' Reads a catch file (CSV, TSV or Excel) and standardizes it.
#'
#' Minimal expected columns: \code{year} and \code{catch}.
#'
#' @param file Path of the file to import.
#' @return A standardized \code{data.frame}, with a \code{type} column
#'   equal to \code{"catch"}.
#' @seealso \code{\link{read_effort}}, \code{\link{import_data}}
#' @examples
#' \dontrun{
#'   catches <- read_catches("data/raw/catches.csv")
#' }
#' @export
read_catches <- function(file){
  
  
  message(
    "Importing catch data..."
  )
  
  
  
  data <- read_generic(file)
  
  
  data <- clean_types(data)
  
  
  
  required <- c(
    
    "year",
    
    "catch"
    
  )
  
  
  
  missing <- setdiff(
    required,
    names(data)
  )
  
  
  
  if(length(missing)>0){
    
    
    stop(
      "Invalid catch data.\n",
      "Missing columns: ",
      paste(
        missing,
        collapse=", "
      ),
      call.=FALSE
    )
    
  }
  
  
  
  data$type <- "catch"
  
  
  
  return(data)
  
}







#============================================================
# Import of fishing effort
#============================================================


#' Import effort data
#'
#' Expected columns:
#'
#' - year
#' - effort
#'
#'
#' @param file file path
#'
#' @return data.frame



#' Import fishing effort data
#'
#' Reads an effort file (CSV, TSV or Excel) and standardizes it.
#'
#' Minimal expected columns: \code{year} and \code{effort}.
#'
#' @param file Path of the file to import.
#' @return A standardized \code{data.frame}, with a \code{type} column
#'   equal to \code{"effort"}.
#' @seealso \code{\link{read_catches}}, \code{\link{import_data}}
#' @examples
#' \dontrun{
#'   effort <- read_effort("data/raw/effort.csv")
#' }
#' @export
read_effort <- function(file){
  
  
  message(
    "Importing effort data..."
  )
  
  
  
  data <- read_generic(file)
  
  
  
  data <- clean_types(data)
  
  
  
  required <- c(
    
    "year",
    
    "effort"
    
  )
  
  
  
  missing <- setdiff(
    required,
    names(data)
  )
  
  
  
  if(length(missing)>0){
    
    
    stop(
      "Invalid effort data.\n",
      "Missing columns: ",
      paste(
        missing,
        collapse=", "
      ),
      call.=FALSE
    )
    
  }
  
  
  
  data$type <- "effort"
  
  
  
  return(data)
  
}







#============================================================
# Import configuration file
#============================================================



#' Read the project configuration
#'
#' @param file yaml file
#'
#' @return list


#' Read the project configuration file
#'
#' Loads a YAML configuration file describing the species, biological
#' parameters, sampling and management targets.
#'
#' @param file Path of the YAML file (default \code{"config.yml"}).
#' @return A nested list containing the configuration.
#' @examples
#' \dontrun{
#'   cfg <- read_config("config.yml")
#'   cfg$biology$Linf
#' }
#' @export
read_config <- function(file="config.yml"){
  
  
  
  if(!file.exists(file)){
    
    
    stop(
      "Configuration file missing: ",
      file,
      call.=FALSE
    )
    
  }
  
  
  
  config <-
    
    yaml::read_yaml(file)
  
  
  
  return(config)
  
}








#============================================================
# Utility function:
# quick summary of the data
#============================================================



#' Summarize an imported dataset
#'
#' Computes basic descriptive statistics of an imported dataset (number
#' of rows, columns, missing values).
#'
#' @param data Imported \code{data.frame}.
#' @return A summary \code{data.frame}.
#' @examples
#' \dontrun{
#'   summarize_import(lengths)
#' }
#' @export
summarize_import <- function(data){
  
  
  
  if(!is.data.frame(data)){
    
    
    stop(
      "The object must be a data.frame"
    )
    
  }
  
  
  
  summary <- list()
  
  
  
  summary$n <- nrow(data)
  
  
  
  summary$variables <-
    
    names(data)
  
  
  
  if("length" %in% names(data)){
    
    
    summary$length_range <-
      
      range(
        data$length,
        na.rm=TRUE
      )
    
    
  }
  
  
  
  if("year" %in% names(data)){
    
    
    summary$years <-
      
      range(
        data$year,
        na.rm=TRUE
      )
    
    
  }
  
  
  
  if("month" %in% names(data)){
    
    
    summary$months <-
      
      sort(
        unique(
          data$month
        )
      )
    
    
  }
  
  
  
  return(summary)
  
}








#============================================================
# Summary display
#============================================================



print_import_summary <- function(summary){
  
  
  
  cat(
    "\n==============================\n",
    " Import summary\n",
    "==============================\n"
  )
  
  
  
  cat(
    "Number of observations: ",
    summary$n,
    "\n"
  )
  
  
  
  cat(
    "Variables: ",
    paste(
      summary$variables,
      collapse=", "
    ),
    "\n"
  )
  
  
  
  if(!is.null(summary$length_range)){
    
    
    cat(
      "Lengths: ",
      paste(
        summary$length_range,
        collapse=" - "
      ),
      "\n"
    )
    
  }
  
  
  
  if(!is.null(summary$years)){
    
    
    cat(
      "Period: ",
      paste(
        summary$years,
        collapse=" - "
      ),
      "\n"
    )
    
  }
  
  
  
  cat(
    "==============================\n\n"
  )
  
  
  
}

############################################################
#
# STANDARD FISHSTOCKDATA OBJECT
#
############################################################



#============================================================
# Creating the FishStockData object
#============================================================



#' Build a FishStockData object
#'
#' Assembles the various already imported datasets into a single
#' \code{FishStockData} object, consumable by the analysis modules.
#'
#' @param lengths Length \code{data.frame} (optional).
#' @param maturity Maturity \code{data.frame} (optional).
#' @param catches Catch \code{data.frame} (optional).
#' @param effort Effort \code{data.frame} (optional).
#' @param config Configuration list (optional).
#' @return An object of class \code{FishStockData}.
#' @seealso \code{\link{import_data}}
#' @examples
#' \dontrun{
#'   stock <- create_stock_object(lengths = lengths_df)
#' }
#' @export
create_stock_object <- function(
    
  lengths=NULL,
  maturity=NULL,
  catches=NULL,
  effort=NULL,
  config=NULL
  
){
  
  
  stock <- list(
    
    lengths=lengths,
    
    maturity=maturity,
    
    catches=catches,
    
    effort=effort,
    
    config=config,
    
    metadata=list(
      
      created=Sys.time(),
      
      R_version=R.version.string
      
    )
    
  )
  
  
  
  class(stock) <- "FishStockData"
  
  
  
  return(stock)
  
}







#============================================================
# Print object
#============================================================



#' Display of a FishStockData object
#'
#' @param x FishStockData object
#' @param ... additional arguments (ignored).
#'
#' @return \code{x}, invisibly (called for its display side effect).
#' @export


print.FishStockData <- function(x,...){
  
  
  
  cat(
    "\n================================\n"
  )
  
  cat(
    " stockflow object\n"
  )
  
  cat(
    "================================\n\n"
  )
  
  
  
  components <- c(
    
    "lengths",
    
    "maturity",
    
    "catches",
    
    "effort",
    
    "config"
    
  )
  
  
  
  for(i in components){
    
    
    if(!is.null(x[[i]])){
      
      
      cat(
        "\u2713 ",
        i,
        "\n"
      )
      
      
    }else{
      
      
      cat(
        "\u25cb ",
        i,
        "(absent)\n"
      )
      
      
    }
    
    
  }
  
  
  
  cat(
    "\nCreated on: ",
    as.character(
      x$metadata$created
    ),
    "\n"
  )
  
  
  
  cat(
    "\n================================\n"
  )
  
  
}








#============================================================
# Import report
#============================================================



#' Generate an import report
#'
#' @param stock FishStockData object
#'
#' @return list



#' Import summary report
#'
#' Displays in the console a summary of the imported datasets (number of
#' rows, columns, period covered) for quick checking.
#'
#' @param stock \code{FishStockData} object.
#' @return The \code{stock} object, invisibly (called for its display
#'   side effect).
#' @examples
#' \dontrun{
#'   import_report(stock)
#' }
#' @export
import_report <- function(stock){
  
  
  
  report <- list()
  
  
  
  report$date <-
    
    Sys.time()
  
  
  
  report$datasets <- list()
  
  
  
  if(!is.null(stock$lengths)){
    
    
    report$datasets$lengths <-
      
      summarize_import(
        stock$lengths
      )
    
    
  }
  
  
  
  if(!is.null(stock$maturity)){
    
    
    report$datasets$maturity <-
      
      summarize_import(
        stock$maturity
      )
    
    
  }
  
  
  
  if(!is.null(stock$catches)){
    
    
    report$datasets$catches <-
      
      summarize_import(
        stock$catches
      )
    
    
  }
  
  
  
  if(!is.null(stock$effort)){
    
    
    report$datasets$effort <-
      
      summarize_import(
        stock$effort
      )
    
    
  }
  
  
  
  return(report)
  
}







#============================================================
# Save object
#============================================================



#' Save FishStockData object
#'
#' @param stock FishStockData object
#' @param file save path


#' Save a FishStockData object
#'
#' Writes the object in RDS format in order to make the pipeline
#' reproducible.
#'
#' @param stock \code{FishStockData} object.
#' @param file Path of the output RDS file.
#' @return The path of the written file (invisible).
#' @seealso \code{\link{load_stock}}
#' @examples
#' \dontrun{
#'   save_stock(stock, "data/processed/FishStockData.rds")
#' }
#' @export
save_stock <- function(
    
  stock,
  
  file="data/processed/FishStockData.rds"
  
){
  
  
  
  if(!inherits(stock,"FishStockData")){
    
    
    stop(
      "The object is not a FishStockData",
      call.=FALSE
    )
    
  }
  
  
  
  saveRDS(
    
    stock,
    
    file
    
  )
  
  
  
  message(
    "Object saved: ",
    file
  )
  
  
}







#============================================================
# Load object
#============================================================



#' Load a FishStockData object
#'
#' Rereads a \code{FishStockData} object previously saved in RDS format.
#'
#' @param file Path of the RDS file.
#' @return The \code{FishStockData} object.
#' @seealso \code{\link{save_stock}}
#' @examples
#' \dontrun{
#'   stock <- load_stock("data/processed/FishStockData.rds")
#' }
#' @export
load_stock <- function(
    
  file="data/processed/FishStockData.rds"
  
){
  
  
  
  stock <- readRDS(file)
  
  
  
  if(!inherits(stock,"FishStockData")){
    
    
    warning(
      "The loaded object is not a FishStockData"
    )
    
    
  }
  
  
  
  return(stock)
  
}







#============================================================
# Main import function
#============================================================



#' Default configuration
#'
#' Builds the configuration list used when no \code{config.yml} file is
#' provided: unknown species, biological parameters and management
#' targets left empty or at their usual values.
#'
#' @return A configuration list (species, biology, sampling,
#'   management).
#' @seealso \code{\link{read_config}}, \code{\link{import_data}}
#' @examples
#'   cfg <- create_default_config()
#'   names(cfg)
#' @export

create_default_config <- function(){
  
  config <- list(
    
    species=list(
      name="unknown"
    ),
    
    biology=list(
      Linf=NA,
      K=NA,
      M=NA
    ),
    
    management=list(
      SPR_target=0.4
    ),
    
    mse=list(
      years=30,
      simulations=1000
    )
    
  )
  
  return(config)
  
}

#' Import all the data of a stock
#'
#' Main entry function of the pipeline: imports the available datasets
#' (lengths, maturity, catches, effort), reads the configuration and
#' assembles a \code{FishStockData} object.
#'
#' Each source is optional: only those provided are imported.
#'
#' @param lengths Path of the lengths file (optional).
#' @param maturity Path of the maturity file (optional).
#' @param catches Path of the catches file (optional).
#' @param effort Path of the effort file (optional).
#' @param config Path of the YAML configuration file.
#' @return A \code{FishStockData} object (list) grouping the imported
#'   datasets and the configuration.
#' @seealso \code{\link{validate_stock}}, \code{\link{prepare_tropfish}}
#' @examples
#' \dontrun{
#'   stock <- import_data(lengths = "data/raw/lengths.csv",
#'                        catches = "data/raw/catches.csv")
#' }
#' @export
import_data <- function(
    
  lengths=NULL,
  
  maturity=NULL,
  
  catches=NULL,
  
  effort=NULL,
  
  config="config.yml"
  
){
  
  
  
  message(
    "=============================="
  )
  
  message(
    "stockflow - Import"
  )
  
  message(
    "=============================="
  )
  
  
  
  
  
  #--------------------------
  # Lengths
  #--------------------------
  
  
  if(!is.null(lengths)){
    
    
    lengths <-
      
      read_lengths(lengths)
    
    
  }
  
  
  
  
  
  
  #--------------------------
  # maturity
  #--------------------------
  
  
  if(!is.null(maturity)){
    
    
    maturity <-
      
      read_maturity(maturity)
    
    
  }
  
  
  
  
  
  
  #--------------------------
  # Catches
  #--------------------------
  
  
  if(!is.null(catches)){
    
    
    catches <-
      
      read_catches(catches)
    
    
  }
  
  
  
  
  
  
  #--------------------------
  # Effort
  #--------------------------
  
  
  if(!is.null(effort)){
    
    
    effort <-
      
      read_effort(effort)
    
    
  }
  
  
  
  
  
  
  #--------------------------
  # Configuration
  #--------------------------
  
  
  if(!is.null(config)){
    
    if(file.exists(config)){
      
      config <- read_config(config)
      
    }else{
      
      warning(
        "config.yml missing. Creating a default configuration."
      )
      
      config <- create_default_config()
      
    }
    
  }
  
  
  
  
  
  #--------------------------
  # Object creation
  #--------------------------
  
  
  stock <-
    
    create_stock_object(
      
      lengths,
      
      maturity,
      
      catches,
      
      effort,
      
      config
      
    )
  
  
  
  
  
  #--------------------------
  # Report
  #--------------------------
  
  
  report <-
    
    import_report(stock)
  
  
  
  stock$metadata$import_report <-
    
    report
  
  
  
  
  
  message(
    "Import completed successfully."
  )
  
  
  
  return(stock)
  
}







############################################################
#
# END OF MODULE import.R
#
############################################################