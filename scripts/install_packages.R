####################################################
##
## stockflow
## Package installation
##
####################################################

packages <- c(
  
  "TropFishR",
  
  "LBSPR",
  
  "LBB",
  
  "FSA",
  
  "FLSR",
  
  "MSEtool",
  
  "DLMtool",
  
  "tidyverse",
  
  "lubridate",
  
  "janitor",
  
  "targets",
  
  "tarchetypes",
  
  "yaml",
  
  "readxl",
  
  "writexl",
  
  "ggplot2",
  
  "patchwork",
  
  "gt",
  
  "flextable",
  
  "quarto",
  
  "here",
  
  "fs"
  
)

install_if_missing <- function(pkg){
  
  if(!require(pkg,character.only=TRUE))
    
    install.packages(pkg)
  
}

invisible(lapply(packages,install_if_missing))