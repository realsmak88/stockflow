############################################################
##
## stockflow
## setup_project.R
##
## Fully initialize the project
##
############################################################

#-----------------------------------------------------------
# Utilities
#-----------------------------------------------------------

create_dir <- function(path) {
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("✓ Folder created: ", path)
  } else {
    message("• Folder already present: ", path)
  }
  
}

create_file <- function(path, content = "") {
  
  if (!file.exists(path)) {
    
    writeLines(content, path)
    
    message("✓ File created: ", path)
    
  } else {
    
    message("• File already present: ", path)
    
  }
  
}

#-----------------------------------------------------------
# Creating folders
#-----------------------------------------------------------

create_project_structure <- function() {
  
  dirs <- c(
    
    "R",
    
    "scripts",
    
    "reports",
    
    "tests",
    
    "inst/extdata",
    
    "data/raw",
    
    "data/processed",
    
    "output",
    
    "output/figures",
    
    "output/tables",
    
    "output/reports"
    
  )
  
  invisible(lapply(dirs, create_dir))
  
}

#-----------------------------------------------------------
# Creating main files
#-----------------------------------------------------------

create_project_files <- function() {
  
  create_file("README.md")
  
  create_file("DESCRIPTION")
  
  create_file("NAMESPACE")
  
  create_file("config.yml")
  
  create_file("_targets.R")
  
  create_file(".gitignore")
  
}

#-----------------------------------------------------------
# Creating R modules
#-----------------------------------------------------------

create_modules <- function() {
  
  modules <- c(
    
    "import",
    
    "validation",
    
    "lfq",
    
    "growth",
    
    "mortality",
    
    "recruitment",
    
    "ypr",
    
    "lbb",
    
    "lbspr",
    
    "diagnostics",
    
    "plots",
    
    "tables",
    
    "scenarios",
    
    "mse"
    
  )
  
  for(m in modules){
    
    create_file(
      
      file.path("R",paste0(m,".R"))
      
    )
    
  }
  
}

#-----------------------------------------------------------
# Summary
#-----------------------------------------------------------

project_summary <- function(){
  
  cat(
    "
====================================

 stockflow

 Project initialized

====================================

"
  )
  
}

#-----------------------------------------------------------
# Execution
#-----------------------------------------------------------

create_project_structure()

create_project_files()

create_modules()

project_summary()
