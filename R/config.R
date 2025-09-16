# ============================================================
# Configuration Management for PGLS Analysis Framework  
# Handles loading and validation of analysis parameters
# ============================================================

source("R/utils.R")

# Load and validate configuration from YAML file
load_config <- function(config_file = "analysis_config.yaml") {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' required for config loading. Install with: install.packages('yaml')")
  }
  
  if (!file.exists(config_file)) {
    stop(sprintf("Configuration file '%s' not found", config_file))
  }
  
  config <- yaml::yaml.load_file(config_file)
  validate_config(config)
  return(config)
}

# Validate configuration parameters
validate_config <- function(config) {
  required_sections <- c("input", "variables", "preprocessing", "models", 
                        "statistics", "visualization", "output")
  
  missing_sections <- setdiff(required_sections, names(config))
  if (length(missing_sections) > 0) {
    stop(sprintf("Missing required config sections: %s", 
                paste(missing_sections, collapse = ", ")))
  }
  
  # Validate file paths
  if (!file.exists(config$input$tree_file)) {
    stop(sprintf("Tree file not found: %s", config$input$tree_file))
  }
  
  if (!file.exists(config$input$data_file)) {
    stop(sprintf("Data file not found: %s", config$input$data_file))
  }
  
  # Validate variable specifications
  if (length(config$variables$response_vars) == 0) {
    stop("At least one response variable must be specified")
  }
  
  if (length(config$variables$predictor_vars) == 0) {
    stop("At least one predictor variable must be specified")
  }
  
  invisible(TRUE)
}

# Initialize package dependencies
init_packages <- function() {
  required_packages <- c("ape", "nlme", "dplyr", "stringr", "ggplot2", 
                        "patchwork", "MuMIn", "tibble", "readr", "phytools", 
                        "purrr", "yaml")
  
  load_packages(required_packages)
}