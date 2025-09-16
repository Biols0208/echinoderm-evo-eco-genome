# ============================================================
# Utility Functions for PGLS Analysis Framework
# Common helper functions used across modules
# ============================================================

# String concatenation operator
`%.%` <- function(x, y) paste0(x, y)

# Load required packages with error checking
load_packages <- function(packages) {
  for (p in packages) { 
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' not installed. Install with: install.packages('%s')", p, p))
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}

# Initialize analysis environment
init_environment <- function(config) {
  # Set locale for proper character handling
  Sys.setlocale("LC_ALL", "en_US.UTF-8")
  options(encoding = "UTF-8")
  
  # Set random seed if requested
  if (config$reproducibility$set_seed) {
    set.seed(config$reproducibility$random_seed)
    if (config$output$verbose) {
      cat(sprintf("Random seed set to: %d\n", config$reproducibility$random_seed))
    }
  }
  
  # Create output directory
  output_dir <- config$output$output_dir
  if (config$output$output_timestamp) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_dir <- paste0(output_dir, "_", timestamp)
  }
  
  if (nchar(config$output$output_prefix) > 0) {
    output_dir <- paste0(config$output$output_prefix, "_", output_dir)
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  return(output_dir)
}

# Add significance stars to p-values
add_significance_stars <- function(p_values, alpha_levels = c(0.001, 0.01, 0.05, 0.1),
                                  symbols = c("***", "**", "*", ".")) {
  stars <- character(length(p_values))
  
  for (i in seq_along(p_values)) {
    pval <- p_values[i]
    if (!is.na(pval)) {
      for (j in seq_along(alpha_levels)) {
        if (pval < alpha_levels[j]) {
          stars[i] <- symbols[j]
          break
        }
      }
    }
  }
  
  return(stars)
}