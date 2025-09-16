# ============================================================
# Report Generation Module for PGLS Analysis Framework
# Creates comprehensive analysis reports and exports results
# ============================================================

source("R/utils.R")

# Generate comprehensive analysis report
generate_comprehensive_report <- function(all_results, config, output_dir) {
  report_file <- file.path(output_dir, "comprehensive_analysis_report.txt")
  
  if (config$output$verbose) cat("\nGenerating comprehensive report...\n")
  
  con <- file(report_file, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  
  # Header
  writeLines(paste(rep("=", 70), collapse = ""), con)
  writeLines("COMPREHENSIVE PGLS/GLS ANALYSIS REPORT", con)
  writeLines(paste(rep("=", 70), collapse = ""), con)
  writeLines("", con)
  
  # Metadata
  writeLines(sprintf("Generated: %s", Sys.time()), con)
  writeLines(sprintf("R version: %s", R.version.string), con)
  writeLines(sprintf("Output directory: %s", normalizePath(output_dir)), con)
  writeLines("", con)
  
  # Configuration summary
  writeLines(paste(rep("-", 70), collapse = ""), con)
  writeLines("ANALYSIS CONFIGURATION", con)
  writeLines(paste(rep("-", 70), collapse = ""), con)
  writeLines(sprintf("Response variables: %s", 
                    paste(config$variables$response_vars, collapse = ", ")), con)
  writeLines(sprintf("Predictor variables: %s", 
                    paste(config$variables$predictor_vars, collapse = ", ")), con)
  writeLines(sprintf("Group variable: %s", 
                    ifelse(is.null(config$variables$group_var), "None", config$variables$group_var)), con)
  writeLines(sprintf("Log transform responses: %s", config$preprocessing$log_transform_response), con)
  writeLines(sprintf("Confidence level: %s", config$statistics$confidence_level), con)
  writeLines("", con)
  
  # Results summary
  writeLines(paste(rep("-", 70), collapse = ""), con)
  writeLines("RESULTS SUMMARY", con)
  writeLines(paste(rep("-", 70), collapse = ""), con)
  
  for (combo_name in names(all_results)) {
    results <- all_results[[combo_name]]$results
    
    if (is.null(results) || nrow(results) == 0) {
      writeLines(sprintf("\n%s: No results", combo_name), con)
      next
    }
    
    writeLines(sprintf("\n%s", combo_name), con)
    writeLines(rep("-", nchar(combo_name)), con)
    
    # Summary statistics
    n_models <- nrow(results)
    n_successful <- sum(!is.na(results$estimate))
    n_significant <- sum(results$p.value < 0.05, na.rm = TRUE)
    
    writeLines(sprintf("Models attempted: %d", n_models), con)
    writeLines(sprintf("Successful fits: %d", n_successful), con)
    writeLines(sprintf("Significant results (p < 0.05): %d", n_significant), con)
    
    # Best model by AICc
    if (n_successful > 0) {
      best_model <- results %>%
        filter(!is.na(AICc)) %>%
        arrange(AICc) %>%
        slice(1)
      
      if (nrow(best_model) > 0) {
        writeLines(sprintf("Best model (by AICc): %s (AICc = %.2f)", 
                          best_model$model[1], best_model$AICc[1]), con)
      }
    }
    
    # Significant effects
    sig_results <- results %>% filter(p.value < 0.05)
    if (nrow(sig_results) > 0) {
      writeLines("\nSignificant effects:", con)
      for (i in 1:nrow(sig_results)) {
        writeLines(sprintf("  %s: estimate = %.4f, p = %.4f %s",
                          sig_results$model[i],
                          sig_results$estimate[i],
                          sig_results$p.value[i],
                          sig_results$stars[i]), con)
      }
    }
  }
  
  # Detailed results table
  writeLines(paste("\n", paste(rep("=", 70), collapse = ""), sep = ""), con)
  writeLines("DETAILED RESULTS TABLE", con)
  writeLines(paste(rep("=", 70), collapse = ""), con)
  
  # Combine all results
  all_results_df <- map_dfr(all_results, ~ .x$results, .id = "analysis")
  
  if (nrow(all_results_df) > 0) {
    # Format for printing
    print_df <- all_results_df %>%
      select(analysis, model, estimate, conf.low, conf.high, 
             p.value, stars, n, AIC, AICc) %>%
      mutate(across(where(is.numeric), ~ round(., 4)))
    
    capture.output(print(print_df, n = Inf), file = con, append = TRUE)
  }
  
  writeLines(paste("\n", paste(rep("=", 70), collapse = ""), sep = ""), con)
  writeLines("END OF REPORT", con)
  writeLines(paste(rep("=", 70), collapse = ""), con)
  
  if (config$output$verbose) {
    cat(sprintf("Report saved to: %s\n", report_file))
  }
}

# Export results to various formats
export_results <- function(results_df, output_file, format = "csv") {
  if (is.null(results_df) || nrow(results_df) == 0) {
    warning("No results to export")
    return(FALSE)
  }
  
  switch(format,
    "csv" = {
      write_csv(results_df, output_file)
    },
    "tsv" = {
      write_tsv(results_df, output_file)
    },
    "xlsx" = {
      if (requireNamespace("writexl", quietly = TRUE)) {
        writexl::write_xlsx(results_df, output_file)
      } else {
        warning("Package 'writexl' not available. Using CSV format instead.")
        write_csv(results_df, gsub("\\.xlsx$", ".csv", output_file))
      }
    },
    {
      write_csv(results_df, output_file)
    }
  )
  
  return(TRUE)
}