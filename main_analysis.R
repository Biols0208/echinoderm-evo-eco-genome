# ============================================================
# Main Analysis Script for PGLS Analysis Framework
# Orchestrates the complete phylogenetic comparative analysis workflow
# ============================================================

# Load all required modules
source("R/config.R")
source("R/data_processing.R")
source("R/model_fitting.R")
source("R/visualization.R")
source("R/reporting.R")

# Main analysis function
run_pgls_analysis <- function(config_file = "analysis_config.yaml") {
  
  # Initialize environment and load configuration
  cat("\n=== Starting PGLS Analysis Framework ===\n")
  config <- load_config(config_file)
  init_packages()
  output_dir <- init_environment(config)
  
  if (config$output$verbose) {
    cat(sprintf("Configuration loaded from: %s\n", config_file))
    cat(sprintf("Output directory: %s\n", output_dir))
  }
  
  # Load and prepare data
  data_loaded <- load_data(config)
  tree_raw <- data_loaded$tree
  data_raw <- data_loaded$data
  
  # Prepare analysis data
  data_prepared <- prepare_analysis_data(data_raw, config)
  
  # Clean tree and create group trees if applicable
  tree_clean <- clean_tree(tree_raw, config$preprocessing$exclude_species)
  group_trees <- NULL
  if (!is.null(config$variables$group_var)) {
    group_trees <- create_group_trees(tree_clean, data_prepared, config)
  }
  
  # Store all analysis results
  all_analysis_results <- list()
  
  # Run analyses for each response-predictor combination
  for (resp_var in config$variables$response_vars) {
    for (pred_var in config$variables$predictor_vars) {
      
      combo_name <- sprintf("%s ~ %s", resp_var, pred_var)
      cat(sprintf("\n%s\n", paste(rep("=", 50), collapse = "")))
      cat(sprintf("Analysis: %s\n", combo_name))
      cat(sprintf("%s\n", paste(rep("=", 50), collapse = "")))
      
      # Run analysis for this combination
      analysis_result <- run_analysis_combination(
        resp_var, pred_var, data_prepared, tree_clean, 
        group_trees, config
      )
      
      if (!is.null(analysis_result$results) && nrow(analysis_result$results) > 0) {
        all_analysis_results[[combo_name]] <- analysis_result
        
        # Save individual results if requested
        if (config$output$save_individual_results) {
          filename <- sprintf("results_%s_vs_%s.%s", 
                             resp_var, pred_var, config$output$export_format)
          export_results(analysis_result$results, 
                        file.path(output_dir, filename),
                        config$output$export_format)
          
          if (config$output$verbose) {
            cat(sprintf("Results saved to: %s\n", filename))
          }
        }
        
        # Generate visualizations if requested
        if (config$visualization$generate_plots) {
          generate_plots_for_combination(analysis_result, combo_name, 
                                        output_dir, config)
        }
        
      } else {
        if (config$output$verbose) {
          cat("No valid results for this combination\n")
        }
      }
    }
  }
  
  # Save combined results and generate reports
  finalize_analysis(all_analysis_results, output_dir, config)
  
  # Final summary
  print_analysis_summary(all_analysis_results, output_dir, config)
  
  return(all_analysis_results)
}

# Run analysis for one response-predictor combination
run_analysis_combination <- function(resp_var, pred_var, data_prepared, 
                                   tree_clean, group_trees, config) {
  
  if (config$output$verbose) {
    cat(sprintf("\n=== Analyzing: %s ~ %s ===\n", resp_var, pred_var))
  }
  
  all_results <- list()
  all_models <- list()
  
  # Determine actual response variable name (with log transform if applicable)
  response_use <- resp_var
  if (config$preprocessing$log_transform_response) {
    log_name <- sprintf("%s_log%g", resp_var, config$preprocessing$log_base)
    if (log_name %in% names(data_prepared)) {
      response_use <- log_name
      if (config$output$verbose) {
        cat(sprintf("Using log%g transformed response: %s\n", 
                   config$preprocessing$log_base, response_use))
      }
    }
  }
  
  # Check data availability
  complete_data <- data_prepared[complete.cases(
    data_prepared[, c(response_use, pred_var, config$variables$species_col)]), ]
  
  if (nrow(complete_data) < config$variables$min_group_size) {
    if (config$output$verbose) {
      cat("Insufficient complete cases for analysis\n")
    }
    return(list(results = NULL, models = NULL))
  }
  
  # Analyze all species
  if (config$output$verbose) cat("\nAnalyzing all species...\n")
  td_all <- prep_tree_data(tree_clean, complete_data, 
                          config$variables$species_col, 
                          config$variables$min_group_size)
  
  if (!is.null(td_all)) {
    for (model_type in names(config$models$model_types)) {
      model_label <- config$models$model_types[[model_type]]
      full_label <- paste("All -", model_label)
      
      if (config$output$verbose) {
        cat(sprintf("  Fitting %s...\n", model_label))
      }
      
      mod <- tryCatch({
        fit_model(td_all$data, td_all$tree, response_use, pred_var, 
                 model_type, config)
      }, error = function(e) {
        warning(sprintf("Failed to fit %s: %s", model_label, e$message))
        NULL
      })
      
      all_models[[full_label]] <- mod
      
      # Extract main model results
      main_results <- extract_model_info(
        mod, td_all$data, full_label, resp_var, pred_var, config
      )
      
      # Add pairwise comparisons if applicable
      if (config$statistics$perform_pairwise && !is.null(mod)) {
        pairwise_res <- perform_pairwise_comparisons(
          mod, pred_var, resp_var, full_label, config
        )
        if (!is.null(pairwise_res)) {
          main_results <- bind_rows(main_results, pairwise_res)
        }
      }
      
      all_results[[full_label]] <- main_results
    }
  }
  
  # Analyze by groups if applicable
  if (!is.null(config$variables$group_var) && !is.null(group_trees)) {
    group_results <- analyze_by_groups(complete_data, group_trees, response_use, pred_var, 
                                      resp_var, config)
    if (!is.null(group_results)) {
      all_results <- c(all_results, group_results$results)
      all_models <- c(all_models, group_results$models)
    }
  }
  
  # Combine and process results
  combined_results <- bind_rows(all_results)
  
  # Add additional columns for visualization
  if (nrow(combined_results) > 0) {
    combined_results <- add_visualization_columns(combined_results, response_use, config)
  }
  
  return(list(results = combined_results, models = all_models))
}

# Analyze by groups helper function
analyze_by_groups <- function(complete_data, group_trees, response_use, pred_var, 
                             resp_var, config) {
  
  group_levels <- levels(complete_data[[config$variables$group_var]])
  group_results <- list()
  group_models <- list()
  
  for (group_name in group_levels) {
    if (config$output$verbose) {
      cat(sprintf("\nAnalyzing group: %s\n", group_name))
    }
    
    # Filter data for this group
    group_data <- complete_data[complete_data[[config$variables$group_var]] == group_name, ]
    
    if (nrow(group_data) < config$variables$min_group_size) {
      if (config$output$verbose) {
        cat("  Insufficient data for this group\n")
      }
      next
    }
    
    # Check if tree exists for this group
    if (!group_name %in% names(group_trees) || is.null(group_trees[[group_name]])) {
      if (config$output$verbose) {
        cat("  No tree available for this group\n")
      }
      next
    }
    
    # Prepare tree and data
    td_group <- prep_tree_data(group_trees[[group_name]], group_data, 
                              config$variables$species_col,
                              config$variables$min_group_size)
    
    if (is.null(td_group)) {
      if (config$output$verbose) {
        cat("  Could not prepare data for this group\n")
      }
      next
    }
    
    # Check predictor variation
    if (length(unique(td_group$data[[pred_var]])) < 2) {
      if (config$output$verbose) {
        cat("  Insufficient predictor variation for this group\n")
      }
      next
    }
    
    # Fit models for this group
    for (model_type in names(config$models$model_types)) {
      model_label <- config$models$model_types[[model_type]]
      full_label <- paste(group_name, "-", model_label)
      
      if (config$output$verbose) {
        cat(sprintf("  Fitting %s...\n", model_label))
      }
      
      mod <- tryCatch({
        fit_model(td_group$data, td_group$tree, response_use, pred_var, 
                 model_type, config)
      }, error = function(e) {
        warning(sprintf("Failed to fit %s for %s: %s", 
                       model_label, group_name, e$message))
        NULL
      })
      
      group_models[[full_label]] <- mod
      group_results[[full_label]] <- extract_model_info(
        mod, td_group$data, full_label, resp_var, pred_var, config
      )
    }
  }
  
  return(list(results = group_results, models = group_models))
}

# Add visualization columns to results
add_visualization_columns <- function(combined_results, response_use, config) {
  combined_results %>%
    mutate(
      # Format lambda
      lambda_fmt = case_when(
        grepl("GLS\\(IID\\)", model) ~ "-",
        grepl("lambda=0", model) ~ "0", 
        grepl("lambda=1", model) ~ "1",
        TRUE ~ ifelse(is.na(lambda), "-", sprintf("%.3f", lambda))
      ),
      # Calculate percentage change for log-transformed responses
      pct_change = if (grepl("log", response_use)) {
        if (config$preprocessing$log_base == 10) {
          (10^estimate - 1) * 100
        } else if (config$preprocessing$log_base == 2) {
          (2^estimate - 1) * 100
        } else {
          (exp(estimate) - 1) * 100
        }
      } else {
        NA_real_
      },
      pct_change_low = if (grepl("log", response_use)) {
        if (config$preprocessing$log_base == 10) {
          (10^conf.low - 1) * 100
        } else if (config$preprocessing$log_base == 2) {
          (2^conf.low - 1) * 100
        } else {
          (exp(conf.low) - 1) * 100
        }
      } else {
        NA_real_
      },
      pct_change_high = if (grepl("log", response_use)) {
        if (config$preprocessing$log_base == 10) {
          (10^conf.high - 1) * 100
        } else if (config$preprocessing$log_base == 2) {
          (2^conf.high - 1) * 100
        } else {
          (exp(conf.high) - 1) * 100
        }
      } else {
        NA_real_
      }
    )
}

# Generate plots for analysis combination
generate_plots_for_combination <- function(analysis_result, combo_name, 
                                         output_dir, config) {
  
  # Generate forest plot
  if (config$visualization$generate_forest_plots) {
    plot_title <- sprintf("Forest Plot: %s", combo_name)
    forest_plot <- create_forest_plot(analysis_result$results, config, plot_title)
    
    if (!is.null(forest_plot)) {
      filename <- sprintf("forest_%s.pdf", 
                         gsub(" ~ ", "_vs_", combo_name))
      
      # Calculate plot height based on number of models
      n_models <- nrow(filter(analysis_result$results, !is.na(estimate)))
      height <- max(config$visualization$forest_plot_min_height, 
                   n_models * config$visualization$forest_plot_height_per_model * 1.2)
      height <- min(height, config$visualization$forest_plot_max_height)
      
      ggsave(
        file.path(output_dir, filename),
        forest_plot,
        width = config$visualization$forest_plot_width,
        height = height,
        device = cairo_pdf
      )
      
      if (config$output$verbose) {
        cat(sprintf("Forest plot saved to: %s\n", filename))
      }
    }
  }
  
  # Generate diagnostic plots for best models
  if (config$visualization$generate_diagnostic_plots && 
      length(analysis_result$models) > 0) {
    
    # Find best model by AICc for each group
    results_by_group <- analysis_result$results %>%
      filter(!is.na(AICc)) %>%
      mutate(group = str_extract(model, "^[^-]+")) %>%
      group_by(group) %>%
      slice_min(AICc, n = 1) %>%
      ungroup()
    
    for (i in 1:nrow(results_by_group)) {
      model_name <- results_by_group$model[i]
      if (model_name %in% names(analysis_result$models)) {
        diag_plot <- create_diagnostic_plots(
          analysis_result$models[[model_name]],
          model_name,
          config
        )
        
        if (!is.null(diag_plot)) {
          filename <- sprintf("diagnostics_%s_%s.pdf",
                             tolower(results_by_group$group[i]),
                             gsub(" ~ ", "_vs_", combo_name))
          
          ggsave(
            file.path(output_dir, filename),
            diag_plot,
            width = config$visualization$diagnostic_plot_width,
            height = config$visualization$diagnostic_plot_height,
            device = cairo_pdf
          )
          
          if (config$output$verbose) {
            cat(sprintf("Diagnostic plot saved to: %s\n", filename))
          }
        }
      }
    }
  }
}

# Finalize analysis with combined results and reports
finalize_analysis <- function(all_analysis_results, output_dir, config) {
  if (length(all_analysis_results) == 0) {
    warning("No analysis results to finalize")
    return()
  }
  
  # Save combined results if requested
  if (config$output$save_combined_results) {
    all_results_df <- map_dfr(all_analysis_results, ~ .x$results, .id = "analysis")
    filename <- sprintf("all_results_combined.%s", config$output$export_format)
    export_results(all_results_df, file.path(output_dir, filename), 
                  config$output$export_format)
    
    if (config$output$verbose) {
      cat(sprintf("\nCombined results saved to: %s\n", filename))
    }
  }
  
  # Generate comprehensive report if requested
  if (config$output$generate_report) {
    generate_comprehensive_report(all_analysis_results, config, output_dir)
  }
}

# Print analysis summary
print_analysis_summary <- function(all_analysis_results, output_dir, config) {
  cat(sprintf("\n%s", paste(rep("=", 70), collapse = "")))
  cat("\nANALYSIS COMPLETE")
  cat(sprintf("\n%s", paste(rep("=", 70), collapse = "")))
  cat(sprintf("\nOutput directory: %s", normalizePath(output_dir)))
  cat(sprintf("\nTotal analyses completed: %d", length(all_analysis_results)))
  
  if (config$visualization$generate_plots) {
    cat("\n\nGenerated visualizations:")
    if (config$visualization$generate_forest_plots) {
      cat("\n- Forest plots for each analysis combination")
    }
    if (config$visualization$generate_diagnostic_plots) {
      cat("\n- Diagnostic plots for best models")
    }
  }
  
  if (config$output$generate_report) {
    cat("\n\nGenerated reports:")
    cat("\n- comprehensive_analysis_report.txt")
  }
  
  cat("\n\nGenerated data files:")
  if (config$output$save_combined_results) {
    cat(sprintf("\n- all_results_combined.%s", config$output$export_format))
  }
  if (config$output$save_individual_results) {
    cat("\n- Individual result files for each combination")
  }
  
  cat(sprintf("\n%s", paste(rep("=", 70), collapse = "")))
  cat("\n")
}

# Run analysis if script is executed directly
if (!interactive()) {
  # Check for command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  config_file <- if (length(args) > 0) args[1] else "analysis_config.yaml"
  
  # Check if config file exists before running analysis
  if (!file.exists(config_file)) {
    cat(sprintf("Error: Configuration file '%s' not found\n", config_file))
    cat("Usage: Rscript main_analysis.R [config_file.yaml]\n")
    cat("Available config files in current directory:\n")
    yaml_files <- list.files(pattern = "\\.yaml$|\\.yml$")
    if (length(yaml_files) > 0) {
      cat(paste("  -", yaml_files), sep = "\n")
    } else {
      cat("  No YAML files found\n")
    }
    quit(status = 1)
  }
  
  # Run the analysis
  cat(sprintf("Starting analysis with config file: %s\n", config_file))
  results <- run_pgls_analysis(config_file)
}
