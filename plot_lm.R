#!/usr/bin/env Rscript

# Complete device cleanup and options setup
options(warn = 1)
suppressMessages({
  while (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }
  options(device = "pdf", bitmapType = "cairo")
})

# Global constants (keeping most functions from working version, focusing on plot fixes)
DEFAULT_OPTS <- list(
  input = "test2.csv",
  xcol = "Size_Value_cm", 
  ycol = "genomesize",
  groupcol = "group",
  by_group = FALSE,
  log = "none",
  log_base = "10",
  outfile = NA_character_,
  width = 7,
  height = 5.5,
  dpi = 300,
  title = NA_character_,
  log_label = "pow10",
  log_step = 1.0,
  log_digits = 1,
  pgls = FALSE,
  pgls_correlation = "brownian",
  tree_file = NA_character_,
  species_col = "species",
  check_phylosig = TRUE,
  make_ultrametric = TRUE,
  help = FALSE
)

# [Include all utility functions from the working version - keeping them the same]
parse_log_base <- function(x) {
  if (tolower(x) %in% c("e", "ln")) return(exp(1))
  val <- suppressWarnings(as.numeric(x))
  if (!is.na(val) && val > 0) return(val)
  return(10)
}

validate_columns <- function(df, xcol, ycol, groupcol = NULL) {
  missing_cols <- c()
  if (!(xcol %in% names(df))) missing_cols <- c(missing_cols, xcol)
  if (!(ycol %in% names(df))) missing_cols <- c(missing_cols, ycol)
  if (!is.null(groupcol) && !(groupcol %in% names(df))) {
    warning(sprintf("Group column '%s' not found, ignoring grouping", groupcol))
    groupcol <- NULL
  }
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns: %s", paste(missing_cols, collapse = ", ")))
  }
  return(groupcol)
}

validate_data_quality <- function(df, xcol, ycol) {
  x_vals <- as.numeric(df[[xcol]])
  y_vals <- as.numeric(df[[ycol]])
  
  x_finite <- sum(is.finite(x_vals), na.rm = TRUE)
  y_finite <- sum(is.finite(y_vals), na.rm = TRUE)
  
  if (x_finite < 3 || y_finite < 3) {
    stop("Insufficient valid data points for regression (need at least 3)")
  }
  
  if (length(unique(x_vals[is.finite(x_vals)])) < 2) {
    stop("X variable has insufficient variation")
  }
  if (length(unique(y_vals[is.finite(y_vals)])) < 2) {
    stop("Y variable has insufficient variation") 
  }
  
  cat(sprintf("Data quality check passed: %d valid X values, %d valid Y values\n", 
              x_finite, y_finite))
}

clean_numeric_data <- function(df, xcol, ycol) {
  suppressWarnings({
    df[[xcol]] <- as.numeric(df[[xcol]])
    df[[ycol]] <- as.numeric(df[[ycol]])
  })
  
  if (anyNA(df[[xcol]]) || anyNA(df[[ycol]])) {
    warning("Non-numeric values found in x or y columns, will be removed")
  }
  
  complete_idx <- stats::complete.cases(df[[xcol]], df[[ycol]])
  cleaned <- df[complete_idx, , drop = FALSE]
  
  if (nrow(cleaned) == 0) {
    stop("No complete cases remaining after cleaning")
  }
  
  return(cleaned)
}

filter_positive_values <- function(df, xcol, ycol, to_log_x, to_log_y) {
  n_original <- nrow(df)
  
  if (to_log_x) {
    positive_x <- !is.na(df[[xcol]]) & df[[xcol]] > 0 & is.finite(df[[xcol]])
    df <- df[positive_x, , drop = FALSE]
    if (sum(!positive_x) > 0) {
      warning(sprintf("Removed %d rows with x <= 0, NA, or infinite for log transformation", sum(!positive_x)))
    }
  }
  
  if (to_log_y) {
    positive_y <- !is.na(df[[ycol]]) & df[[ycol]] > 0 & is.finite(df[[ycol]])
    df <- df[positive_y, , drop = FALSE]
    if (sum(!positive_y) > 0) {
      warning(sprintf("Removed %d rows with y <= 0, NA, or infinite for log transformation", sum(!positive_y)))
    }
  }
  
  if (nrow(df) == 0) {
    stop("No valid data points remaining after filtering")
  }
  
  cat(sprintf("Filtered data: %d -> %d rows for log transformation\n", n_original, nrow(df)))
  return(df)
}

# Enhanced argument parsing (same as working version)
parse_args <- function(args) {
  opts <- DEFAULT_OPTS
  
  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    
    if (arg %in% c("-h", "--help")) {
      opts$help <- TRUE
      i <- i + 1
      next
    }
    
    if (!grepl("^--", arg)) {
      warning(sprintf("Ignoring positional argument: %s", arg))
      i <- i + 1
      next
    }
    
    if (grepl("=", arg)) {
      kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
      key <- kv[[1]]
      val <- paste(kv[-1], collapse = "=")
    } else {
      key <- sub("^--", "", arg)
      if (i == length(args) || grepl("^--", args[[i+1]])) {
        val <- TRUE
      } else {
        val <- args[[i+1]]
        i <- i + 1
      }
    }
    
    key_norm <- gsub("-", "_", key)
    
    key_mapping <- list(
      "in" = "input", "x" = "xcol", "y" = "ycol", "group" = "groupcol",
      "out" = "outfile", "log_axis" = "log", "log_base" = "log_base",
      "tree" = "tree_file", "species" = "species_col", "correlation" = "pgls_correlation"
    )
    
    if (key_norm %in% names(key_mapping)) {
      key_norm <- key_mapping[[key_norm]]
    }
    
    if (key_norm %in% names(opts)) {
      opts[[key_norm]] <- convert_value(val, opts[[key_norm]])
    } else {
      warning(sprintf("Unknown option: --%s", key))
    }
    
    i <- i + 1
  }
  
  valid_corr <- c("brownian", "pagel", "grafen", "auto")
  if (!(opts$pgls_correlation %in% valid_corr)) {
    warning(sprintf("Invalid pgls_correlation '%s', using 'brownian'", opts$pgls_correlation))
    opts$pgls_correlation <- "brownian"
  }
  
  if (is.na(opts$outfile) || nchar(opts$outfile) == 0) {
    opts$outfile <- generate_output_filename(opts)
  } else if (!grepl("\\.pdf$", tolower(opts$outfile))) {
    opts$outfile <- paste0(opts$outfile, ".pdf")
  }
  
  return(opts)
}

convert_value <- function(val, template) {
  if (is.logical(template)) {
    if (is.logical(val)) return(val)
    return(tolower(as.character(val)) %in% c("1", "true", "t", "yes", "y"))
  } else if (is.numeric(template)) {
    return(suppressWarnings(as.numeric(val)))
  } else {
    return(as.character(val))
  }
}

generate_output_filename <- function(opts) {
  suffix <- if (tolower(opts$log) %in% c("x", "y", "xy", "both")) {
    paste0("_log", opts$log)
  } else {
    ""
  }
  bygroup_suffix <- if (isTRUE(opts$by_group)) "_bygroup" else ""
  pgls_suffix <- if (isTRUE(opts$pgls)) {
    paste0("_pgls_", opts$pgls_correlation)
  } else {
    ""
  }
  return(sprintf("lm_%s_vs_%s%s%s%s.pdf", opts$ycol, opts$xcol, suffix, bygroup_suffix, pgls_suffix))
}

# Tree preparation (same as working version)
prepare_tree_for_pgls <- function(tree, make_ultrametric = TRUE) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("ape package required for tree operations")
  }
  
  cat("===== Tree Preparation =====\n")
  cat(sprintf("Original tree: %d tips, ultrametric = %s\n", 
              length(tree$tip.label), ape::is.ultrametric(tree)))
  
  tree_fixed <- tree
  
  if (make_ultrametric && !ape::is.ultrametric(tree_fixed)) {
    cat("Making tree ultrametric...\n")
    tryCatch({
      tree_fixed <- ape::chronos(tree_fixed, quiet = TRUE)
      cat("Successfully converted to ultrametric tree\n")
    }, error = function(e) {
      warning(sprintf("Failed to make tree ultrametric with chronos: %s", e$message))
      tryCatch({
        tree_fixed <- ape::compute.brlen(tree_fixed, method = "Grafen")
        cat("Used Grafen method to add branch lengths\n")
      }, error = function(e2) {
        warning(sprintf("All ultrametric conversion methods failed: %s", e2$message))
        cat("Proceeding with original tree (may cause PGLS issues)\n")
        tree_fixed <- tree
      })
    })
  }
  
  # Ensure positive branch lengths
  if (any(tree_fixed$edge.length <= 0, na.rm = TRUE)) {
    cat("Fixing non-positive branch lengths...\n")
    tree_fixed$edge.length[tree_fixed$edge.length <= 0] <- 1e-6
  }
  
  # Check for NAs in branch lengths
  if (any(is.na(tree_fixed$edge.length))) {
    cat("Fixing NA branch lengths...\n")
    tree_fixed$edge.length[is.na(tree_fixed$edge.length)] <- 1e-6
  }
  
  cat(sprintf("Final tree: ultrametric = %s, min branch length = %.2e\n", 
              ape::is.ultrametric(tree_fixed), min(tree_fixed$edge.length, na.rm = TRUE)))
  
  return(tree_fixed)
}

# PGLS data preparation (same as working version)
prepare_pgls_data <- function(dfc, tree, species_col) {
  tryCatch({
    if (!(species_col %in% names(dfc))) {
      warning(sprintf("Species column '%s' not found in data", species_col))
      return(NULL)
    }
    
    if (is.null(tree) || !inherits(tree, "phylo")) {
      warning("Invalid phylogenetic tree")
      return(NULL)
    }
    
    # Clean species names in data
    data_species <- as.character(dfc[[species_col]])
    data_species <- gsub("^\\s+|\\s+$", "", data_species)
    data_species <- gsub("\\s+", "_", data_species)
    
    # Clean species names in tree
    tree_species <- tree$tip.label
    tree_species <- gsub("^\\s+|\\s+$", "", tree_species)
    
    cat(sprintf("Data species (first 5): %s\n", paste(head(data_species, 5), collapse = ", ")))
    cat(sprintf("Tree species (first 5): %s\n", paste(head(tree_species, 5), collapse = ", ")))
    
    # Find common species
    common_species <- intersect(data_species, tree_species)
    if (length(common_species) < 3) {
      warning(sprintf("Insufficient overlap between data and tree. Found %d common species (need ≥3).", 
                     length(common_species)))
      return(NULL)
    }
    
    cat(sprintf("Data-tree overlap: %d species in common\n", length(common_species)))
    cat(sprintf("Species only in data: %d\n", length(setdiff(data_species, tree_species))))
    cat(sprintf("Species only in tree: %d\n", length(setdiff(tree_species, data_species))))
    
    # Update species names in original data
    dfc[[species_col]] <- data_species
    
    # Filter data to common species
    dfc_filtered <- dfc[data_species %in% common_species, , drop = FALSE]
    
    # Handle duplicates
    duplicated_species <- duplicated(dfc_filtered[[species_col]])
    if (any(duplicated_species)) {
      n_dup <- sum(duplicated_species)
      warning(sprintf("Removing %d duplicate species entries (keeping first occurrence)", n_dup))
      dfc_filtered <- dfc_filtered[!duplicated_species, , drop = FALSE]
    }
    
    # CRITICAL: Set row names to species names for PGLS
    rownames(dfc_filtered) <- as.character(dfc_filtered[[species_col]])
    
    if (nrow(dfc_filtered) < 3) {
      warning("Insufficient data for PGLS after filtering")
      return(NULL)
    }
    
    cat(sprintf("Prepared PGLS dataset: %d species\n", nrow(dfc_filtered)))
    return(as.data.frame(dfc_filtered))
    
  }, error = function(e) {
    warning(sprintf("Error preparing PGLS data: %s", e$message))
    return(NULL)
  })
}

# WORKING PGLS model fitting (same as working version)
fit_pgls_model_working <- function(data, formula_str, tree, correlation_type = "auto") {
  if (!requireNamespace("ape", quietly = TRUE) || !requireNamespace("nlme", quietly = TRUE)) {
    stop("PGLS requires 'ape' and 'nlme' packages")
  }
  
  # Ensure data is properly formatted
  data <- as.data.frame(data)
  species_names <- rownames(data)
  
  if (length(species_names) == 0 || any(is.na(species_names)) || any(species_names == "")) {
    stop("Invalid species names in data row names")
  }
  
  # Prune tree to match data
  cat(sprintf("Pruning tree from %d to %d tips\n", length(tree$tip.label), length(species_names)))
  tree_pruned <- ape::keep.tip(tree, species_names)
  
  # CRITICAL: Ensure tree and data are in exact same order
  data_ordered <- data[tree_pruned$tip.label, , drop = FALSE]
  
  cat(sprintf("Final data for PGLS: %d species\n", nrow(data_ordered)))
  
  models <- list()
  
  # Define correlation structures to try (based on debug success)
  if (correlation_type == "auto") {
    corr_types <- c("brownian", "pagel", "grafen")
  } else {
    corr_types <- correlation_type
  }
  
  for (corr_type in corr_types) {
    cat(sprintf("Attempting %s correlation model...\n", corr_type))
    
    tryCatch({
      # Based on debug findings: basic corBrownian works!
      cor_struct <- switch(corr_type,
        "brownian" = {
          # This is what worked in the debug session
          ape::corBrownian(phy = tree_pruned)
        },
        "pagel" = {
          ape::corPagel(value = 0.5, phy = tree_pruned, fixed = FALSE)
        },
        "grafen" = {
          ape::corGrafen(value = 1, phy = tree_pruned, fixed = FALSE)
        }
      )
      
      if (!is.null(cor_struct)) {
        # Use the same settings that worked in debug
        control_params <- nlme::glsControl(
          maxIter = 200,
          msMaxIter = 400,
          tolerance = 1e-6,
          msTol = 1e-7,
          msVerbose = FALSE
        )
        
        # The key insight: warnings are OK, we just need to catch actual errors
        fit <- nlme::gls(stats::as.formula(formula_str), 
                        data = data_ordered, 
                        correlation = cor_struct, 
                        method = "ML",
                        control = control_params,
                        na.action = na.omit)
        
        # Check if model actually fitted (even with warnings)
        if (!is.null(fit) && !is.null(fit$coefficients)) {
          models[[corr_type]] <- fit
          cat(sprintf("✓ Successfully fitted %s correlation model\n", corr_type))
        } else {
          cat(sprintf("✗ %s model returned NULL coefficients\n", corr_type))
        }
      }
    }, error = function(e) {
      cat(sprintf("✗ Failed to fit %s model: %s\n", corr_type, e$message))
    })
  }
  
  if (length(models) == 0) {
    cat("All PGLS models failed. Attempting fallback OLS model...\n")
    fallback_model <- lm(stats::as.formula(formula_str), data = data_ordered)
    warning("All PGLS models failed to fit. Returning OLS model instead.")
    return(fallback_model)
  }
  
  # Model selection if multiple models fitted
  if (length(models) > 1) {
    aics <- sapply(models, AIC)
    best_model <- models[[which.min(aics)]]
    cat(sprintf("Model selection results:\n"))
    for (i in seq_along(aics)) {
      cat(sprintf("  %s: AIC = %.2f %s\n", 
                  names(aics)[i], aics[i], 
                  ifelse(i == which.min(aics), "(BEST)", "")))
    }
    cat(sprintf("Selected best model: %s\n", names(which.min(aics))))
    return(best_model)
  } else {
    return(models[[1]])
  }
}

# Enhanced phylogenetic signal checking (same as working version)
check_phylogenetic_signal <- function(data, tree, ycol, xcol = NULL) {
  cat("===== Phylogenetic Signal Analysis =====\n")
  
  if (!requireNamespace("ape", quietly = TRUE)) {
    warning("Cannot check phylogenetic signal: 'ape' package not available")
    return(list(lambda_y = NA, lambda_x = NA))
  }
  
  results <- list()
  
  # Check Y variable
  tryCatch({
    if (requireNamespace("phytools", quietly = TRUE)) {
      phylosig_y <- phytools::phylosig(tree, setNames(data[[ycol]], rownames(data)), method = "lambda")
      results$lambda_y <- phylosig_y$lambda
      cat(sprintf("Phylogenetic signal in %s (λ): %.3f (p = %.3f)\n", 
                  ycol, phylosig_y$lambda, phylosig_y$P))
    } else {
      signal_y <- ape::Moran.I(data[[ycol]], ape::cophenetic.phylo(tree))
      results$lambda_y <- signal_y$observed
      cat(sprintf("Phylogenetic signal in %s (Moran's I): %.3f (p = %.3f)\n", 
                  ycol, signal_y$observed, signal_y$p.value))
    }
  }, error = function(e) {
    warning(sprintf("Failed to calculate phylogenetic signal for %s: %s", ycol, e$message))
    results$lambda_y <- NA
  })
  
  # Check X variable if provided
  if (!is.null(xcol)) {
    tryCatch({
      if (requireNamespace("phytools", quietly = TRUE)) {
        phylosig_x <- phytools::phylosig(tree, setNames(data[[xcol]], rownames(data)), method = "lambda")
        results$lambda_x <- phylosig_x$lambda
        cat(sprintf("Phylogenetic signal in %s (λ): %.3f (p = %.3f)\n", 
                    xcol, phylosig_x$lambda, phylosig_x$P))
      } else {
        signal_x <- ape::Moran.I(data[[xcol]], ape::cophenetic.phylo(tree))
        results$lambda_x <- signal_x$observed
        cat(sprintf("Phylogenetic signal in %s (Moran's I): %.3f (p = %.3f)\n", 
                    xcol, signal_x$observed, signal_x$p.value))
      }
    }, error = function(e) {
      warning(sprintf("Failed to calculate phylogenetic signal for %s: %s", xcol, e$message))
      results$lambda_x <- NA
    })
  }
  
  # Provide recommendations
  avg_signal <- mean(c(results$lambda_y, results$lambda_x), na.rm = TRUE)
  if (!is.na(avg_signal)) {
    if (avg_signal < 0.1) {
      cat("RECOMMENDATION: Low phylogenetic signal detected. PGLS may not be necessary.\n")
    } else if (avg_signal > 0.7) {
      cat("RECOMMENDATION: Strong phylogenetic signal detected. PGLS is highly recommended.\n")
    } else {
      cat("RECOMMENDATION: Moderate phylogenetic signal. PGLS may provide some improvement.\n")
    }
  }
  
  cat("\n")
  return(results)
}

# Log axis handling functions (same as working version)
create_log_scale <- function(axis = "x", data_range, log_base, opts) {
  step <- if (is.na(opts$log_step) || opts$log_step <= 0) 1.0 else opts$log_step
  
  if (abs(log_base - 10) < 1e-8) {
    lo <- floor(log10(data_range[1]))
    hi <- ceiling(log10(data_range[2]))
    breaks <- 10^(seq(lo, hi, by = step))
  } else {
    lo <- floor(log(data_range[1], base = log_base))
    hi <- ceiling(log(data_range[2], base = log_base))
    breaks <- log_base^(seq(lo, hi, by = step))
  }
  
  labels <- create_log_labels(breaks, log_base, opts)
  
  base_name <- if (abs(log_base - exp(1)) < 1e-8) "e" else as.character(log_base)
  axis_label <- switch(tolower(opts$log_label),
    "pow10" = paste0(" (log", base_name, ")"),
    "log10" = paste0("log", base_name, "()"),
    paste0(" (log", base_name, ")")
  )
  
  return(list(breaks = breaks, labels = labels, axis_label = axis_label))
}

create_log_labels <- function(breaks, log_base, opts) {
  switch(tolower(opts$log_label),
    "pow10" = {
      if (abs(log_base - 10) < 1e-8) {
        parse(text = paste0("10^", formatC(log10(breaks), format = 'f', digits = 0)))
      } else if (abs(log_base - exp(1)) < 1e-8) {
        parse(text = paste0("e^", formatC(log(breaks), format = 'f', digits = 0)))
      } else {
        parse(text = paste0(log_base, "^", formatC(log(breaks, base = log_base), format = 'f', digits = 0)))
      }
    },
    "log10" = {
      digits <- if (is.na(opts$log_digits)) 1 else as.integer(round(opts$log_digits))
      if (abs(log_base - 10) < 1e-8) {
        formatC(log10(breaks), format = 'f', digits = digits)
      } else {
        formatC(log(breaks, base = log_base), format = 'f', digits = digits)
      }
    },
    breaks
  )
}

# Enhanced model fitting and summary (same as working version)
fit_and_summarize_models <- function(dfc, xcol, ycol, groupcol, by_group, formula_str, 
                                   use_pgls = FALSE, tree = NULL, species_col = NULL, 
                                   pgls_correlation = "brownian", check_phylosig = TRUE) {
  cat("===== Model Summary =====\n")
  
  if (use_pgls && !is.null(tree) && !is.null(species_col)) {
    # Prepare tree
    tree_fixed <- prepare_tree_for_pgls(tree, make_ultrametric = TRUE)
    
    # Prepare data for PGLS
    dfc_pgls <- prepare_pgls_data(dfc, tree_fixed, species_col)
    if (is.null(dfc_pgls)) {
      warning("Failed to prepare PGLS data. Falling back to standard linear regression.")
      use_pgls <- FALSE
    } else {
      # Check phylogenetic signal if requested
      if (check_phylosig) {
        phylosig_results <- check_phylogenetic_signal(dfc_pgls, tree_fixed, ycol, xcol)
      }
    }
  }
  
  if (use_pgls && !is.null(tree) && !is.null(species_col)) {
    # Overall PGLS model
    fit_overall <- fit_pgls_model_working(dfc_pgls, formula_str, tree_fixed, pgls_correlation)
    cat("Overall PGLS model:\n")
    
    if (inherits(fit_overall, "gls")) {
      cat("Note: PGLS corrects for phylogenetic non-independence in residuals\n")
      print(summary(fit_overall))
      
      # Extract lambda parameter if Pagel model was used
      if (inherits(fit_overall$modelStruct$corStruct, "corPagel")) {
        lambda_est <- as.numeric(coef(fit_overall$modelStruct$corStruct, unconstrained = FALSE))
        cat(sprintf("Estimated Pagel's λ: %.3f\n", lambda_est))
      }
    } else {
      cat("Note: PGLS failed, using OLS model\n")
      print(summary(fit_overall))
    }
    cat("\n")
    
  } else {
    # Standard linear regression
    fit_overall <- lm(stats::as.formula(formula_str), data = dfc)
    cat("Overall model (OLS):\n")
    print(summary(fit_overall))
    cat("\n")
  }
}

# Extract statistical information from model (same as working version)
extract_model_stats <- function(fit, xcol, ycol, method = "OLS") {
  if (inherits(fit, "gls")) {
    # PGLS model
    coeffs <- summary(fit)$tTable
    intercept <- coeffs["(Intercept)", "Value"]
    slope <- coeffs[xcol, "Value"]
    slope_pvalue <- coeffs[xcol, "p-value"]
    r_squared <- NA
    n <- nrow(fit$data)
    
    # Try to calculate pseudo R²
    tryCatch({
      y_obs <- fit$data[[ycol]]
      y_pred <- fitted(fit)
      ss_res <- sum((y_obs - y_pred)^2)
      ss_tot <- sum((y_obs - mean(y_obs))^2)
      r_squared <- 1 - (ss_res / ss_tot)
    }, error = function(e) {
      r_squared <- NA
    })
  } else {
    # OLS model
    coeffs <- summary(fit)$coefficients
    intercept <- coeffs["(Intercept)", "Estimate"]
    slope <- coeffs[xcol, "Estimate"]
    slope_pvalue <- coeffs[xcol, "Pr(>|t|)"]
    r_squared <- summary(fit)$r.squared
    n <- nrow(fit$model)
  }
  
  # Format p-value
  if (slope_pvalue < 0.001) {
    p_text <- "p < 0.001"
  } else {
    p_text <- sprintf("p = %.3f", slope_pvalue)
  }
  
  # Create significance stars
  if (slope_pvalue < 0.001) {
    sig_stars <- "***"
  } else if (slope_pvalue < 0.01) {
    sig_stars <- "**"
  } else if (slope_pvalue < 0.05) {
    sig_stars <- "*"
  } else {
    sig_stars <- ""
  }
  
  # Format equation
  if (slope >= 0) {
    equation <- sprintf("%s = %.4f + %.4f × %s", ycol, intercept, slope, xcol)
  } else {
    equation <- sprintf("%s = %.4f - %.4f × %s", ycol, intercept, abs(slope), xcol)
  }
  
  return(list(
    equation = equation,
    slope = slope,
    intercept = intercept,
    pvalue = slope_pvalue,
    p_text = p_text,
    sig_stars = sig_stars,
    r_squared = r_squared,
    n = n,
    method = method
  ))
}

# ð§ CORRECTED PLOTTING FUNCTIONS - THE KEY FIX! ð§

# NEW: Function to calculate PGLS prediction line data
calculate_pgls_predictions <- function(fit, data, xcol, ycol) {
  tryCatch({
    # Get the range of X values
    x_range <- range(data[[xcol]], na.rm = TRUE)
    x_pred <- seq(x_range[1], x_range[2], length.out = 100)
    
    # Create prediction data frame
    pred_data <- data.frame(x_pred)
    names(pred_data) <- xcol
    
    # Get PGLS predictions
    y_pred <- predict(fit, newdata = pred_data)
    
    # Create data frame for plotting
    pred_df <- data.frame(
      x = x_pred,
      y = as.numeric(y_pred),
      type = "PGLS"
    )
    names(pred_df)[1:2] <- c(xcol, ycol)
    
    return(pred_df)
  }, error = function(e) {
    warning(sprintf("Failed to calculate PGLS predictions: %s", e$message))
    return(NULL)
  })
}

# NEW: Function to calculate OLS prediction line data
calculate_ols_predictions <- function(fit, data, xcol, ycol) {
  tryCatch({
    # Get the range of X values  
    x_range <- range(data[[xcol]], na.rm = TRUE)
    x_pred <- seq(x_range[1], x_range[2], length.out = 100)
    
    # Create prediction data frame
    pred_data <- data.frame(x_pred)
    names(pred_data) <- xcol
    
    # Get OLS predictions
    y_pred <- predict(fit, newdata = pred_data)
    
    # Create data frame for plotting
    pred_df <- data.frame(
      x = x_pred,
      y = as.numeric(y_pred),
      type = "OLS"
    )
    names(pred_df)[1:2] <- c(xcol, ycol)
    
    return(pred_df)
  }, error = function(e) {
    warning(sprintf("Failed to calculate OLS predictions: %s", e$message))
    return(NULL)
  })
}

# CORRECTED plotting with proper PGLS vs OLS line distinction
plot_with_ggplot2 <- function(dfc, opts, xcol, ycol, groupcol, by_group, 
                              to_log_x, to_log_y, log_base, formula_str, 
                              use_pgls = FALSE, tree = NULL, species_col = NULL) {
  
  # Ensure clean device state
  while (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required but not available")
  }
  
  suppressPackageStartupMessages({
    library(ggplot2)
  })
  
  # Create plot based on PGLS vs OLS
  if (use_pgls && !is.null(tree) && !is.null(species_col)) {
    # Create comparison plot with both OLS and PGLS
    p <- create_pgls_comparison_plot_corrected(dfc, opts, xcol, ycol, groupcol, by_group,
                                              to_log_x, to_log_y, log_base, formula_str,
                                              tree, species_col)
  } else {
    # Standard single plot
    p <- create_single_plot(dfc, opts, xcol, ycol, groupcol, by_group,
                           to_log_x, to_log_y, log_base, formula_str)
  }
  
  # Save plot with robust device management
  tryCatch({
    if (inherits(p, "gtable") || inherits(p, "grob")) {
      # For grid objects (like arrangeGrob output)
      grDevices::pdf(file = opts$outfile, width = opts$width, height = opts$height)
      grid::grid.draw(p)
      grDevices::dev.off()
    } else {
      # For ggplot objects - use ggsave with explicit device
      ggplot2::ggsave(filename = opts$outfile, plot = p, 
                     width = opts$width, height = opts$height, 
                     dpi = opts$dpi, units = "in", device = "pdf")
    }
    
    # Final device cleanup
    while (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
    
  }, error = function(e) {
    # Emergency cleanup
    while (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
    stop(sprintf("Failed to save plot: %s", e$message))
  })
  
  return(invisible(p))
}

# Create single standard plot (same as before)
create_single_plot <- function(dfc, opts, xcol, ycol, groupcol, by_group,
                              to_log_x, to_log_y, log_base, formula_str) {
  
  # Create base plot with data
  p <- ggplot(dfc, aes(x = .data[[xcol]], y = .data[[ycol]]))
  
  # Add points
  if (!is.null(groupcol)) {
    p <- p + geom_point(aes(color = .data[[groupcol]]), shape = 16, size = 2.8, alpha = 0.9)
  } else {
    p <- p + geom_point(shape = 16, size = 2.8, alpha = 0.9)
  }
  
  # Fit model and extract statistics
  fit <- lm(stats::as.formula(formula_str), data = dfc)
  model_stats <- extract_model_stats(fit, xcol, ycol, "OLS")
  
  # Add regression lines with explicit formula to avoid warnings
  if (by_group && !is.null(groupcol)) {
    p <- p + geom_smooth(aes(color = .data[[groupcol]]), method = "lm", 
                        formula = y ~ x, se = TRUE, alpha = 0.3)
  } else {
    p <- p + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, 
                        color = "black", alpha = 0.3)
  }
  
  # Configure axes
  p <- configure_axes(p, dfc, xcol, ycol, to_log_x, to_log_y, log_base, opts)
  
  # Add statistical information to plot
  if (!is.null(model_stats)) {
    # Create statistics text
    stats_text <- paste0(
      model_stats$equation, model_stats$sig_stars, "\n",
      model_stats$p_text, "\n",
      sprintf("n = %d", model_stats$n)
    )
    
    if (!is.na(model_stats$r_squared)) {
      stats_text <- paste0(stats_text, sprintf("\nR² = %.3f", model_stats$r_squared))
    }
    
    # Add text annotation
    p <- p + annotate("text", 
                      x = -Inf, y = Inf, 
                      label = stats_text,
                      hjust = -0.1, vjust = 1.1,
                      size = 3.5, 
                      fontface = "plain",
                      color = "black")
  }
  
  # Add theme and labels
  title_text <- if (!is.na(opts$title)) opts$title else paste("Linear regression:", formula_str)
  
  p <- p +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(t = 10, r = 18, b = 14, l = 18)
    ) +
    labs(
      title = title_text,
      color = if (is.null(groupcol)) NULL else groupcol
    )
  
  # Add log ticks if needed
  if (to_log_x || to_log_y) {
    sides <- paste0(if (to_log_x) "b" else "", if (to_log_y) "l" else "")
    p <- p + annotation_logticks(sides = sides) + coord_cartesian(clip = "off")
  }
  
  return(p)
}

# ð§ CORRECTED PGLS comparison plot - THE MAIN FIX! ð§
create_pgls_comparison_plot_corrected <- function(dfc, opts, xcol, ycol, groupcol, by_group,
                                                 to_log_x, to_log_y, log_base, formula_str,
                                                 tree, species_col) {
  # Check for required packages
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    warning("gridExtra package not available. Creating standard plot only.")
    return(create_single_plot(dfc, opts, xcol, ycol, groupcol, by_group,
                             to_log_x, to_log_y, log_base, formula_str))
  }
  
  suppressPackageStartupMessages({
    library(gridExtra)
    library(grid)
  })
  
  # Prepare tree and PGLS data
  tree_fixed <- prepare_tree_for_pgls(tree, make_ultrametric = TRUE)
  dfc_pgls <- prepare_pgls_data(dfc, tree_fixed, species_col)
  if (is.null(dfc_pgls)) {
    warning("Failed to prepare PGLS data. Creating standard plot only.")
    return(create_single_plot(dfc, opts, xcol, ycol, groupcol, by_group,
                             to_log_x, to_log_y, log_base, formula_str))
  }
  
  # Fit both models on the SAME data for fair comparison
  formula_str_simple <- paste(ycol, "~", xcol)
  
  # Fit OLS model
  ols_fit <- lm(stats::as.formula(formula_str_simple), data = dfc_pgls)
  ols_stats <- extract_model_stats(ols_fit, xcol, ycol, "OLS")
  
  # Fit PGLS model  
  pgls_fit <- fit_pgls_model_working(dfc_pgls, formula_str_simple, tree_fixed, opts$pgls_correlation)
  pgls_stats <- extract_model_stats(pgls_fit, xcol, ycol, 
                                   ifelse(inherits(pgls_fit, "gls"), "PGLS", "OLS"))
  
  # Calculate prediction lines for BOTH methods
  ols_pred <- calculate_ols_predictions(ols_fit, dfc_pgls, xcol, ycol)
  pgls_pred <- if (inherits(pgls_fit, "gls")) {
    calculate_pgls_predictions(pgls_fit, dfc_pgls, xcol, ycol)
  } else {
    # If PGLS failed, use OLS predictions but label them appropriately
    pred <- calculate_ols_predictions(pgls_fit, dfc_pgls, xcol, ycol)
    if (!is.null(pred)) pred$type <- "PGLS (failed, using OLS)"
    pred
  }
  
  # Create OLS plot (left panel) with TRUE OLS line
  p1 <- create_corrected_regression_panel(dfc_pgls, xcol, ycol, groupcol, by_group,
                                         to_log_x, to_log_y, log_base, opts,
                                         "OLS (Ordinary Least Squares)", ols_stats, ols_pred)
  
  # Create PGLS plot (right panel) with TRUE PGLS line
  p2 <- create_corrected_regression_panel(dfc_pgls, xcol, ycol, groupcol, by_group,
                                         to_log_x, to_log_y, log_base, opts,
                                         paste("PGLS", paste0("(", opts$pgls_correlation, ")")), 
                                         pgls_stats, pgls_pred)
  
  # Combine plots using arrangeGrob
  main_title <- if (!is.na(opts$title)) opts$title else paste("Regression Comparison:", formula_str)
  
  combined_plot <- arrangeGrob(p1, p2, ncol = 2, 
                              top = textGrob(main_title, gp = gpar(fontsize = 14, fontface = "bold")))
  
  return(combined_plot)
}

# ð§ CORRECTED individual regression panel - SHOWS ACTUAL FITTED LINES! ð§
create_corrected_regression_panel <- function(data, xcol, ycol, groupcol, by_group,
                                             to_log_x, to_log_y, log_base, opts, 
                                             panel_title, model_stats, pred_data) {
  
  # Base plot with points
  p <- ggplot(data, aes(x = .data[[xcol]], y = .data[[ycol]]))
  
  # Add points
  if (!is.null(groupcol)) {
    p <- p + geom_point(aes(color = .data[[groupcol]]), shape = 16, size = 2.5, alpha = 0.8)
  } else {
    p <- p + geom_point(shape = 16, size = 2.5, alpha = 0.8)
  }
  
  # ð§ KEY FIX: Add the ACTUAL fitted line from the model predictions ð§
  if (!is.null(pred_data)) {
    p <- p + geom_line(data = pred_data, 
                      aes(x = .data[[xcol]], y = .data[[ycol]]), 
                      color = "red", linewidth = 1.2, alpha = 0.8, inherit.aes = FALSE)
    
    # Add a subtitle to indicate which type of line this is
    line_type <- if (grepl("PGLS", panel_title)) {
      ifelse(grepl("failed", pred_data$type[1]), "PGLS (fallback to OLS)", "PGLS fitted line")
    } else {
      "OLS fitted line"
    }
    p <- p + labs(subtitle = line_type)
  }
  
  # ð§ OPTIONAL: Also add confidence bands using geom_smooth for comparison ð§
  # This will show the difference between manual prediction and geom_smooth
  if (by_group && !is.null(groupcol)) {
    p <- p + geom_smooth(aes(color = .data[[groupcol]]), method = "lm", 
                        formula = y ~ x, se = TRUE, alpha = 0.2, linewidth = 0.5)
  } else {
    p <- p + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, 
                        color = "blue", alpha = 0.2, linewidth = 0.5)
  }
  
  # Configure axes
  p <- configure_axes(p, data, xcol, ycol, to_log_x, to_log_y, log_base, opts)
  
  # Add statistical information to plot
  if (!is.null(model_stats)) {
    # Create statistics text
    stats_text <- paste0(
      model_stats$equation, model_stats$sig_stars, "\n",
      model_stats$p_text, "\n",
      sprintf("n = %d", model_stats$n)
    )
    
    if (!is.na(model_stats$r_squared)) {
      stats_text <- paste0(stats_text, sprintf("\nR² = %.3f", model_stats$r_squared))
    }
    
    # Add method type
    stats_text <- paste0(stats_text, sprintf("\nMethod: %s", model_stats$method))
    
    # Add text annotation
    p <- p + annotate("text", 
                      x = -Inf, y = Inf, 
                      label = stats_text,
                      hjust = -0.1, vjust = 1.1,
                      size = 3, 
                      fontface = "plain",
                      color = "black")
  }
  
  # Add theme and labels
  p <- p +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "red"),
      plot.margin = margin(t = 5, r = 10, b = 10, l = 10),
      legend.position = if (!is.null(groupcol)) "bottom" else "none",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    ) +
    labs(
      title = panel_title,
      color = if (is.null(groupcol)) NULL else groupcol
    )
  
  # Add log ticks if needed
  if (to_log_x || to_log_y) {
    sides <- paste0(if (to_log_x) "b" else "", if (to_log_y) "l" else "")
    p <- p + annotation_logticks(sides = sides) + coord_cartesian(clip = "off")
  }
  
  return(p)
}

configure_axes <- function(p, dfc, xcol, ycol, to_log_x, to_log_y, log_base, opts) {
  # X-axis configuration
  if (to_log_x) {
    x_range <- range(dfc[[xcol]], na.rm = TRUE)
    x_scale <- create_log_scale("x", x_range, log_base, opts)
    
    if (abs(log_base - 10) < 1e-8) {
      p <- p + scale_x_log10(breaks = x_scale$breaks, labels = x_scale$labels,
                            expand = expansion(mult = c(0.03, 0.06)))
    } else if (requireNamespace("scales", quietly = TRUE)) {
      p <- p + scale_x_continuous(trans = scales::log_trans(base = log_base),
                                 breaks = x_scale$breaks, labels = x_scale$labels,
                                 expand = expansion(mult = c(0.03, 0.06)))
    } else {
      p <- p + scale_x_log10(breaks = x_scale$breaks, labels = x_scale$labels,
                            expand = expansion(mult = c(0.03, 0.06)))
    }
    
    x_label <- paste0(xcol, x_scale$axis_label)
  } else {
    x_label <- xcol
  }
  
  # Y-axis configuration
  if (to_log_y) {
    y_range <- range(dfc[[ycol]], na.rm = TRUE)
    y_scale <- create_log_scale("y", y_range, log_base, opts)
    
    if (abs(log_base - 10) < 1e-8) {
      p <- p + scale_y_log10(breaks = y_scale$breaks, labels = y_scale$labels,
                            expand = expansion(mult = c(0.03, 0.06)))
    } else if (requireNamespace("scales", quietly = TRUE)) {
      p <- p + scale_y_continuous(trans = scales::log_trans(base = log_base),
                                 breaks = y_scale$breaks, labels = y_scale$labels,
                                 expand = expansion(mult = c(0.03, 0.06)))
    } else {
      p <- p + scale_y_log10(breaks = y_scale$breaks, labels = y_scale$labels,
                            expand = expansion(mult = c(0.03, 0.06)))
    }
    
    y_label <- paste0(ycol, y_scale$axis_label)
  } else {
    y_label <- ycol
  }
  
  p <- p + labs(x = x_label, y = y_label)
  return(p)
}

# Help function
print_usage <- function() {
  cat("
Corrected PGLS Plotting Script - Shows ACTUAL Different Fitted Lines

Usage: Rscript plot_lm_pgls_corrected.R [OPTIONS]

Basic Options:
  --input PATH           CSV input file (default: test2.csv)
  --xcol NAME            X-axis column name (default: Size_Value_cm)
  --ycol NAME            Y-axis column name (default: genomesize)
  --groupcol NAME        Grouping column for colors (default: group)
  --by-group [BOOL]      Fit separate lines per group (default: false)
  --outfile PATH         Output PDF path (auto-generated if omitted)
  --title TEXT           Plot title

PGLS Options:
  --pgls [BOOL]          Use PGLS correction (default: false)
  --pgls-correlation TYPE Correlation structure: brownian|pagel|grafen|auto (default: brownian)
  --tree-file PATH       Phylogenetic tree file (Newick format)
  --species-col NAME     Species column name (default: species)
  --check-phylosig [BOOL] Check phylogenetic signal (default: true)
  --make-ultrametric [BOOL] Automatically make tree ultrametric (default: true)

Output Options:
  --width NUM            Figure width in inches (default: 7)
  --height NUM           Figure height in inches (default: 5.5)
  --dpi NUM              DPI for embedded rasters (default: 300)

Example:
  Rscript plot_lm.R --input genome_Phenotype.csv --xcol DNA --ycol genomesize \\
                                   --pgls true --tree-file ref.tree --species-col Species \\
                                   --width 15 --outfile corrected_output.pdf --pgls-correlation auto

")
}

# Enhanced main function (same structure as working version)
main <- function() {
  # Comprehensive device cleanup
  while (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }
  
  args <- commandArgs(trailingOnly = TRUE)
  opts <- parse_args(args)
  
  if (isTRUE(opts$help)) {
    print_usage()
    quit(status = 0)
  }
  
  # Validate input file
  if (!file.exists(opts$input)) {
    stop(sprintf("Input file not found: %s", opts$input))
  }
  
  cat("===== Data Loading and Validation =====\n")
  
  # Read and validate data
  df <- read.csv(opts$input, stringsAsFactors = FALSE, check.names = FALSE)
  cat(sprintf("Loaded data: %d rows, %d columns\n", nrow(df), ncol(df)))
  
  opts$groupcol <- validate_columns(df, opts$xcol, opts$ycol, opts$groupcol)
  
  # Load phylogenetic tree if PGLS is requested
  tree <- NULL
  if (isTRUE(opts$pgls)) {
    if (is.na(opts$tree_file) || !file.exists(opts$tree_file)) {
      stop("PGLS requires a valid tree file. Use --tree-file option.")
    }
    if (!requireNamespace("ape", quietly = TRUE)) {
      stop("PGLS requires the 'ape' package. Please install it: install.packages('ape')")
    }
    
    cat("Loading phylogenetic tree...\n")
    tree <- ape::read.tree(opts$tree_file)
    cat(sprintf("Loaded phylogenetic tree with %d tips\n", length(tree$tip.label)))
    
    # Validate species column
    if (!(opts$species_col %in% names(df))) {
      stop(sprintf("Species column '%s' not found in data", opts$species_col))
    }
  }
  
  # Clean and prepare data
  cat("Cleaning and preparing data...\n")
  dfc <- clean_numeric_data(df, opts$xcol, opts$ycol)
  validate_data_quality(dfc, opts$xcol, opts$ycol)
  
  # Determine log transformations
  to_log_x <- tolower(opts$log) %in% c("x", "xy", "both")
  to_log_y <- tolower(opts$log) %in% c("y", "xy", "both")
  log_base <- parse_log_base(opts$log_base)
  
  # Filter data for log transformations
  if (to_log_x || to_log_y) {
    dfc <- filter_positive_values(dfc, opts$xcol, opts$ycol, to_log_x, to_log_y)
  }
  
  # Create formula string for display
  logfun_name <- if (abs(log_base - exp(1)) < 1e-8) "log" else paste0("log", log_base)
  lhs <- if (to_log_y) paste0(logfun_name, "(", opts$ycol, ")") else opts$ycol
  rhs <- if (to_log_x) paste0(logfun_name, "(", opts$xcol, ")") else opts$xcol
  formula_str <- paste(lhs, "~", rhs)
  
  cat("===== Creating Plot =====\n")
  
  # Create plot
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_with_ggplot2(dfc, opts, opts$xcol, opts$ycol, opts$groupcol, 
                      isTRUE(opts$by_group), to_log_x, to_log_y, log_base, formula_str,
                      isTRUE(opts$pgls), tree, opts$species_col)
  } else {
    stop("ggplot2 package is required. Please install it: install.packages('ggplot2')")
  }
  
  cat("===== Statistical Analysis =====\n")
  
  # Print model summaries
  fit_and_summarize_models(dfc, opts$xcol, opts$ycol, opts$groupcol, 
                          isTRUE(opts$by_group), formula_str, 
                          isTRUE(opts$pgls), tree, opts$species_col,
                          opts$pgls_correlation, isTRUE(opts$check_phylosig))
  
  cat(sprintf("===== Output =====\n"))
  cat(sprintf("Plot saved to: %s\n", opts$outfile))
  
  # Final cleanup
  while (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }
  
  # Remove any accidentally created Rplots.pdf
  if (file.exists("Rplots.pdf")) {
    file.remove("Rplots.pdf")
    cat("Removed unwanted Rplots.pdf\n")
  }
  
  cat("ð Analysis completed with CORRECTED PGLS vs OLS visualization! ð\n")
  return(invisible(TRUE))
}

# Enhanced error handling wrapper
tryCatch({
  # Ensure clean environment
  while (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }
  
  # Run main function
  main()
  
}, error = function(e) {
  # Emergency cleanup
  while (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }
  
  # Remove any accidentally created files
  if (file.exists("Rplots.pdf")) {
    file.remove("Rplots.pdf")
  }
  
  cat(sprintf("Error: %s\n", e$message))
  cat("\nFor help, run: Rscript plot_lm_pgls_corrected.R --help\n")
  quit(status = 1)
  
}, finally = {
  # Final cleanup
  while (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }
})
