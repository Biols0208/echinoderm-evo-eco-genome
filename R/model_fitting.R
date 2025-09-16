# ============================================================
# Model Fitting Module for PGLS Analysis Framework
# Handles PGLS/GLS model fitting and pairwise comparisons
# ============================================================

source("R/utils.R")

# Prepare tree and data for analysis
prep_tree_data <- function(phy, df, species_col, min_group_size = 3) {
  if (is.null(phy) || is.null(df)) return(NULL)
  
  # Find common species
  common_species <- intersect(phy$tip.label, df[[species_col]])
  
  if (length(common_species) < min_group_size) {
    warning(sprintf("Too few common species (%d < %d)", 
                   length(common_species), min_group_size))
    return(NULL)
  }
  
  # Subset tree and data
  phy_use <- drop.tip(phy, setdiff(phy$tip.label, common_species))
  df_use <- df[df[[species_col]] %in% common_species, ]
  
  # Order data to match tree
  df_use <- df_use[match(phy_use$tip.label, df_use[[species_col]]), ]
  
  return(list(tree = phy_use, data = df_use))
}

# Fit PGLS/GLS models with different correlation structures
fit_model <- function(df, phy, response_var, predictor_var, model_type, config) {
  # Create formula
  form <- as.formula(paste(response_var, "~", predictor_var))
  
  # Set correlation structure based on model type
  if (model_type %in% c("pgls_ml", "pgls_lambda1", "pgls_lambda0")) {
    if (config$models$correlation_structure == "brownian") {
      cor_struct <- corBrownian(phy = phy, 
                               form = as.formula(paste("~", config$variables$species_col)))
    } else {
      # Pagel correlation
      lambda_val <- switch(model_type,
        pgls_ml = 0.5,
        pgls_lambda1 = 1,
        pgls_lambda0 = 0
      )
      fixed_lambda <- model_type != "pgls_ml"
      
      cor_struct <- corPagel(value = lambda_val, phy = phy, 
                            fixed = fixed_lambda, 
                            form = as.formula(paste("~", config$variables$species_col)))
    }
  } else {
    cor_struct <- NULL
  }
  
  # Set variance structure for heteroscedasticity
  wts <- NULL
  if (config$models$allow_heteroscedasticity && length(unique(df[[predictor_var]])) > 1) {
    wts <- varIdent(form = as.formula(paste("~ 1 |", predictor_var)))
  }
  
  # Fit model with specified optimization method
  model <- tryCatch({
    gls(form, data = df, method = config$models$optimization_method,
        correlation = cor_struct, weights = wts, na.action = na.omit)
  }, error = function(e) {
    if (config$models$allow_heteroscedasticity && !is.null(wts)) {
      # Try without heteroscedasticity
      if (config$output$debug_mode) {
        cat(sprintf("Retrying without heteroscedasticity: %s\n", e$message))
      }
      gls(form, data = df, method = config$models$optimization_method,
          correlation = cor_struct, weights = NULL, na.action = na.omit)
    } else {
      if (config$output$stop_on_error) {
        stop(e)
      } else {
        warning(e$message)
        NULL
      }
    }
  })
  
  return(model)
}

# Extract comprehensive model information with enhanced completeness checking
extract_model_info <- function(mod, df, label, response_var, predictor_var, config,
                              predictor_term = NULL) {
  
  if (is.null(mod)) {
    return(tibble(
      model = label,
      response = response_var,
      predictor = predictor_var,
      term = predictor_term,
      estimate = NA, se = NA, conf.low = NA, conf.high = NA,
      p.value = NA, stars = "",
      n = nrow(df), lambda = NA,
      AIC = NA, AICc = NA, BIC = NA, R2 = NA, logLik = NA
    ))
  }
  
  # Get coefficient table
  coef_table <- summary(mod)$tTable
  
  # Find relevant terms (exclude intercept)
  if (is.null(predictor_term)) {
    terms <- rownames(coef_table)[grepl(predictor_var, rownames(coef_table))]
  } else {
    terms <- predictor_term
  }
  
  # For categorical predictors, add reference level information if verbose
  if (config$output$debug_mode && length(unique(df[[predictor_var]])) > 2) {
    # Get all levels of the predictor
    all_levels <- levels(factor(df[[predictor_var]]))
    
    # Find reference level (the one not in coefficient terms)
    coef_levels <- gsub(paste0("^", predictor_var), "", terms)
    ref_level <- setdiff(all_levels, coef_levels)[1]
    
    if (!is.null(ref_level) && config$output$verbose) {
      cat(sprintf("  Reference level for %s: %s (coefficient = 0)\n", 
                 predictor_var, ref_level))
    }
  }
  
  # Extract info for each term
  results <- map_dfr(terms, function(term) {
    if (!term %in% rownames(coef_table)) {
      return(tibble(
        model = label, response = response_var, predictor = predictor_var,
        term = term, estimate = NA, se = NA, conf.low = NA, conf.high = NA,
        p.value = NA, stars = "", n = nrow(df), lambda = NA,
        AIC = NA, AICc = NA, BIC = NA, R2 = NA, logLik = NA
      ))
    }
    
    # Extract coefficients
    est <- coef_table[term, "Value"]
    se <- coef_table[term, "Std.Error"]
    pval <- coef_table[term, "p-value"]
    
    # Calculate confidence intervals
    t_crit <- qt(0.5 + config$statistics$confidence_level/2, 
                 df = nrow(df) - length(coef(mod)))
    conf.low <- est - t_crit * se
    conf.high <- est + t_crit * se
    
    # Try to get better CIs from intervals()
    tryCatch({
      ints <- intervals(mod, level = config$statistics$confidence_level)$coef
      if (term %in% rownames(ints)) {
        conf.low <- ints[term, "lower"]
        conf.high <- ints[term, "upper"]
      }
    }, error = function(e) {
      # Keep the normal approximation CIs if intervals() fails
      NULL
    })
    
    # Add significance stars
    stars <- add_significance_stars(pval, config$statistics$alpha_levels, 
                                   config$statistics$signif_symbols)
    
    # Extract lambda
    lambda <- tryCatch({
      coef(mod$modelStruct$corStruct, unconstrained = FALSE)
    }, error = function(e) NA_real_)
    
    # Calculate R-squared if requested
    r2 <- if (config$models$calculate_r2) {
      tryCatch({
        1 - var(residuals(mod)) / var(df[[response_var]], na.rm = TRUE)
      }, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    
    tibble(
      model = label,
      response = response_var,
      predictor = predictor_var,
      term = term,
      estimate = est,
      se = se,
      conf.low = conf.low,
      conf.high = conf.high,
      p.value = pval,
      stars = stars,
      n = nrow(df),
      lambda = lambda,
      AIC = AIC(mod),
      AICc = MuMIn::AICc(mod),
      BIC = BIC(mod),
      R2 = r2,
      logLik = as.numeric(logLik(mod))
    )
  })
  
  return(results)
}

# Perform pairwise comparisons for multi-level factors
perform_pairwise_comparisons <- function(mod, predictor_var, response_var, 
                                        label_prefix, config) {
  if (!config$statistics$perform_pairwise || is.null(mod)) return(NULL)
  
  # Get model data and check levels
  model_data <- mod$data
  n_levels <- length(unique(model_data[[predictor_var]]))
  
  if (n_levels <= 2) return(NULL)  # No need for pairwise if only 2 levels
  
  if (config$output$verbose) {
    cat(sprintf("  Performing pairwise comparisons for %d levels...\n", n_levels))
  }
  
  # Get unique levels
  levels_list <- sort(unique(as.character(model_data[[predictor_var]])))
  
  # Initialize results
  pairwise_results <- list()
  
  # Method 1: Using emmeans (most robust)
  if (config$statistics$pairwise_method == "emmeans" && 
      requireNamespace("emmeans", quietly = TRUE)) {
    tryCatch({
      # Create emmeans object
      emm <- emmeans(mod, specs = predictor_var)
      
      # Get all pairwise comparisons
      pairs_result <- pairs(emm, adjust = ifelse(config$statistics$adjust_p_values, 
                                               config$statistics$p_adjust_method, "none"))
      
      # Convert to tibble
      pairs_df <- as.data.frame(pairs_result)
      
      # Process each comparison
      for (i in 1:nrow(pairs_df)) {
        comparison_name <- as.character(pairs_df$contrast[i])
        
        pairwise_results[[i]] <- tibble(
          model = paste0(label_prefix, " [", comparison_name, "]"),
          response = response_var,
          predictor = predictor_var,
          term = paste0(predictor_var, "_pairwise"),
          comparison = comparison_name,
          estimate = pairs_df$estimate[i],
          se = pairs_df$SE[i],
          conf.low = pairs_df$estimate[i] - 1.96 * pairs_df$SE[i],
          conf.high = pairs_df$estimate[i] + 1.96 * pairs_df$SE[i],
          p.value = pairs_df$p.value[i],
          p.adjusted = if (config$statistics$adjust_p_values) pairs_df$p.value[i] else NA_real_,
          stars = add_significance_stars(pairs_df$p.value[i]),
          n = nrow(model_data),
          lambda = NA_real_,
          AIC = NA_real_,
          AICc = NA_real_,
          BIC = NA_real_,
          R2 = NA_real_,
          logLik = NA_real_
        )
      }
      
    }, error = function(e) {
      if (config$output$debug_mode) warning("emmeans method failed: ", e$message)
      config$statistics$pairwise_method <<- "manual"
    })
  }
  
  # Manual calculation (fallback method)
  if (config$statistics$pairwise_method == "manual" || length(pairwise_results) == 0) {
    # Get fitted values by group
    fitted_by_group <- model_data %>%
      mutate(fitted_vals = fitted(mod)) %>%
      group_by(.data[[predictor_var]]) %>%
      summarise(
        mean_fitted = mean(fitted_vals, na.rm = TRUE),
        se_fitted = sd(fitted_vals, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )
    
    # Create all pairwise combinations
    idx <- 1
    for (i in 1:(length(levels_list)-1)) {
      for (j in (i+1):length(levels_list)) {
        level1 <- levels_list[i]
        level2 <- levels_list[j]
        comparison_name <- paste(level2, "-", level1)
        
        # Get means and SEs
        mean1 <- fitted_by_group$mean_fitted[fitted_by_group[[predictor_var]] == level1]
        mean2 <- fitted_by_group$mean_fitted[fitted_by_group[[predictor_var]] == level2]
        se1 <- fitted_by_group$se_fitted[fitted_by_group[[predictor_var]] == level1]
        se2 <- fitted_by_group$se_fitted[fitted_by_group[[predictor_var]] == level2]
        
        # Calculate difference and pooled SE
        diff <- mean2 - mean1
        se_diff <- sqrt(se1^2 + se2^2)
        
        # Calculate t-statistic and p-value (approximate)
        t_stat <- diff / se_diff
        df_approx <- nrow(model_data) - length(levels_list)
        p_val <- 2 * pt(abs(t_stat), df = df_approx, lower.tail = FALSE)
        
        # Apply multiple comparison correction if needed
        if (config$statistics$adjust_p_values) {
          n_comparisons <- choose(length(levels_list), 2)
          if (config$statistics$p_adjust_method == "bonferroni") {
            p_val <- min(1, p_val * n_comparisons)
          }
        }
        
        pairwise_results[[idx]] <- tibble(
          model = paste0(label_prefix, " [", comparison_name, "]"),
          response = response_var,
          predictor = predictor_var,
          term = paste0(predictor_var, "_pairwise"),
          comparison = comparison_name,
          estimate = diff,
          se = se_diff,
          conf.low = diff - 1.96 * se_diff,
          conf.high = diff + 1.96 * se_diff,
          p.value = p_val,
          p.adjusted = if (config$statistics$adjust_p_values) p_val else NA_real_,
          stars = add_significance_stars(p_val),
          n = nrow(model_data),
          lambda = NA_real_,
          AIC = NA_real_,
          AICc = NA_real_,
          BIC = NA_real_,
          R2 = NA_real_,
          logLik = NA_real_
        )
        
        idx <- idx + 1
      }
    }
  }
  
  # Combine results
  if (length(pairwise_results) > 0) {
    final_results <- bind_rows(pairwise_results)
    
    if (config$output$verbose) {
      cat(sprintf("  Computed %d pairwise comparisons\n", nrow(final_results)))
    }
    
    return(final_results)
  }
  
  return(NULL)
}