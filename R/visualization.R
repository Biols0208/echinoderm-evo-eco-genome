# ============================================================
# Visualization Module for PGLS Analysis Framework
# Creates forest plots and diagnostic plots
# ============================================================

source("R/utils.R")

# Create enhanced forest plot with configuration options
create_forest_plot <- function(results_df, config, title = NULL) {
  if (is.null(results_df) || nrow(results_df) == 0) {
    if (config$output$verbose) cat("No results for forest plot\n")
    return(NULL)
  }
  
  # Filter valid results
  valid_results <- results_df %>% filter(!is.na(estimate))
  
  if (nrow(valid_results) == 0) {
    if (config$output$verbose) cat("No valid estimates for forest plot\n")
    return(NULL)
  }
  
  # Create display labels for multi-level predictors
  if ("term" %in% names(valid_results)) {
    valid_results <- valid_results %>%
      mutate(
        # Extract the predictor level from term
        predictor_level = gsub(paste0("^", predictor), "", term),
        # Create display label
        model_display = ifelse(
          !is.na(term) & term != "" & !grepl("Intercept", term),
          paste0(model, " [", predictor_level, "]"),
          model
        )
      )
  } else {
    valid_results$model_display <- valid_results$model
  }
  
  # Set model order to avoid overlaps
  valid_results$model_display <- factor(valid_results$model_display, 
                                        levels = rev(unique(valid_results$model_display)))
  
  # Set theme based on configuration
  if (config$visualization$use_custom_theme) {
    theme_forest <- theme_bw(base_size = config$visualization$base_size)
  } else {
    theme_forest <- theme_minimal(base_size = config$visualization$base_size)
  }
  
  # Apply grid line configuration
  theme_forest <- theme_forest +
    theme(
      axis.title.y = element_blank(),
      legend.position = "none",
      text = element_text(family = config$visualization$font_family)
    )
  
  if (config$visualization$grid_lines == "none") {
    theme_forest <- theme_forest + theme(panel.grid = element_blank())
  } else if (config$visualization$grid_lines == "minimal") {
    theme_forest <- theme_forest + theme(panel.grid.minor = element_blank())
  }
  
  # Set color based on scheme
  point_color <- switch(config$visualization$forest_color_scheme,
    "viridis" = "#440154FF",
    "grayscale" = "gray30",
    "darkblue"  # default
  )
  
  # Main forest plot
  p_main <- ggplot(valid_results, aes(x = estimate, y = model_display)) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                   height = config$visualization$forest_ci_thickness) +
    geom_point(size = config$visualization$forest_point_size, color = point_color) +
    labs(
      x = "Effect Size",
      title = ifelse(is.null(title), "Forest Plot", title),
      y = NULL
    ) +
    theme_forest +
    theme(
      plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
      axis.text.y = element_text(size = config$visualization$base_size * 0.9)
    )
  
  # Add percentage change labels if configured and applicable
  if (config$visualization$show_percentage_change && "pct_change" %in% names(valid_results) && 
      any(!is.na(valid_results$pct_change))) {
    p_main <- p_main +
      geom_text(
        aes(x = conf.high,
            label = sprintf("  %.1f%% [%.1f, %.1f]",
                          pct_change, pct_change_low, pct_change_high)),
        hjust = 0, vjust = 0.5, 
        size = config$visualization$base_size * 0.25, 
        family = config$visualization$font_family,
        data = filter(valid_results, !is.na(pct_change))
      ) +
      coord_cartesian(clip = "off") +
      theme(plot.margin = margin(5.5, 80, 5.5, 5.5, "pt"))
  }
  
  # Create summary table
  table_df <- valid_results %>%
    mutate(
      lambda_fmt = case_when(
        grepl("GLS\\(IID\\)", model) ~ "-",
        grepl("lambda=0", model) ~ "0",
        grepl("lambda=1", model) ~ "1",
        TRUE ~ ifelse(is.na(lambda), "-", sprintf("%.3f", lambda))
      )
    ) %>%
    select(model_display, stars, n, lambda_fmt) %>%
    rename(Model = model_display, Signif = stars, N = n, Lambda = lambda_fmt)
  
  # Create table as separate plot
  p_table <- ggplot(table_df, aes(y = Model)) +
    geom_blank() +
    annotate("text", x = 1, y = nlevels(table_df$Model), 
             label = "Signif", vjust = -1.5, fontface = "bold", size = 3.5) +
    annotate("text", x = 2, y = nlevels(table_df$Model), 
             label = "N", vjust = -1.5, fontface = "bold", size = 3.5) +
    annotate("text", x = 3, y = nlevels(table_df$Model), 
             label = "Lambda", vjust = -1.5, fontface = "bold", size = 3.5) +
    geom_text(aes(x = 1, label = Signif), size = 3.5) +
    geom_text(aes(x = 2, label = N), size = 3.5) +
    geom_text(aes(x = 3, label = Lambda), size = 3.5) +
    scale_x_continuous(limits = c(0.5, 3.5), breaks = NULL) +
    scale_y_discrete(limits = levels(table_df$Model)) +
    theme_void() +
    theme(plot.margin = margin(28, 5.5, 5.5, 0, "pt"))
  
  # Combine plots
  combined <- p_main + p_table + plot_layout(widths = c(3, 1))
  
  return(combined)
}

# Create diagnostic plots for model assessment
create_diagnostic_plots <- function(model, model_name, config) {
  if (is.null(model)) return(NULL)
  
  tryCatch({
    # Extract residuals and fitted values
    res <- residuals(model, type = "normalized")
    fit <- fitted(model)
    
    plots <- list()
    
    # Residuals vs Fitted
    plots$resid_fitted <- ggplot(data.frame(fit = fit, res = res), aes(x = fit, y = res)) +
      geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "loess", se = FALSE, color = "red", size = 0.8) +
      labs(
        x = "Fitted values",
        y = "Normalized residuals",
        title = paste(model_name, "- Residuals vs Fitted")
      ) +
      theme_minimal(base_size = config$visualization$base_size)
    
    # Q-Q plot
    qq <- qqnorm(res, plot.it = FALSE)
    plots$qq <- ggplot(data.frame(theoretical = qq$x, sample = qq$y),
                 aes(x = theoretical, y = sample)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, alpha = 0.5) +
      geom_point(alpha = 0.6) +
      labs(
        x = "Theoretical Quantiles",
        y = "Sample Quantiles",
        title = paste(model_name, "- Normal Q-Q")
      ) +
      theme_minimal(base_size = config$visualization$base_size)
    
    # Combine plots
    return(plots$resid_fitted + plots$qq)
    
  }, error = function(e) {
    if (config$output$debug_mode) {
      warning(sprintf("Could not create diagnostics for %s: %s", model_name, e$message))
    }
    NULL
  })
}