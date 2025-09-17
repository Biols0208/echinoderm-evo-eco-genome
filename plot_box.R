#!/usr/bin/env Rscript

# Enhanced box/violin plot script with advanced statistical options and visualization features
# Author: Enhanced version with multiple statistical tests and visualization options

# Enhanced default options
DEFAULT_OPTS <- list(
  input = "test.csv",
  xcol = "Reproductive_Mode",
  ycol = "genomesize", 
  groupcol = "group",
  type = "box",
  log = FALSE,
  log_base = "10",
  log_label = "pow10",
  log_step = 1.0,
  log_digits = 1,
  outfile = NA_character_,
  width = 8,
  height = 5,
  dpi = 300,
  title = NA_character_,
  rotate_x = FALSE,
  jitter_width = 0.15,
  point_size = 2,
  point_alpha = 0.8,
  help = FALSE,
  
  # New enhanced options
  show_mean = FALSE,
  show_median = FALSE,
  show_n = FALSE,
  show_outliers = TRUE,
  show_ci = FALSE,
  ci_level = 0.95,
  color_palette = "viridis",
  theme_style = "bw",
  stat_test = "auto",
  pairwise = FALSE,
  effect_size = FALSE,
  normality_test = FALSE,
  export_stats = FALSE,
  facet_col = NA_character_,
  interactive = FALSE,
  violin_trim = TRUE,
  box_notch = FALSE,
  density_plot = FALSE,
  trend_line = FALSE,
  verbose = FALSE
)

# Load required libraries
load_libraries <- function() {
  required_packages <- c("ggplot2", "dplyr", "scales")
  optional_packages <- c("plotly", "viridis", "RColorBrewer", "car", "multcomp", "htmlwidgets")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Required package '%s' not found. Please install it.", pkg))
    }
  }
  
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(scales)
  })
  
  # Load optional packages
  for (pkg in optional_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
}

# Enhanced statistical testing
perform_comprehensive_tests <- function(df, xcol, ycol, opts) {
  cat("===== Comprehensive Statistical Analysis =====\n")
  
  x_groups <- unique(df[[xcol]])
  n_groups <- length(x_groups)
  results <- list()
  
  if (n_groups < 2) {
    cat("Only one group found, no statistical tests performed.\n\n")
    return(invisible(results))
  }
  
  # Normality tests
  if (opts$normality_test) {
    cat("\n--- Normality Tests ---\n")
    for (group in x_groups) {
      group_data <- df[df[[xcol]] == group, ycol]
      if (length(group_data) >= 3) {
        tryCatch({
          sw_test <- shapiro.test(group_data)
          cat(sprintf("%s: Shapiro-Wilk W = %.4f, p = %.3e\n", 
                     group, sw_test$statistic, sw_test$p.value))
          results[[paste0("normality_", group)]] <- sw_test
        }, error = function(e) {
          cat(sprintf("%s: Normality test failed\n", group))
        })
      }
    }
  }
  
  # Homogeneity of variance test
  if (requireNamespace("car", quietly = TRUE) && n_groups > 1) {
    tryCatch({
      levene_test <- car::leveneTest(df[[ycol]], df[[xcol]])
      cat(sprintf("\nLevene's test for homogeneity: F = %.3f, p = %.3e\n", 
                 levene_test$`F value`[1], levene_test$`Pr(>F)`[1]))
      results$levene <- levene_test
    }, error = function(e) {
      cat("Levene's test failed:", e$message, "\n")
    })
  }
  
  # Main statistical tests
  cat("\n--- Group Comparison Tests ---\n")
  
  # Kruskal-Wallis test (non-parametric)
  tryCatch({
    kw_test <- kruskal.test(df[[ycol]], df[[xcol]])
    cat(sprintf("Kruskal-Wallis: H = %.3f, p = %.3e\n", 
               kw_test$statistic, kw_test$p.value))
    results$kruskal_wallis <- kw_test
  }, error = function(e) {
    cat("Kruskal-Wallis test failed:", e$message, "\n")
  })
  
  # ANOVA (parametric)
  tryCatch({
    aov_test <- aov(df[[ycol]] ~ df[[xcol]])
    aov_summary <- summary(aov_test)
    f_stat <- aov_summary[[1]][["F value"]][1]
    p_val <- aov_summary[[1]][["Pr(>F)"]][1]
    cat(sprintf("ANOVA: F = %.3f, p = %.3e\n", f_stat, p_val))
    results$anova <- aov_summary
    
    # Effect size (eta-squared)
    if (opts$effect_size) {
      ss_total <- sum(aov_summary[[1]][["Sum Sq"]])
      ss_group <- aov_summary[[1]][["Sum Sq"]][1]
      eta_squared <- ss_group / ss_total
      cat(sprintf("Effect size (η²): %.3f\n", eta_squared))
      results$eta_squared <- eta_squared
    }
  }, error = function(e) {
    cat("ANOVA test failed:", e$message, "\n")
  })
  
  # Welch's ANOVA (for unequal variances)
  tryCatch({
    welch_test <- oneway.test(df[[ycol]] ~ df[[xcol]])
    cat(sprintf("Welch's ANOVA: F = %.3f, p = %.3e\n", 
               welch_test$statistic, welch_test$p.value))
    results$welch <- welch_test
  }, error = function(e) {
    cat("Welch's ANOVA failed:", e$message, "\n")
  })
  
  # Pairwise comparisons
  if (opts$pairwise && n_groups > 2) {
    cat("\n--- Pairwise Comparisons ---\n")
    
    # Pairwise t-tests with Bonferroni correction
    tryCatch({
      pairwise_t <- pairwise.t.test(df[[ycol]], df[[xcol]], 
                                   p.adjust.method = "bonferroni")
      cat("Pairwise t-tests (Bonferroni corrected):\n")
      print(pairwise_t$p.value)
      results$pairwise_t <- pairwise_t
    }, error = function(e) {
      cat("Pairwise t-tests failed:", e$message, "\n")
    })
    
    # Dunn's test for multiple comparisons (non-parametric)
    tryCatch({
      pairwise_wilcox <- pairwise.wilcox.test(df[[ycol]], df[[xcol]], 
                                             p.adjust.method = "bonferroni")
      cat("\nPairwise Wilcoxon tests (Bonferroni corrected):\n")
      print(pairwise_wilcox$p.value)
      results$pairwise_wilcox <- pairwise_wilcox
    }, error = function(e) {
      cat("Pairwise Wilcoxon tests failed:", e$message, "\n")
    })
    
    # Tukey HSD (if ANOVA is significant)
    if (requireNamespace("multcomp", quietly = TRUE)) {
      tryCatch({
        aov_fit <- aov(df[[ycol]] ~ df[[xcol]])
        tukey_test <- multcomp::glht(aov_fit, linfct = multcomp::mcp(`df[[xcol]]` = "Tukey"))
        cat("\nTukey HSD results:\n")
        print(summary(tukey_test))
        results$tukey <- tukey_test
      }, error = function(e) {
        cat("Tukey HSD failed:", e$message, "\n")
      })
    }
  }
  
  # Export statistical results
  if (opts$export_stats) {
    stats_file <- sub("\\.[^.]*$", "_statistics.txt", opts$outfile)
    capture.output(
      {
        cat("Statistical Analysis Results\n")
        cat("============================\n\n")
        for (name in names(results)) {
          cat(sprintf("\n%s:\n", name))
          print(results[[name]])
        }
      },
      file = stats_file
    )
    cat(sprintf("\nStatistical results exported to: %s\n", stats_file))
  }
  
  cat("\n")
  return(invisible(results))
}

# Enhanced data summary with more statistics
enhanced_data_summary <- function(df, xcol, ycol, groupcol, opts) {
  cat("===== Enhanced Data Summary =====\n")
  cat(sprintf("Total observations: %d\n", nrow(df)))
  cat(sprintf("X variable (%s): %d unique values\n", xcol, length(unique(df[[xcol]]))))
  
  # Y variable statistics
  y_stats <- summary(df[[ycol]])
  cat(sprintf("Y variable (%s):\n", ycol))
  cat(sprintf("  Range: [%.3e, %.3e]\n", min(df[[ycol]], na.rm = TRUE), max(df[[ycol]], na.rm = TRUE)))
  cat(sprintf("  Mean ± SD: %.3e ± %.3e\n", mean(df[[ycol]], na.rm = TRUE), sd(df[[ycol]], na.rm = TRUE)))
  cat(sprintf("  Median [IQR]: %.3e [%.3e, %.3e]\n", 
             median(df[[ycol]], na.rm = TRUE), 
             quantile(df[[ycol]], 0.25, na.rm = TRUE),
             quantile(df[[ycol]], 0.75, na.rm = TRUE)))
  
  if (!is.null(groupcol)) {
    groups <- unique(df[[groupcol]])
    cat(sprintf("Groups (%s): %d levels - %s\n", groupcol, length(groups), paste(groups, collapse = ", "))) 
  }
  
  # Detailed summary by categories
  cat("\nDetailed summary by categories:\n")
  summary_stats <- aggregate(df[[ycol]],
            by = list(xcol = df[[xcol]]),
            FUN = function(x) {
					c(n = length(x),
						mean = mean(x, na.rm = TRUE),
						median = median(x, na.rm = TRUE),
						sd = sd(x, na.rm = TRUE),
                        min = min(x, na.rm = TRUE),
                        max = max(x, na.rm = TRUE),
                        q25 = quantile(x, 0.25, na.rm = TRUE),
                        q75 = quantile(x, 0.75, na.rm = TRUE))}
    )
  
  print(as.data.frame(summary_stats))
  cat("\n")
  
  return(summary_stats)
}

# Enhanced plotting function with new features
create_enhanced_plot <- function(df, opts, xcol, ycol, groupcol, use_log, log_base, summary_stats) {
  load_libraries()
  
  # Convert x column to factor for proper ordering
  df$.xcol_factor <- as.factor(df[[xcol]])
  
  # Base plot
  p <- ggplot(df, aes(x = .xcol_factor, y = .data[[ycol]]))
  
  # Add main plot element
  if (opts$type == "box") {
    p <- p + geom_boxplot(
      outlier.shape = if (opts$show_outliers) 19 else NA,
      width = 0.65,
      color = "#333333", 
      fill = NA, 
      alpha = 0.7,
      notch = opts$box_notch
    )
  } else if (opts$type == "violin") {
    p <- p + geom_violin(
      trim = opts$violin_trim, 
      width = 0.85,
      color = "#333333", 
      fill = NA, 
      alpha = 0.7
    )
    
    # Add median lines for violin plots
    if (opts$show_median) {
      p <- p + stat_summary(fun = median, geom = "crossbar", 
                           width = 0.5, color = "red", size = 0.3)
    }
  }
  
  # Add density plot overlay
  if (opts$density_plot && opts$type == "box") {
    p <- p + geom_violin(alpha = 0.3, color = NA, fill = "lightblue")
  }
  
  # Add confidence intervals
  if (opts$show_ci) {
    p <- p + stat_summary(
      fun.data = function(x) {
        m <- mean(x, na.rm = TRUE)
        se <- sd(x, na.rm = TRUE) / sqrt(length(x))
        ci <- qt(1 - (1 - opts$ci_level) / 2, length(x) - 1) * se
        data.frame(y = m, ymin = m - ci, ymax = m + ci)
      },
      geom = "errorbar", width = 0.2, color = "red", size = 1
    )
  }
  
  # Add mean points
  if (opts$show_mean) {
    p <- p + stat_summary(fun = mean, geom = "point", 
                         shape = 18, size = 4, color = "red")
  }
  
  # Add jittered points with enhanced coloring
  if (!is.null(groupcol)) {
    # Choose color palette
    if (opts$color_palette == "viridis" && requireNamespace("viridis", quietly = TRUE)) {
      p <- p + geom_jitter(aes(color = .data[[groupcol]]),
                          width = opts$jitter_width, height = 0,
                          size = opts$point_size, alpha = opts$point_alpha, shape = 16) +
               viridis::scale_color_viridis(discrete = TRUE)
    } else if (opts$color_palette == "brewer" && requireNamespace("RColorBrewer", quietly = TRUE)) {
      p <- p + geom_jitter(aes(color = .data[[groupcol]]),
                          width = opts$jitter_width, height = 0,
                          size = opts$point_size, alpha = opts$point_alpha, shape = 16) +
               scale_color_brewer(type = "qual", palette = "Set2")
    } else {
      p <- p + geom_jitter(aes(color = .data[[groupcol]]),
                          width = opts$jitter_width, height = 0,
                          size = opts$point_size, alpha = opts$point_alpha, shape = 16)
    }
  } else {
    p <- p + geom_jitter(width = opts$jitter_width, height = 0,
                        size = opts$point_size, alpha = opts$point_alpha, 
                        shape = 16, color = "steelblue")
  }
  
  # Add sample size annotations
  if (opts$show_n) {
    n_labels <- summary_stats %>%
      mutate(label = paste("n =", n))
    
    p <- p + geom_text(data = n_labels, 
                      aes(x = .data[[xcol]], y = max, label = label),
                      vjust = -0.5, size = 3, inherit.aes = FALSE)
  }
  
  # Add trend line
  if (opts$trend_line) {
    # Convert factor to numeric for trend
    df$.xcol_numeric <- as.numeric(df$.xcol_factor)
    p <- p + geom_smooth(aes(x = .xcol_numeric), method = "lm", 
                        se = TRUE, color = "orange", alpha = 0.3)
  }
  
  # Configure y-axis (log scale if requested)
  y_label <- ycol
  if (use_log) {
    y_range <- range(df[[ycol]], na.rm = TRUE)
    y_scale <- create_log_scale(y_range, log_base, opts)
    
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
    
    # Add log tick marks
    p <- p + annotation_logticks(sides = "l") + 
      coord_cartesian(clip = "off") +
      theme(plot.margin = margin(t = 10, r = 12, b = 10, l = 18))
  }
  
  # Apply theme
  base_theme <- switch(opts$theme_style,
    "bw" = theme_bw(base_size = 12),
    "minimal" = theme_minimal(base_size = 12),
    "classic" = theme_classic(base_size = 12),
    "dark" = theme_dark(base_size = 12),
    theme_bw(base_size = 12)
  )
  
  # Configure labels and theme
  title_text <- if (!is.na(opts$title)) opts$title else NULL
  
  p <- p +
    labs(x = xcol, y = y_label, 
         color = if (is.null(groupcol)) NULL else groupcol,
         title = title_text) +
    base_theme +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = if (is.null(groupcol)) "none" else "right"
    )
  
  # Rotate x-axis labels if requested
  if (isTRUE(opts$rotate_x)) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # Add faceting if specified
  if (!is.na(opts$facet_col) && opts$facet_col %in% names(df)) {
    p <- p + facet_wrap(~ .data[[opts$facet_col]], scales = "free_y")
  }
  
  return(p)
}

# Function to create interactive plot
make_interactive <- function(p, opts) {
  if (opts$interactive && requireNamespace("plotly", quietly = TRUE)) {
    if (opts$verbose) cat("Creating interactive plot...\n")
    return(plotly::ggplotly(p, tooltip = c("x", "y", "colour")))
  }
  return(p)
}

# Enhanced argument parsing with new options
parse_enhanced_args <- function(args) {
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
    
    # Parse key-value pairs
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
    
    # Normalize key name
    key_norm <- gsub("-", "_", key)
    
    # Extended key mapping
    key_mapping <- list(
      "i" = "input", "in" = "input", "x" = "xcol", "y" = "ycol", 
      "group" = "groupcol", "o" = "outfile", "out" = "outfile",
      "log10" = "log", "log_axis" = "log", "palette" = "color_palette",
      "theme" = "theme_style", "stats" = "stat_test", "v" = "verbose"
    )
    
    if (key_norm %in% names(key_mapping)) {
      key_norm <- key_mapping[[key_norm]]
    }
    
    # Assign value with type conversion
    if (key_norm %in% names(opts)) {
      opts[[key_norm]] <- convert_value(val, opts[[key_norm]])
    } else {
      warning(sprintf("Unknown option: --%s", key))
    }
    
    i <- i + 1
  }
  
  # Validate and set defaults for new options
  opts$type <- tolower(opts$type)
  if (!opts$type %in% c("box", "violin")) {
    warning("Invalid plot type, defaulting to 'box'")
    opts$type <- "box"
  }
  
  opts$color_palette <- tolower(opts$color_palette)
  if (!opts$color_palette %in% c("viridis", "brewer", "default")) {
    opts$color_palette <- "default"
  }
  
  opts$theme_style <- tolower(opts$theme_style)
  if (!opts$theme_style %in% c("bw", "minimal", "classic", "dark")) {
    opts$theme_style <- "bw"
  }
  
  # Generate output filename if not specified
  if (is.na(opts$outfile) || nchar(opts$outfile) == 0) {
    opts$outfile <- generate_enhanced_filename(opts)
  } else if (!grepl("\\.(pdf|png|jpg|jpeg|svg)$", tolower(opts$outfile))) {
    opts$outfile <- paste0(opts$outfile, ".pdf")
  }
  
  return(opts)
}

generate_enhanced_filename <- function(opts) {
  suffix <- ""
  if (isTRUE(opts$log)) suffix <- paste0(suffix, "_log")
  if (isTRUE(opts$show_mean)) suffix <- paste0(suffix, "_mean")
  if (isTRUE(opts$show_ci)) suffix <- paste0(suffix, "_ci")
  if (isTRUE(opts$pairwise)) suffix <- paste0(suffix, "_pairwise")
  
  base_name <- sprintf("%s_%s_by_%s%s", opts$ycol, opts$type, opts$xcol, suffix)
  return(paste0(base_name, ".pdf"))
}

# Enhanced help function
print_enhanced_usage <- function() {
  cat("\nEnhanced Box/Violin Plot with Advanced Statistical Analysis\n\n")
  cat("Usage: Rscript plot_box_enhanced.R [OPTIONS]\n\n")
  cat("Basic Options:\n")
  cat("  --input PATH           CSV input file (default: test.csv)\n")
  cat("  --xcol NAME            X-axis column name (default: Reproductive_Mode)\n")
  cat("  --ycol NAME            Y-axis column name (default: genomesize)\n")
  cat("  --groupcol NAME        Grouping column for point colors (default: group)\n")
  cat("  --type box|violin      Plot type (default: box)\n")
  cat("  --outfile PATH         Output file path (auto-generated if omitted)\n\n")
  cat("Visualization Options:\n")
  cat("  --show-mean [BOOL]     Show mean points (default: false)\n")
  cat("  --show-median [BOOL]   Show median lines (violin plots) (default: false)\n")
  cat("  --show-n [BOOL]        Show sample sizes (default: false)\n")
  cat("  --show-ci [BOOL]       Show confidence intervals (default: false)\n")
  cat("  --ci-level NUM         Confidence level (default: 0.95)\n")
  cat("  --show-outliers [BOOL] Show outliers in box plots (default: true)\n")
  cat("  --color-palette STR    Color palette: viridis|brewer|default (default: viridis)\n")
  cat("  --theme-style STR      Theme: bw|minimal|classic|dark (default: bw)\n")
  cat("  --box-notch [BOOL]     Add notches to box plots (default: false)\n")
  cat("  --violin-trim [BOOL]   Trim violin plot tails (default: true)\n")
  cat("  --density-plot [BOOL]  Add density overlay to box plots (default: false)\n")
  cat("  --trend-line [BOOL]    Add trend line (default: false)\n")
  cat("  --facet-col NAME       Column for faceting (optional)\n")
  cat("  --interactive [BOOL]   Create interactive plot (default: false)\n\n")
  cat("Statistical Options:\n")
  cat("  --stat-test STR        Statistical test: auto|anova|kruskal (default: auto)\n")
  cat("  --pairwise [BOOL]      Perform pairwise comparisons (default: false)\n")
  cat("  --effect-size [BOOL]   Calculate effect sizes (default: false)\n")
  cat("  --normality-test [BOOL] Test normality assumptions (default: false)\n")
  cat("  --export-stats [BOOL]  Export statistical results (default: false)\n\n")
  cat("Log Scale Options:\n")
  cat("  --log [BOOL]           Use log scale on y-axis (default: false)\n")
  cat("  --log-base 10|e        Log base (default: 10)\n")
  cat("  --log-label TYPE       Log labels: pow10|log10|plain (default: pow10)\n\n")
  cat("Other Options:\n")
  cat("  --width NUM            Figure width in inches (default: 8)\n")
  cat("  --height NUM           Figure height in inches (default: 5)\n")
  cat("  --dpi NUM              DPI for raster formats (default: 300)\n")
  cat("  --title TEXT           Plot title\n")
  cat("  --rotate-x [BOOL]      Rotate x-axis labels 45 degrees (default: false)\n")
  cat("  --verbose [BOOL]       Verbose output (default: false)\n")
  cat("  -h, --help             Show this help\n\n")
  cat("Examples:\n")
  cat("  # Basic enhanced plot with means and confidence intervals\n")
  cat("  Rscript plot_box_enhanced.R --input data.csv --show-mean --show-ci\n")
  cat("  # Comprehensive statistical analysis\n")
  cat("  Rscript plot_box_enhanced.R --input data.csv --pairwise --effect-size --normality-test\n")
  cat("  # Interactive violin plot with custom theme\n")
  cat("  Rscript plot_box_enhanced.R --type violin --interactive --theme minimal\n")
  cat("  # Faceted plot with sample sizes\n")
  cat("  Rscript plot_box_enhanced.R --facet-col treatment --show-n --export-stats\n\n")
}

# Utility functions (reuse from original script)
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
    warning(sprintf("Group column '%s' not found, points will not be colored by groups", groupcol))
    groupcol <- NULL
  }
  
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", "))) 
  }
  
  return(groupcol)
}

clean_numeric_data <- function(df, ycol) {
  suppressWarnings({
    df[[ycol]] <- as.numeric(df[[ycol]])
  })
  
  if (anyNA(df[[ycol]])) {
    warning(sprintf("Non-numeric values found in '%s' column, will be removed", ycol))
  }
  
  complete_idx <- !is.na(df[[ycol]])
  return(df[complete_idx, , drop = FALSE])
}

filter_positive_values <- function(df, ycol, use_log) {
  if (!use_log) return(df)
  
  positive_idx <- df[[ycol]] > 0
  df_filtered <- df[positive_idx, , drop = FALSE]
  
  n_removed <- sum(!positive_idx)
  if (n_removed > 0) {
    warning(sprintf("Removed %d rows with %s <= 0 for log transformation", n_removed, ycol))
  }
  
  if (nrow(df_filtered) == 0) {
    stop("No valid data points remaining after filtering for log transformation")
  }
  
  return(df_filtered)
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

create_log_scale <- function(data_range, log_base, opts) {
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

# Main enhanced function
main_enhanced <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- parse_enhanced_args(args)
  
  if (isTRUE(opts$help)) {
    print_enhanced_usage()
    quit(status = 0)
  }
  
  if (opts$verbose) cat("Starting enhanced box plot generation...\n")
  
  # Validate input file
  if (!file.exists(opts$input)) {
    stop(sprintf("Input file not found: %s", opts$input))
  }
  
  # Read and validate data
  if (opts$verbose) cat("Reading data...\n")
  df <- tryCatch({
    read.csv(opts$input, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }, error = function(e) {
    stop(sprintf("Failed to read CSV: %s", e$message))
  })
  
  # Validate columns
  opts$groupcol <- validate_columns(df, opts$xcol, opts$ycol, opts$groupcol)
  
  # Clean and prepare data
  if (opts$verbose) cat("Cleaning data...\n")
  df <- clean_numeric_data(df, opts$ycol)
  
  # Parse log settings
  use_log <- isTRUE(opts$log)
  log_base <- parse_log_base(opts$log_base)
  
  # Filter data for log transformation if needed
  df <- filter_positive_values(df, opts$ycol, use_log)
  
  # Enhanced data summary
  if (opts$verbose) cat("Generating data summary...\n")
  summary_stats <- enhanced_data_summary(df, opts$xcol, opts$ycol, opts$groupcol, opts)
  
  # Comprehensive statistical analysis
  if (opts$verbose) cat("Performing statistical tests...\n")
  stats_results <- perform_comprehensive_tests(df, opts$xcol, opts$ycol, opts)
  
  # Create enhanced plot
  if (opts$verbose) cat("Creating plot...\n")
  p <- create_enhanced_plot(df, opts, opts$xcol, opts$ycol, opts$groupcol, 
                           use_log, log_base, summary_stats)
  
  # Make interactive if requested
  p <- make_interactive(p, opts)
  
  # Save plot
  cat(sprintf("Saving plot to: %s\n", opts$outfile))
  if (opts$interactive && inherits(p, "plotly")) {
    # Save interactive plot as HTML
    html_file <- sub("\\.[^.]*$", ".html", opts$outfile)
    if (requireNamespace("htmlwidgets", quietly = TRUE)) {
      htmlwidgets::saveWidget(p, html_file)
      cat(sprintf("Interactive plot saved as: %s\n", html_file))
    } else {
      warning("htmlwidgets package not available, saving as static plot")
      ggsave(filename = opts$outfile, plot = p, width = opts$width, 
             height = opts$height, dpi = opts$dpi, units = "in")
    }
  } else {
    ggsave(filename = opts$outfile, plot = p, width = opts$width, 
           height = opts$height, dpi = opts$dpi, units = "in")
  }
  
  if (opts$verbose) cat("Analysis complete!\n")
  return(invisible(TRUE))
}

# Error handling wrapper
tryCatch({
  main_enhanced()
}, error = function(e) {
  cat(sprintf("Error: %s\n", e$message))
  print_enhanced_usage()
  quit(status = 1)
})
