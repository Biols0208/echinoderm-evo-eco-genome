#!/usr/bin/env Rscript

# box/violin plot script with flexible column support
# Description: Creates box plots or violin plots with overlaid points colored by groups

# Global constants
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
  help = FALSE
)

# Utility functions
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
  
  # Remove rows with NA in y column
  complete_idx <- !is.na(df[[ycol]])
  return(df[complete_idx, , drop = FALSE])
}

filter_positive_values <- function(df, ycol, use_log) {
  if (!use_log) return(df)
  
  n_original <- nrow(df)
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

# Argument parsing (unified with plot_lm.R style)
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
    
    # Map common aliases
    key_mapping <- list(
      "i" = "input", "in" = "input", "x" = "xcol", "y" = "ycol", 
      "group" = "groupcol", "o" = "outfile", "out" = "outfile",
      "log10" = "log", "log_axis" = "log"
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
  
  # Validate plot type
  opts$type <- tolower(opts$type)
  if (!opts$type %in% c("box", "violin")) {
    warning("Invalid plot type, defaulting to 'box'")
    opts$type <- "box"
  }
  
  # Generate output filename if not specified
  if (is.na(opts$outfile) || nchar(opts$outfile) == 0) {
    opts$outfile <- generate_output_filename(opts)
  } else if (!grepl("\\.(pdf|png|jpg|jpeg|svg)$", tolower(opts$outfile))) {
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
  suffix <- if (isTRUE(opts$log)) "_log10" else ""
  base_name <- sprintf("%s_%s_by_%s%s", opts$ycol, opts$type, opts$xcol, suffix)
  return(paste0(base_name, ".pdf"))
}

# Log scale configuration (shared with plot_lm.R)
create_log_scale <- function(data_range, log_base, opts) {
  step <- if (is.na(opts$log_step) || opts$log_step <= 0) 1.0 else opts$log_step
  
  # Generate breaks
  if (abs(log_base - 10) < 1e-8) {
    lo <- floor(log10(data_range[1]))
    hi <- ceiling(log10(data_range[2]))
    breaks <- 10^(seq(lo, hi, by = step))
  } else {
    lo <- floor(log(data_range[1], base = log_base))
    hi <- ceiling(log(data_range[2], base = log_base))
    breaks <- log_base^(seq(lo, hi, by = step))
  }
  
  # Generate labels
  labels <- create_log_labels(breaks, log_base, opts)
  
  # Generate axis label text
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
    breaks  # plain labels
  )
}

# Data summary function
print_data_summary <- function(df, xcol, ycol, groupcol) {
  cat("===== Data Summary =====\n")
  cat(sprintf("Total observations: %d\n", nrow(df)))
  cat(sprintf("X variable (%s): %d unique values\n", xcol, length(unique(df[[xcol]]))))
  cat(sprintf("Y variable (%s): range [%.2e, %.2e]\n", ycol, min(df[[ycol]], na.rm = TRUE), max(df[[ycol]], na.rm = TRUE)))
  
  if (!is.null(groupcol)) {
    groups <- unique(df[[groupcol]])
    cat(sprintf("Groups (%s): %d levels - %s\n", groupcol, length(groups), paste(groups, collapse = ", ")))
  }
  
  # Summary by x categories
  cat("\nSummary by categories:\n")
  x_summary <- aggregate(df[[ycol]], by = list(df[[xcol]]), 
                        FUN = function(x) c(n = length(x), mean = mean(x, na.rm = TRUE), 
                                          median = median(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))
  colnames(x_summary) <- c(xcol, "statistics")
  print(x_summary)
  cat("\n")
}

# Main plotting function
create_plot <- function(df, opts, xcol, ycol, groupcol, use_log, log_base) {
  suppressPackageStartupMessages(library(ggplot2))
  
  # Convert x column to factor for proper ordering
  df$.xcol_factor <- as.factor(df[[xcol]])
  
  # Base plot
  p <- ggplot(df, aes(x = .xcol_factor, y = .data[[ycol]]))
  
  # Add box plot or violin plot
  if (opts$type == "box") {
    p <- p + geom_boxplot(outlier.shape = NA, width = 0.65, 
                         color = "#333333", fill = NA, alpha = 0.7)
  } else if (opts$type == "violin") {
    p <- p + geom_violin(trim = FALSE, width = 0.85, 
                        color = "#333333", fill = NA, alpha = 0.7)
  }
  
  # Add jittered points
  if (!is.null(groupcol)) {
    p <- p + geom_jitter(aes(color = .data[[groupcol]]),
                        width = opts$jitter_width, height = 0,
                        size = opts$point_size, alpha = opts$point_alpha, shape = 16)
  } else {
    p <- p + geom_jitter(width = opts$jitter_width, height = 0,
                        size = opts$point_size, alpha = opts$point_alpha, shape = 16)
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
      # Fallback to log10
      p <- p + scale_y_log10(breaks = y_scale$breaks, labels = y_scale$labels,
                            expand = expansion(mult = c(0.03, 0.06)))
    }
    
    y_label <- paste0(ycol, y_scale$axis_label)
    
    # Add log tick marks
    p <- p + annotation_logticks(sides = "l") + 
      coord_cartesian(clip = "off") +
      theme(plot.margin = margin(t = 10, r = 12, b = 10, l = 18))
  }
  
  # Configure theme and labels
  title_text <- if (!is.na(opts$title)) opts$title else NULL
  
  p <- p +
    labs(x = xcol, y = y_label, 
         color = if (is.null(groupcol)) NULL else groupcol,
         title = title_text) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Rotate x-axis labels if requested
  if (isTRUE(opts$rotate_x)) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  return(p)
}

# Statistical tests (optional enhancement)
perform_statistical_tests <- function(df, xcol, ycol) {
  cat("===== Statistical Tests =====\n")
  
  # Check if we have enough groups for testing
  x_groups <- unique(df[[xcol]])
  if (length(x_groups) < 2) {
    cat("Only one group found, no statistical tests performed.\n\n")
    return()
  }
  
  # Kruskal-Wallis test (non-parametric)
  tryCatch({
    kw_test <- kruskal.test(df[[ycol]], df[[xcol]])
    cat(sprintf("Kruskal-Wallis test: H = %.3f, p-value = %.3e\n", 
                kw_test$statistic, kw_test$p.value))
  }, error = function(e) {
    cat("Kruskal-Wallis test failed:", e$message, "\n")
  })
  
  # ANOVA (parametric, if assumptions are met)
  tryCatch({
    aov_test <- aov(df[[ycol]] ~ df[[xcol]])
    aov_summary <- summary(aov_test)
    f_stat <- aov_summary[[1]][["F value"]][1]
    p_val <- aov_summary[[1]][["Pr(>F)"]][1]
    cat(sprintf("ANOVA: F = %.3f, p-value = %.3e\n", f_stat, p_val))
  }, error = function(e) {
    cat("ANOVA test failed:", e$message, "\n")
  })
  
  cat("\n")
}

# Help function
print_usage <- function() {
  cat("
Box plot or violin plot with overlaid points colored by groups.

Usage: Rscript plot_box.R [OPTIONS]

Options:
  --input PATH           CSV input file (default: test.csv)
  --xcol NAME            X-axis column name (default: Reproductive_Mode)
  --ycol NAME            Y-axis column name (default: genomesize)
  --groupcol NAME        Grouping column for point colors (default: group)
  --type box|violin      Plot type (default: box)
  --log [BOOL]           Use log scale on y-axis (default: false)
  --log-base 10|e        Log base (default: 10)
  --log-label TYPE       Log labels: pow10|log10|plain (default: pow10)
  --log-step NUM         Log break step in decades (default: 1)
  --log-digits NUM       Digits for log10 labels (default: 1)
  --outfile PATH         Output file path (auto-generated if omitted)
  --width NUM            Figure width in inches (default: 8)
  --height NUM           Figure height in inches (default: 5)
  --dpi NUM              DPI for raster formats (default: 300)
  --title TEXT           Plot title
  --rotate-x [BOOL]      Rotate x-axis labels 45 degrees (default: false)
  --jitter-width NUM     Jitter width for points (default: 0.15)
  --point-size NUM       Point size (default: 2)
  --point-alpha NUM      Point transparency (default: 0.8)
  -h, --help             Show this help

Examples:
  Rscript plot_box.R --input data.csv --xcol treatment --ycol response
  Rscript plot_box.R --input data.csv --type violin --log true
  Rscript plot_box.R --input data.csv --xcol species --rotate-x true

")
}

# Main function
main <- function() {
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
  
  # Read and validate data
  df <- tryCatch({
    read.csv(opts$input, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }, error = function(e) {
    stop(sprintf("Failed to read CSV: %s", e$message))
  })
  
  # Validate columns
  opts$groupcol <- validate_columns(df, opts$xcol, opts$ycol, opts$groupcol)
  
  # Clean and prepare data
  df <- clean_numeric_data(df, opts$ycol)
  
  # Parse log settings
  use_log <- isTRUE(opts$log)
  log_base <- parse_log_base(opts$log_base)
  
  # Filter data for log transformation if needed
  df <- filter_positive_values(df, opts$ycol, use_log)
  
  # Print data summary
  print_data_summary(df, opts$xcol, opts$ycol, opts$groupcol)
  
  # Perform statistical tests
  perform_statistical_tests(df, opts$xcol, opts$ycol)
  
  # Create and save plot
  p <- create_plot(df, opts, opts$xcol, opts$ycol, opts$groupcol, use_log, log_base)
  
  cat(sprintf("Saving plot to: %s\n", opts$outfile))
  ggsave(filename = opts$outfile, plot = p, width = opts$width, 
         height = opts$height, dpi = opts$dpi, units = "in")
  
  return(invisible(TRUE))
}

# Error handling wrapper
tryCatch({
  main()
}, error = function(e) {
  cat(sprintf("Error: %s\n", e$message))
  print_usage()
  quit(status = 1)
})
