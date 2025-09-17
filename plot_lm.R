#!/usr/bin/env Rscript

# linear regression plotting script
# Description: Creates scatter plots with linear regression lines and 95% CI

# Global constants
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
    warning(sprintf("Group column '%s' not found, ignoring grouping", groupcol))
    groupcol <- NULL
  }
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns: %s", paste(missing_cols, collapse = ", ")))
  }
  return(groupcol)
}

clean_numeric_data <- function(df, xcol, ycol) {
  suppressWarnings({
    df[[xcol]] <- as.numeric(df[[xcol]])
    df[[ycol]] <- as.numeric(df[[ycol]])
  })
  
  if (anyNA(df[[xcol]]) || anyNA(df[[ycol]])) {
    warning("Non-numeric values found in x or y columns, will be removed")
  }
  
  # Remove incomplete cases
  complete_idx <- stats::complete.cases(df[[xcol]], df[[ycol]])
  return(df[complete_idx, , drop = FALSE])
}

filter_positive_values <- function(df, xcol, ycol, to_log_x, to_log_y) {
  n_original <- nrow(df)
  
  if (to_log_x) {
    positive_x <- df[[xcol]] > 0
    df <- df[positive_x, , drop = FALSE]
    if (sum(!positive_x) > 0) {
      warning(sprintf("Removed %d rows with x <= 0 for log transformation", sum(!positive_x)))
    }
  }
  
  if (to_log_y) {
    positive_y <- df[[ycol]] > 0
    df <- df[positive_y, , drop = FALSE]
    if (sum(!positive_y) > 0) {
      warning(sprintf("Removed %d rows with y <= 0 for log transformation", sum(!positive_y)))
    }
  }
  
  if (nrow(df) == 0) {
    stop("No valid data points remaining after filtering")
  }
  
  return(df)
}

# Argument parsing (simplified and more robust)
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
      "in" = "input", "x" = "xcol", "y" = "ycol", "group" = "groupcol",
      "out" = "outfile", "log_axis" = "log", "log_base" = "log_base"
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
  
  # Generate output filename if not specified
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
  return(sprintf("lm_%s_vs_%s%s%s.pdf", opts$ycol, opts$xcol, suffix, bygroup_suffix))
}

# Log axis handling (consolidated function)
create_log_scale <- function(axis = "x", data_range, log_base, opts) {
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

# Model fitting and summary
fit_and_summarize_models <- function(dfc, xcol, ycol, groupcol, by_group, formula_str) {
  cat("===== Model Summary =====\n")
  
  # Overall model
  fit_overall <- lm(stats::as.formula(formula_str), data = dfc)
  cat("Overall model:\n")
  print(summary(fit_overall))
  cat("\n")
  
  # Group-specific models
  if (by_group && !is.null(groupcol)) {
    groups <- unique(dfc[[groupcol]])
    for (g in groups) {
      sub_data <- dfc[dfc[[groupcol]] == g, , drop = FALSE]
      if (nrow(sub_data) >= 2) {
        fit_group <- lm(stats::as.formula(formula_str), data = sub_data)
        cat(sprintf("Group %s = %s:\n", groupcol, g))
        print(summary(fit_group))
        cat("\n")
      } else {
        warning(sprintf("Insufficient data for group %s = %s", groupcol, g))
      }
    }
  }
}

# ggplot2 plotting function (optimized)
plot_with_ggplot2 <- function(dfc, opts, xcol, ycol, groupcol, by_group, 
                              to_log_x, to_log_y, log_base, formula_str) {
  suppressPackageStartupMessages(library(ggplot2))
  
  # Base plot
  p <- ggplot(dfc, aes(x = .data[[xcol]], y = .data[[ycol]]))
  
  # Add points
  if (!is.null(groupcol)) {
    p <- p + geom_point(aes(color = .data[[groupcol]]), shape = 16, size = 2.8, alpha = 0.9)
  } else {
    p <- p + geom_point(shape = 16, size = 2.8, alpha = 0.9)
  }
  
  # Add regression lines
  if (by_group && !is.null(groupcol)) {
    p <- p + geom_smooth(aes(color = .data[[groupcol]]), method = "lm", se = TRUE)
  } else {
    p <- p + geom_smooth(method = "lm", se = TRUE, color = "black")
  }
  
  # Configure axes
  p <- configure_axes(p, dfc, xcol, ycol, to_log_x, to_log_y, log_base, opts)
  
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
  
  # Save plot
  ggsave(filename = opts$outfile, plot = p, width = opts$width, 
         height = opts$height, dpi = opts$dpi, units = "in")
  
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

# Base R plotting function (simplified)
plot_with_base_r <- function(dfc, opts, xcol, ycol, groupcol, by_group,
                            to_log_x, to_log_y, log_base, formula_str) {
  grDevices::pdf(file = opts$outfile, width = opts$width, height = opts$height)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  x <- dfc[[xcol]]
  y <- dfc[[ycol]]
  
  # Configure colors
  if (!is.null(groupcol)) {
    grp <- factor(dfc[[groupcol]])
    cols <- grDevices::rainbow(length(levels(grp)))
    col_vec <- cols[as.integer(grp)]
  } else {
    grp <- NULL
    col_vec <- "black"
  }
  
  # Create base plot
  log_str <- paste0(if (to_log_x) "x" else "", if (to_log_y) "y" else "")
  title_text <- if (!is.na(opts$title)) opts$title else paste("Linear regression:", formula_str)
  
  graphics::plot(x, y, pch = 16, col = col_vec, 
                xlab = xcol, ylab = ycol, main = title_text, log = log_str)
  
  # Add regression lines and confidence intervals
  add_regression_lines(x, y, grp, cols, by_group, to_log_x, to_log_y, log_base)
  
  # Add legend if needed
  if (!is.null(grp)) {
    graphics::legend("topleft", legend = levels(grp), 
                    col = if (by_group) cols else unique(col_vec), 
                    pch = 16, bty = "n")
  }
}

add_regression_lines <- function(x, y, grp, cols, by_group, to_log_x, to_log_y, log_base) {
  logfun <- if (abs(log_base - exp(1)) < 1e-8) log else function(v) log(v, base = log_base)
  expfun <- if (abs(log_base - exp(1)) < 1e-8) exp else function(v) log_base^v
  
  fit_line <- function(xx, yy, col_line = "black") {
    x_fit <- if (to_log_x) logfun(xx) else xx
    y_fit <- if (to_log_y) logfun(yy) else yy
    
    fit <- lm(y_fit ~ x_fit)
    newx <- seq(min(xx, na.rm = TRUE), max(xx, na.rm = TRUE), length.out = 200)
    newx_fit <- if (to_log_x) logfun(newx) else newx
    
    pred <- predict(fit, newdata = data.frame(x_fit = newx_fit), interval = "confidence")
    y_draw <- if (to_log_y) expfun(pred[, "fit"]) else pred[, "fit"]
    y_lwr <- if (to_log_y) expfun(pred[, "lwr"]) else pred[, "lwr"]
    y_upr <- if (to_log_y) expfun(pred[, "upr"]) else pred[, "upr"]
    
    graphics::lines(newx, y_draw, col = col_line, lwd = 2)
    graphics::polygon(c(newx, rev(newx)), c(y_lwr, rev(y_upr)), 
                     border = NA, col = grDevices::adjustcolor(col_line, alpha.f = 0.2))
    
    return(fit)
  }
  
  if (by_group && !is.null(grp)) {
    levs <- levels(grp)
    for (i in seq_along(levs)) {
      idx <- which(grp == levs[i])
      if (length(idx) >= 2) {
        fit_line(x[idx], y[idx], cols[i])
      }
    }
  } else {
    fit_line(x, y)
  }
}

# Help function
print_usage <- function() {
  cat("
Scatter plot with linear regression and 95% confidence intervals.

Usage: Rscript plot_lm.R [OPTIONS]

Options:
  --input PATH           CSV input file (default: test2.csv)
  --xcol NAME            X-axis column name (default: Size_Value_cm)
  --ycol NAME            Y-axis column name (default: genomesize)
  --groupcol NAME        Grouping column for colors (default: group)
  --by-group [BOOL]      Fit separate lines per group (default: false)
  --log none|x|y|xy      Apply log transform to axes (default: none)
  --log-base 10|e        Log base (default: 10)
  --log-label TYPE       Log labels: pow10|log10|plain (default: pow10)
  --log-step NUM         Log break step in decades (default: 1)
  --log-digits NUM       Digits for log10 labels (default: 1)
  --outfile PATH         Output PDF path (auto-generated if omitted)
  --width NUM            Figure width in inches (default: 7)
  --height NUM           Figure height in inches (default: 5.5)
  --dpi NUM              DPI for embedded rasters (default: 300)
  --title TEXT           Plot title
  -h, --help             Show this help

Examples:
  Rscript plot_lm.R --input data.csv --xcol size --ycol weight
  Rscript plot_lm.R --input data.csv --log xy --by-group true
  Rscript plot_lm.R --input data.csv --log-base e --log-label log10

")
}

# Main function (simplified and more robust)
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
  df <- read.csv(opts$input, stringsAsFactors = FALSE, check.names = FALSE)
  opts$groupcol <- validate_columns(df, opts$xcol, opts$ycol, opts$groupcol)
  
  # Clean and prepare data
  dfc <- clean_numeric_data(df, opts$xcol, opts$ycol)
  
  # Determine log transformations
  to_log_x <- tolower(opts$log) %in% c("x", "xy", "both")
  to_log_y <- tolower(opts$log) %in% c("y", "xy", "both")
  log_base <- parse_log_base(opts$log_base)
  
  # Filter data for log transformations
  dfc <- filter_positive_values(dfc, opts$xcol, opts$ycol, to_log_x, to_log_y)
  
  # Create formula string for display
  logfun_name <- if (abs(log_base - exp(1)) < 1e-8) "log" else "log10"
  lhs <- if (to_log_y) paste0(logfun_name, "(", opts$ycol, ")") else opts$ycol
  rhs <- if (to_log_x) paste0(logfun_name, "(", opts$xcol, ")") else opts$xcol
  formula_str <- paste(lhs, "~", rhs)
  
  # Create plot
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_with_ggplot2(dfc, opts, opts$xcol, opts$ycol, opts$groupcol, 
                      isTRUE(opts$by_group), to_log_x, to_log_y, log_base, formula_str)
  } else {
    message("ggplot2 not available, using base R graphics")
    plot_with_base_r(dfc, opts, opts$xcol, opts$ycol, opts$groupcol,
                     isTRUE(opts$by_group), to_log_x, to_log_y, log_base, formula_str)
  }
  
  # Print model summaries
  fit_and_summarize_models(dfc, opts$xcol, opts$ycol, opts$groupcol, 
                          isTRUE(opts$by_group), formula_str)
  
  cat(sprintf("Plot saved to: %s\n", opts$outfile))
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
