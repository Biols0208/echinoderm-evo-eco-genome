# ============================================================
# Data Processing Module for PGLS Analysis Framework
# Handles data loading, cleaning, transformation, and tree processing
# ============================================================

source("R/utils.R")

# Load data files with configurable options
load_data <- function(config) {
  if (config$output$verbose) cat("Loading data files...\n")
  
  # Load tree
  tree <- tryCatch({
    read.tree(config$input$tree_file)
  }, error = function(e) {
    stop(sprintf("Error loading tree file '%s': %s", config$input$tree_file, e$message))
  })
  
  # Load data with configurable options
  data <- tryCatch({
    if (config$input$data_format == "xlsx" && requireNamespace("readxl", quietly = TRUE)) {
      readxl::read_excel(config$input$data_file)
    } else {
      read.csv(config$input$data_file, 
               na.strings = c("", "NA"), 
               sep = config$input$data_separator,
               dec = config$input$data_decimal,
               fileEncoding = config$input$data_encoding, 
               stringsAsFactors = FALSE)
    }
  }, error = function(e) {
    stop(sprintf("Error loading data file '%s': %s", config$input$data_file, e$message))
  })
  
  if (config$output$verbose) {
    cat(sprintf("Tree loaded: %d tips\n", length(tree$tip.label)))
    cat(sprintf("Data loaded: %d rows, %d columns\n", nrow(data), ncol(data)))
  }
  
  return(list(tree = tree, data = data))
}

# Detect outliers using different methods
detect_outliers <- function(x, method = "IQR", threshold = 3) {
  outliers <- switch(method,
    "IQR" = {
      Q1 <- quantile(x, 0.25, na.rm = TRUE)
      Q3 <- quantile(x, 0.75, na.rm = TRUE)
      IQR <- Q3 - Q1
      (x < Q1 - threshold * IQR) | (x > Q3 + threshold * IQR)
    },
    "zscore" = {
      z <- abs((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
      z > threshold
    },
    "mad" = {
      med <- median(x, na.rm = TRUE)
      mad <- median(abs(x - med), na.rm = TRUE)
      abs(x - med) / mad > threshold
    },
    rep(FALSE, length(x))
  )
  
  return(outliers)
}

# Comprehensive data preparation function
prepare_analysis_data <- function(data, config) {
  if (config$output$verbose) cat("\n=== Data Preparation ===\n")
  
  # Check required columns
  required_cols <- c(config$variables$species_col, 
                    config$variables$response_vars, 
                    config$variables$predictor_vars)
  
  if (!is.null(config$variables$group_var)) {
    required_cols <- c(required_cols, config$variables$group_var)
  }
  
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns in data: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Prepare species names
  data[[config$variables$species_col]] <- gsub(" ", "_", trimws(data[[config$variables$species_col]]))
  
  # Handle predictors (continuous vs categorical)
  for (i in seq_along(config$variables$predictor_vars)) {
    pred <- config$variables$predictor_vars[i]
    
    if (i <= length(config$variables$predictor_is_continuous) && 
        config$variables$predictor_is_continuous[i]) {
      # Continuous predictor
      data[[pred]] <- as.numeric(data[[pred]])
      
      if (config$preprocessing$center_predictors) {
        data[[pred]] <- scale(data[[pred]], center = TRUE, scale = FALSE)[,1]
        if (config$output$verbose) cat(sprintf("Centered continuous predictor: %s\n", pred))
      }
      
      if (config$preprocessing$scale_predictors) {
        data[[pred]] <- scale(data[[pred]], center = FALSE, scale = TRUE)[,1]
        if (config$output$verbose) cat(sprintf("Scaled continuous predictor: %s\n", pred))
      }
    } else {
      # Categorical predictor
      data[[pred]] <- factor(data[[pred]])
      levels_count <- length(levels(data[[pred]]))
      
      if (config$output$verbose) {
        cat(sprintf("Predictor '%s': %d levels\n", pred, levels_count))
      }
    }
  }
  
  # Convert group variable to factor if provided
  if (!is.null(config$variables$group_var)) {
    data[[config$variables$group_var]] <- factor(data[[config$variables$group_var]])
    
    # Filter groups if specified
    if (!is.null(config$variables$groups_to_analyze)) {
      data <- data[data[[config$variables$group_var]] %in% config$variables$groups_to_analyze, ]
      data[[config$variables$group_var]] <- droplevels(data[[config$variables$group_var]])
    }
    
    if (config$output$verbose) {
      cat(sprintf("Groups in '%s': %s\n", config$variables$group_var, 
                  paste(levels(data[[config$variables$group_var]]), collapse = ", ")))
    }
  }
  
  # Handle outliers in response variables
  if (config$preprocessing$detect_outliers) {
    for (resp in config$variables$response_vars) {
      outliers <- detect_outliers(data[[resp]], 
                                 config$preprocessing$outlier_method,
                                 config$preprocessing$outlier_threshold)
      n_outliers <- sum(outliers, na.rm = TRUE)
      
      if (n_outliers > 0) {
        if (config$output$verbose) {
          cat(sprintf("Detected %d outliers in '%s'\n", n_outliers, resp))
        }
        
        if (config$preprocessing$remove_outliers) {
          data[[resp]][outliers] <- NA
          if (config$output$verbose) {
            cat(sprintf("Removed %d outliers from '%s'\n", n_outliers, resp))
          }
        }
      }
    }
  }
  
  # Log transform response variables if requested
  if (config$preprocessing$log_transform_response) {
    for (resp in config$variables$response_vars) {
      log_name <- sprintf("%s_log%g", resp, config$preprocessing$log_base)
      
      if (all(data[[resp]] > 0, na.rm = TRUE)) {
        if (config$preprocessing$log_base == 10) {
          data[[log_name]] <- log10(data[[resp]])
        } else if (config$preprocessing$log_base == 2) {
          data[[log_name]] <- log2(data[[resp]])
        } else if (config$preprocessing$log_base == exp(1)) {
          data[[log_name]] <- log(data[[resp]])
        } else {
          data[[log_name]] <- log(data[[resp]]) / log(config$preprocessing$log_base)
        }
        
        if (config$output$verbose) {
          cat(sprintf("Created log%g transformed variable: %s\n", 
                     config$preprocessing$log_base, log_name))
        }
      } else {
        warning(sprintf("Cannot log transform '%s' - contains non-positive values", resp))
      }
    }
  }
  
  # Remove NAs if requested
  initial_n <- nrow(data)
  
  if (config$preprocessing$exclude_na_response) {
    for (resp in config$variables$response_vars) {
      data <- data[!is.na(data[[resp]]), ]
    }
  }
  
  if (config$preprocessing$exclude_na_predictor) {
    for (pred in config$variables$predictor_vars) {
      data <- data[!is.na(data[[pred]]), ]
    }
  }
  
  final_n <- nrow(data)
  if (config$output$verbose && final_n < initial_n) {
    cat(sprintf("Removed %d rows with NA values\n", initial_n - final_n))
  }
  
  return(data)
}

# Clean tree by removing specified species
clean_tree <- function(tree, exclude_species = NULL) {
  if (is.null(exclude_species) || length(exclude_species) == 0) {
    return(tree)
  }
  
  cat(sprintf("Removing %d species from tree...\n", length(exclude_species)))
  tree_clean <- tryCatch({
    drop.tip(tree, exclude_species)
  }, error = function(e) {
    warning(sprintf("Could not remove species: %s", e$message))
    tree
  })
  
  tree_clean$node.label <- NULL
  return(tree_clean)
}

# Create subset trees for groups
create_group_trees <- function(tree, data, config) {
  group_var <- config$variables$group_var
  if (is.null(group_var)) return(NULL)
  
  if (config$output$verbose) cat("\nCreating trees for each group...\n")
  
  groups <- split(data[[config$variables$species_col]], data[[group_var]])
  trees <- list()
  
  # Filter by groups_to_analyze if specified
  if (!is.null(config$variables$groups_to_analyze)) {
    groups <- groups[names(groups) %in% config$variables$groups_to_analyze]
  }
  
  for (group_name in names(groups)) {
    species_in_group <- groups[[group_name]]
    species_to_keep <- intersect(tree$tip.label, species_in_group)
    
    if (length(species_to_keep) >= config$variables$min_group_size) {
      trees[[group_name]] <- tryCatch({
        drop.tip(tree, setdiff(tree$tip.label, species_to_keep))
      }, error = function(e) {
        warning(sprintf("Could not create tree for %s: %s", group_name, e$message))
        NULL
      })
      
      if (!is.null(trees[[group_name]])) {
        trees[[group_name]]$node.label <- NULL
        if (config$output$verbose) {
          cat(sprintf("  %s: %d species\n", group_name, length(species_to_keep)))
        }
      }
    } else {
      if (config$output$verbose) {
        cat(sprintf("  %s: Too few species (%d < %d)\n", 
                   group_name, length(species_to_keep), config$variables$min_group_size))
      }
    }
  }
  
  return(trees)
}