# This script contains the auxiliary functions necessary to explore and process
# the datasets in this project.

# ==============================================================================
# UTILS FOR "1. Setup & data ingestion"
# ==============================================================================

#' @title Downloads and types datasets from Kaggle
#'
#' @description Downloads a public dataset (ZIP format) from Kaggle via URL,
#' unzips it in a temporary folder, dynamically searches for the CSV file
#' regardless of subfolders, moves it to the project directory, and types it.
#'
#' @param direct_url String. Direct link to the ZIP file on Kaggle.
#' @param output_name String. Base name to save the final file.
#' @param dataset_type String. Logical identifier ("cvd" or "cancer").
#' @param dest_dir String. Destination folder relative to the notebook (default "../data/raw").
#' @param kaggle_user String. Kaggle username for API authentication.
#' @param kaggle_key String. Kaggle API key (from kaggle.json).
#'
#' @return A dataframe with the loaded and typed data.
download_and_type_kaggle <- function(direct_url, output_name, dataset_type, dest_dir = "../data/raw", kaggle_user = NULL, kaggle_key = NULL) {
  
  # Create the destination directory if it does not exist yet
  if (!base::dir.exists(dest_dir)) {
    base::dir.create(dest_dir, recursive = TRUE)
  }
  
  # Define the final CSV path
  csv_path <- base::file.path(dest_dir, base::paste0(output_name, ".csv"))
  
  # Start the download only if the file does not exist locally
  if (!base::file.exists(csv_path)) {
    base::message(base::paste("Downloading and extracting dataset:", output_name))
    
    # Create temporary paths in the system
    temp_zip <- base::tempfile(fileext = ".zip")
    temp_dir <- base::file.path(base::tempdir(), output_name)
    base::dir.create(temp_dir, showWarnings = FALSE)
    
    # Perform authenticated HTTP GET request using httr to handle Kaggle's API requirements
    if (!base::is.null(kaggle_user) && !base::is.null(kaggle_key)) {
      if (!base::requireNamespace("httr", quietly = TRUE)) {
        utils::install.packages("httr")
      }
      
      response <- httr::GET(
        url = direct_url,
        httr::authenticate(user = kaggle_user, password = kaggle_key, type = "basic"),
        httr::config(followlocation = TRUE),
        httr::write_disk(temp_zip, overwrite = TRUE)
      )
      
      if (httr::http_error(response)) {
        base::unlink(temp_zip)
        base::unlink(temp_dir, recursive = TRUE)
        base::stop(base::paste("HTTP download failed with status:", httr::status_code(response)))
      }
    } else {
      # Fallback to base download if no credentials are provided
      utils::download.file(url = direct_url, destfile = temp_zip, mode = "wb", quiet = TRUE)
    }
    
    # Unzip the file into the temporary folder
    extracted_files <- utils::unzip(temp_zip, exdir = temp_dir)
    
    # Dynamically search for the CSV file, regardless of the subfolder it is in
    extracted_csv_path <- extracted_files[base::grepl("\\.csv$", extracted_files, ignore.case = TRUE)]
    
    # Stop the process if no CSV file was found inside the ZIP
    if (base::length(extracted_csv_path) == 0) {
      base::unlink(temp_zip)
      base::unlink(temp_dir, recursive = TRUE)
      base::stop("Critical error: no .csv file was found inside the ZIP.")
    }
    
    # Move the CSV to the final destination and rename it
    base::file.copy(from = extracted_csv_path[1], to = csv_path, overwrite = TRUE)
    
    # Clean up the temporary files
    base::unlink(temp_zip)
    base::unlink(temp_dir, recursive = TRUE)
    
  } else {
    base::message(base::paste("The file", output_name, "already exists locally. Skipping download."))
  }
  
  # The Kaggle cancer CSV has a header with 33 columns (the 33rd is empty) but
  # all data rows have only 32 fields. readr interprets this mismatch as a
  # parsing error and silently discards the last row of the file. The fix is
  # to read the header manually, discard the empty column, and pass the clean
  # names to readr ourselves.
  if (dataset_type == "cancer") {
    
    # Read only the first line to extract the real column names
    raw_header <- readr::read_lines(csv_path, n_max = 1)
    clean_col_names <- base::strsplit(raw_header, ",")[[1]] |>
      base::trimws() |>
      # Remove quotes and drop the trailing empty column
      (\(x) x[base::nchar(base::gsub('"', '', x)) > 0])() |>
      (\(x) base::gsub('"', '', x))()
    
    # Read the file skipping the original header and using the clean names instead
    df <- readr::read_csv(
      csv_path,
      col_names       = clean_col_names,
      skip            = 1,
      show_col_types  = FALSE,
      skip_empty_rows = TRUE
    )
    
  } else {
    df <- readr::read_csv(csv_path, show_col_types = FALSE, skip_empty_rows = TRUE)
  }
  
  # Format all column names to snake_case
  df <- janitor::clean_names(df)
  
  # Type the variables according to the nature of the dataset
  if (dataset_type == "cvd") {
    df <- df |>
      dplyr::mutate(
        sex                     = base::as.factor(sex),
        blood_pressure_category = base::as.factor(blood_pressure_category),
        smoking_status          = base::as.factor(smoking_status),
        diabetes_status         = base::as.factor(diabetes_status),
        physical_activity_level = base::as.factor(physical_activity_level),
        family_history_of_cvd   = base::as.factor(family_history_of_cvd),
        cvd_risk_level          = base::as.factor(cvd_risk_level),
        blood_pressure_mm_hg    = base::as.character(blood_pressure_mm_hg)
      )
  } else if (dataset_type == "cancer") {
    df <- df |>
      dplyr::mutate(
        id        = base::as.character(id),
        diagnosis = base::as.factor(diagnosis)
      ) |>
      # This select() is no longer strictly necessary, but we keep it as a safety net
      dplyr::select(-dplyr::any_of(c("...33", "unnamed_32", "x33")))
  }
  
  base::return(df)
}

# ==============================================================================
# UTILS FOR "2. Initial cleaning"
# ==============================================================================

#' @title Pre-filters domain variables and negative values
#'
#' @description Removes specific columns based purely on the domain knowledge
#' of the dataset, to avoid data leakage and unnecessary dimensional redundancy.
#' Additionally applies a biological impossibility filter: since all raw
#' numeric variables (ages, geometric distances, concentrations) have a lower
#' physiological limit of zero, any negative value is forced to null (NA).
#'
#' @param df Dataframe. The raw dataset to process.
#' @param dataset_type String. Logical identifier ("cvd" or "cancer") to apply the culling rules.
#'
#' @return A dataframe purged of structural noise and negative artifacts.
remove_domain_variables <- function(df, dataset_type) {
  
  # Remove CVD-specific variables
  if (dataset_type == "cvd") {
    result <- df |>
      dplyr::select(-dplyr::any_of(c(
        "cvd_risk_score",
        "blood_pressure_category",
        "blood_pressure_mm_hg",
        "height_cm"
      )))
    
    # Remove Cancer-specific variables
  } else if (dataset_type == "cancer") {
    result <- df |>
      dplyr::select(-dplyr::any_of(c("id")))
    
  } else {
    result <- df
  }
  
  # Iterate over every numeric column in the resulting dataframe
  result <- result |>
    dplyr::mutate(dplyr::across(
      tidyselect::where(base::is.numeric),
      # If the value is negative, convert it to NA; if positive or zero, keep it as is
      ~ base::ifelse(.x < 0, NA, .x)
    ))
  
  base::return(result)
}

# ------------------------------------------------------------------------------

#' @title Body Mass Index (BMI) sanity check and dimensional reduction
#'
#' @description Evaluates the mathematical congruence between weight, height,
#' and the reported BMI. Corrects inconsistencies, calculates the BMI when the
#' original value is missing but the base metrics exist, and removes the
#' originating columns (`weight_kg` and `height_m`) due to redundancy.
#'
#' @param df Dataframe. The dataset to audit.
#' @param tolerance Numeric. Acceptable discrepancy threshold between the calculated
#' and the reported BMI (default 0.5).
#' @param correct Logical. If TRUE, corrects inconsistent BMIs and fills
#' calculable NAs; if FALSE, only fills NAs without modifying the original values.
#'
#' @return A dataframe with a consolidated BMI column and without the base metrics.
verify_bmi_consistency <- function(df, tolerance = 0.5, correct = TRUE) {
  
  # Extract the relevant columns
  weight <- df$weight_kg
  height <- df$height_m
  original_bmi <- df$bmi
  
  # Calculate the theoretical BMI using the standard formula
  calculated_bmi <- weight / (height^2)
  
  # If the original BMI is NA but can be calculated from weight/height, rescue it
  recovered_bmi <- dplyr::coalesce(original_bmi, calculated_bmi)
  
  # Determine the absolute difference between the recovered/original and the calculated BMI
  absolute_difference <- base::abs(recovered_bmi - calculated_bmi)
  
  # Count the real discrepancies that exceed the tolerance
  discrepancies <- base::sum(absolute_difference > tolerance, na.rm = TRUE)
  total_valid <- base::sum(!base::is.na(absolute_difference))
  
  # Print the analytical report to the console
  base::cat("\n--- BMI consistency audit ---\n")
  base::cat(base::sprintf("Evaluated records (with sufficient data): %d\n", total_valid))
  base::cat(base::sprintf("Discrepancies detected (> %s): %d\n", tolerance, discrepancies))
  
  # Manage the BMI column and impute missing values accordingly
  if (correct) {
    base::cat("Action: correcting inconsistencies and filling calculable missing BMIs.\n\n")
    # If the calculation fails due to missing weight/height (NA), keep the original value
    df$bmi <- base::ifelse(base::is.na(calculated_bmi), original_bmi, calculated_bmi)
  } else {
    base::cat("Action: only filling NAs where possible, without overriding original values.\n\n")
    df$bmi <- recovered_bmi
  }
  
  # Remove the weight and height columns to avoid redundancy and multicollinearity
  df <- df |>
    dplyr::select(-dplyr::any_of(c("weight_kg", "height_m")))
  
  base::return(df)
}

# ==============================================================================
# UTILS FOR "4. Initial exploratory data analysis"
# ==============================================================================

#' @title Visualizes numeric densities
#'
#' @description Extracts the numeric variables from a dataframe and generates a
#' panel with the density plot of each of them.
#'
#' @param df Dataframe. The dataset to visualize.
#' @param dataset_name String. Dataset name used to customize the plot title.
#' @param fill_color String. HEX code for the plot fill color.
#'
#' @return A ggplot2 object.
plot_numeric_densities <- function(df, dataset_name = "", fill_color = "#2980b9") {
  
  # Keep only the numeric variables
  df_num <- df |>
    dplyr::select(tidyselect::where(base::is.numeric))
  
  if (base::ncol(df_num) == 0) {
    base::return(base::message("No numeric variables to visualize."))
  }
  
  # Transform the dataframe to long format so ggplot2 can facet by variable
  df_long <- tidyr::pivot_longer(
    df_num, cols = dplyr::everything(), names_to = "variable", values_to = "value"
  ) |> tidyr::drop_na()
  
  # Build the dynamic plot title
  plot_title <- base::paste0(
    "Density curves - Continuous variables",
    base::ifelse(dataset_name == "", "", base::paste0(" (", dataset_name, ")"))
  )
  
  # Build the density plot with one facet per variable
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = value)) +
    ggplot2::geom_density(fill = fill_color, alpha = 0.7, color = "#1abc9c") +
    ggplot2::facet_wrap(~ variable, scales = "free", ncol = 5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Value",
      y = "Density",
      title = plot_title
    ) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "#ecf0f1", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(face = "bold", size = 14, margin = ggplot2::margin(b = 15))
    )
  
  base::return(p)
}

# ------------------------------------------------------------------------------

#' @title Visualizes categorical variables (bar charts)
#'
#' @description Extracts the categorical variables from a dataframe and generates
#' a grid of ordered bar charts, each showing the percentage and count per category.
#'
#' @param df Dataframe. The dataset to visualize.
#' @param dataset_name String. Dataset name used to customize the plot title.
#'
#' @return A ggplot2 object.
plot_categorical_bars <- function(df, dataset_name = "") {
  
  # Keep only the categorical variables (factor or character)
  df_cat <- df |>
    dplyr::select(tidyselect::where(~ base::is.factor(.x) || base::is.character(.x)))
  
  if (base::ncol(df_cat) == 0) {
    base::return(base::message("No categorical variables to visualize in this dataset."))
  }
  
  # Transform the dataframe to long format so ggplot2 can facet by variable
  df_cat_long <- tidyr::pivot_longer(
    df_cat, cols = dplyr::everything(), names_to = "variable", values_to = "category"
  ) |> tidyr::drop_na()
  
  # Calculate the count and proportion of each category within each variable
  df_cat_summary <- df_cat_long |>
    dplyr::group_by(variable, category) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(variable) |>
    dplyr::mutate(
      prop = count / base::sum(count),
      label = base::paste0(base::round(prop * 100, 1), "%\n(n=", count, ")")
    ) |>
    dplyr::ungroup()
  
  # Build the dynamic plot title
  plot_title <- base::paste0(
    "Distribution of categorical variables",
    base::ifelse(dataset_name == "", "", base::paste0(" (", dataset_name, ")"))
  )
  
  # Build the bar chart with one facet per variable
  p <- ggplot2::ggplot(df_cat_summary, ggplot2::aes(x = stats::reorder(category, -count), y = prop, fill = category)) +
    ggplot2::geom_col(color = "black", alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.3, fontface = "bold", size = 3) +
    ggplot2::facet_wrap(~ variable, scales = "free_x", ncol = 3) +
    ggplot2::scale_y_continuous(
      labels = function(x) base::paste0(x * 100, "%"),
      expand = ggplot2::expansion(mult = base::c(0, 0.3))
    ) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Categories",
      y = "Proportion",
      title = plot_title
    ) +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_rect(fill = "#e8f8f5", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      axis.text.x = ggplot2::element_text(face = "bold", size = 10),
      plot.title = ggplot2::element_text(face = "bold", size = 14, margin = ggplot2::margin(b = 15))
    )
  
  base::return(p)
}

# ==============================================================================
# UTILS FOR "5. Missing and implausible values"
# ==============================================================================

#' @title Comparative audit of missing values
#'
#' @description Compares the total number of missing values (NAs) in a dataframe
#' before and after applying an imputation (or censoring) process, and reports
#' the results and the final state of the matrix via the console.
#'
#' @param df_before Dataframe. The original dataset (with NAs).
#' @param df_after Dataframe. The processed dataset (imputed).
#' @param dataset_name String. Descriptive name of the dataset for the report.
#'
#' @return None. Prints the report directly to the console.
audit_missing_values <- function(df_before, df_after, dataset_name) {
  
  # Calculate the total number of NAs in both matrices
  nulls_before <- base::sum(base::is.na(df_before))
  nulls_after <- base::sum(base::is.na(df_after))
  
  # Print the formatted report to the console
  base::cat(base::sprintf("\n--- Missing values audit: %s ---\n", dataset_name))
  base::cat(base::sprintf("Total NAs before imputation: %d\n", nulls_before))
  base::cat(base::sprintf("Total NAs after imputation: %d\n", nulls_after))
  
  # Evaluate the final state of the matrix
  if (nulls_after == 0) {
    base::cat("Status: SUCCESS. The matrix is completely dense and suitable for mathematical modeling.\n\n")
  } else {
    base::cat("Status: WARNING. There are still residual null values in the matrix.\n\n")
  }
  
  base::return(base::invisible())
}

# ==============================================================================
# UTILS FOR "6. Statistical tests"
# ==============================================================================

#' @title Evaluates normality and homoscedasticity assumptions
#'
#' @description Automatically extracts all numeric variables from a dataset and
#' calculates the p-values for the normality test (Kolmogorov-Smirnov with
#' Lilliefors correction) and the homoscedasticity test (Levene, median-centered)
#' for each of them.
#'
#' @param df Dataframe. The training dataset to evaluate.
#' @param target String. Name of the target column (treated as a factor).
#'
#' @return A dataframe summarizing the variables and their respective p-values.
evaluate_initial_assumptions <- function(df, target) {
  
  # Identify the strictly numeric variables
  numeric_vars <- df |>
    dplyr::select(tidyselect::where(base::is.numeric)) |>
    base::names()
  
  # Make sure the target variable is treated as a factor
  target_vec <- base::as.factor(df[[target]])
  
  # Initialize an empty dataframe to store the results
  results <- base::data.frame(
    variable = base::character(),
    p_val_ks_lilliefors = base::numeric(),
    p_val_levene = base::numeric(),
    stringsAsFactors = FALSE
  )
  
  # Iterate over each numeric variable and calculate its p-values
  for (var in numeric_vars) {
    var_values <- df[[var]]
    
    # Kolmogorov-Smirnov test with Lilliefors correction (requires nortest)
    p_ks <- base::tryCatch(
      {
        nortest::lillie.test(stats::na.omit(var_values))$p.value
      }, 
      error = function(e) { 
        NA_real_ 
      }
    )
    
    # Levene's test centered on the median (requires car)
    p_levene <- base::tryCatch(
      {
        car::leveneTest(var_values ~ target_vec, center = stats::median)$`Pr(>F)`[1]
      }, 
      error = function(e) { 
        NA_real_ 
      }
    )
    
    # Append the results of this variable to the results dataframe
    results <- base::rbind(results, base::data.frame(
      variable = var,
      p_val_ks_lilliefors = p_ks,
      p_val_levene = p_levene,
      stringsAsFactors = FALSE
    ))
  }
  
  # Round the p-values to 5 decimal places
  results$p_val_ks_lilliefors <- base::round(results$p_val_ks_lilliefors, 5)
  results$p_val_levene <- base::round(results$p_val_levene, 5)
  
  # Return the final results dataframe
  base::return(base::as.data.frame(results))
}

# ------------------------------------------------------------------------------

#' @title Comparative visualization of Q-Q plots
#'
#' @description Visually compares the fit to a normal distribution of the numeric
#' variables before and after applying a mathematical transformation.
#'
#' @param df_orig Dataframe. The original dataset (raw).
#' @param df_trans Dataframe. The transformed dataset.
#' @param dataset_name String. Dataset name used to customize the plot title.
#' @param color_orig String. HEX code for the original data points.
#' @param color_trans String. HEX code for the transformed data points.
#'
#' @return A ggplot2 object.
plot_qq_comparisons <- function(df_orig, df_trans, dataset_name = "",
                                color_orig = "#e74c3c", color_trans = "#2980b9") {
  
  # Keep only the numeric variables in both dataframes
  df_num_orig <- df_orig |> dplyr::select(tidyselect::where(base::is.numeric))
  df_num_trans <- df_trans |> dplyr::select(tidyselect::where(base::is.numeric))
  
  if (base::ncol(df_num_orig) == 0) {
    base::return(base::message("No numeric variables to visualize."))
  }
  
  # Transform both dataframes to long format, tagging their state
  df_long_orig <- tidyr::pivot_longer(
    df_num_orig, cols = dplyr::everything(), names_to = "variable", values_to = "value"
  ) |>
    tidyr::drop_na() |>
    dplyr::mutate(state = "1. Original")
  
  df_long_trans <- tidyr::pivot_longer(
    df_num_trans, cols = dplyr::everything(), names_to = "variable", values_to = "value"
  ) |>
    tidyr::drop_na() |>
    dplyr::mutate(state = "2. Transformed")
  
  # Combine both long dataframes into a single one for the visualization
  df_final <- base::rbind(df_long_orig, df_long_trans)
  
  # Build the dynamic plot title
  plot_title <- base::paste0(
    "Comparative Q-Q plots: before vs after transformation",
    base::ifelse(dataset_name == "", "", base::paste0(" (", dataset_name, ")"))
  )
  
  # Build the Q-Q plot with facets for each variable and state
  p <- ggplot2::ggplot(df_final, ggplot2::aes(sample = value, color = state)) +
    ggplot2::stat_qq(alpha = 0.5, size = 1) +
    # Theoretical reference line, drawn in dashed black to stand out
    ggplot2::stat_qq_line(color = "black", linewidth = 0.5, linetype = "dashed") +
    ggplot2::facet_wrap(variable ~ state, scales = "free", ncol = 4) +
    ggplot2::scale_color_manual(values = base::c("1. Original" = color_orig, "2. Transformed" = color_trans)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Theoretical quantiles (Normal)",
      y = "Sample quantiles",
      title = plot_title,
      color = "Data state:"
    ) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#ecf0f1", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      axis.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(face = "bold", size = 14, margin = ggplot2::margin(b = 15))
    )
  
  base::return(p)
}

# ------------------------------------------------------------------------------

#' @title Association test for numeric variables (Wilcoxon / Kruskal-Wallis)
#'
#' @description Runs a non-parametric test for each numeric predictor variable
#' against the target variable, automatically adapting to binary or multiclass targets.
#'
#' @param df Dataframe. Transformed training dataset.
#' @param target_col String. Name of the target variable.
#'
#' @return A dataframe with the test results, ordered by ascending p-value.
test_numeric_associations <- function(df, target_col) {
  
  # Extract only the numeric predictors, excluding the target
  num_vars <- base::setdiff(base::names(df)[base::sapply(df, base::is.numeric)], target_col)
  
  # Determine the number of classes in the target, to choose the appropriate test
  target_vector <- base::as.factor(df[[target_col]])
  num_classes <- base::length(base::levels(target_vector))
  
  results <- purrr::map_dfr(num_vars, function(var) {
    
    if (num_classes == 2) {
      # Binary target: Wilcoxon test
      test <- stats::wilcox.test(df[[var]] ~ target_vector, exact = FALSE)
      test_name <- "Wilcoxon (Mann-Whitney)"
      statistic_val <- test$statistic
    } else {
      # Multiclass target (more than 2 levels): Kruskal-Wallis test
      test <- stats::kruskal.test(df[[var]] ~ target_vector)
      test_name <- "Kruskal-Wallis"
      statistic_val <- test$statistic
    }
    
    base::data.frame(
      variable = var,
      test = test_name,
      statistic = statistic_val,
      p_value = test$p.value
    )
  })
  
  # Sort the results by discriminative power (ascending p-value)
  results <- results |> dplyr::arrange(p_value)
  base::return(results)
}

# ------------------------------------------------------------------------------

#' @title Association test for categorical variables (Chi-squared)
#'
#' @description Runs a Pearson's Chi-squared test of independence for each
#' categorical predictor variable against the target variable.
#'
#' @param df Dataframe. Transformed training dataset.
#' @param target_col String. Name of the target variable.
#'
#' @return A dataframe with the test results, ordered by ascending p-value.
test_categorical_associations <- function(df, target_col) {
  
  # Extract the categorical variables (factors or characters), excluding the target
  is_cat <- base::sapply(df, function(x) base::is.factor(x) || base::is.character(x))
  cat_vars <- base::setdiff(base::names(df)[is_cat], target_col)
  
  if (base::length(cat_vars) == 0) {
    base::return(base::message("No categorical predictor variables to evaluate."))
  }
  
  results <- purrr::map_dfr(cat_vars, function(var) {
    # Temporarily suppress warnings caused by low expected counts in the Chi-squared test
    test <- base::suppressWarnings(stats::chisq.test(df[[var]], df[[target_col]]))
    
    base::data.frame(
      variable = var,
      test = "Pearson's Chi-squared",
      statistic_x2 = test$statistic,
      degrees_of_freedom = test$parameter,
      p_value = test$p.value
    )
  })
  
  results <- results |> dplyr::arrange(p_value)
  base::return(results)
}

# ------------------------------------------------------------------------------

#' @title Univariate visualization against the target (numeric predictors)
#'
#' @description Generates a panel of boxplots to visually compare the discriminative
#' capacity of the numeric variables across the classes of the target variable.
#'
#' @param df Dataframe. The data to visualize.
#' @param target_col String. Name of the categorical target variable.
#' @param title String. Main title of the plot.
#'
#' @return A ggplot2 object with faceted panels.
plot_target_relationship_num <- function(df, target_col, title = "Variable distribution by target") {
  
  # Select the target and the numeric variables, and pivot to long format
  df_plot <- df |>
    dplyr::select(dplyr::all_of(target_col), tidyselect::where(base::is.numeric)) |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of(target_col),
      names_to = "variable",
      values_to = "value"
    )
  
  # Build the boxplot grid, one facet per variable
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data[[target_col]], y = value, fill = .data[[target_col]])) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 1.2) +
    ggplot2::facet_wrap(~ variable, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = title, x = "", y = "Value (transformed)")
  
  base::return(p)
}

# ------------------------------------------------------------------------------

#' @title Univariate visualization against the target (categorical predictors)
#'
#' @description Generates a panel of 100% stacked bar charts to visually compare
#' the proportion of predictor categories across the classes of the target variable.
#'
#' @param df Dataframe. The data to visualize.
#' @param target_col String. Name of the categorical target variable.
#' @param title String. Main title of the plot.
#'
#' @return A ggplot2 object with faceted panels, or NULL (invisibly) if no categorical predictors exist.
plot_target_relationship_cat <- function(df, target_col, title = "Category proportion by target") {
  
  # Identify the categorical variables, excluding the target itself
  cat_cols <- df |>
    dplyr::select(tidyselect::where(~ base::is.factor(.x) || base::is.character(.x))) |>
    dplyr::select(-dplyr::all_of(target_col)) |>
    base::colnames()
  
  # Exit silently if there are no categorical predictors to plot (e.g. Breast Cancer)
  if (base::length(cat_cols) == 0) {
    base::message("Note: no extra categorical predictors to plot in this dataset.")
    base::return(base::invisible(NULL))
  }
  
  # Pivot the dataframe to long format for ggplot2
  df_plot <- df |>
    dplyr::select(dplyr::all_of(base::c(target_col, cat_cols))) |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of(target_col),
      names_to = "variable",
      values_to = "category"
    ) |>
    tidyr::drop_na() # Temporarily drop NAs so they do not distort the visual proportions
  
  # Build the stacked proportion bar grid, one facet per variable
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data[[target_col]], fill = category)) +
    ggplot2::geom_bar(position = "fill", color = "white", alpha = 0.85) +
    ggplot2::facet_wrap(~ variable, scales = "free_x") +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    ) +
    ggplot2::labs(title = title, x = "", y = "Proportion", fill = "Category")
  
  base::return(p)
}

# ==============================================================================
# UTILS FOR "8. Outliers"
# ==============================================================================

#' @title Visual audit of outliers via boxplots
#'
#' @description Generates a panel of boxplots for all numeric variables in the
#' dataframe. Red dots mark the outliers detected by the IQR method (Tukey).
#'
#' @param df Dataframe. The data to audit.
#' @param title String. Plot title.
#'
#' @return A ggplot2 object with all the boxplots.
audit_outliers_visually <- function(df, title = "Visual audit of outliers") {
  
  # Select only the numeric variables and pivot to long format
  df_plot <- df |>
    dplyr::select(tidyselect::where(base::is.numeric)) |>
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "variable",
      values_to = "value"
    ) |>
    tidyr::drop_na()
  
  # Build the boxplot panel, one facet per variable
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(y = value)) +
    ggplot2::geom_boxplot(
      fill = "#3498db",
      alpha = 0.5,
      outlier.colour = "red",
      outlier.shape = 16,
      outlier.size = 2
    ) +
    ggplot2::facet_wrap(~ variable, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 8)
    ) +
    ggplot2::labs(title = title, y = "Value (transformed scale)")
  
  base::return(p)
}

# ------------------------------------------------------------------------------

#' @title Audits the truncation impact (bidirectional Tukey method)
#'
#' @description Automatically detects all numeric variables and calculates how many
#' values fall outside the Tukey limits (Q3 + factor * IQR and Q1 - factor * IQR).
#'
#' @param df Dataframe. The data to audit.
#' @param factor Numeric. IQR multiplier (typically 1.5).
#' @param target_col String. Name of the target column (treated as a factor), used to
#' calculate the denominator of the positive class. If NULL, the percentage column is omitted.
#' @param positive_class String. Level of the target considered the positive class
#' (default "M" for Malignant in the cancer dataset).
#'
#' @return A dataframe with the count of detected outliers per variable.
audit_tukey_impact <- function(df, factor = 1.5, target_col = "diagnosis", positive_class = "M") {
  
  numeric_vars <- df |>
    dplyr::select(tidyselect::where(base::is.numeric)) |>
    base::colnames()
  
  # Calculate the total size of the positive class dynamically from the data itself
  if (!base::is.null(target_col) && target_col %in% base::colnames(df)) {
    total_positive_class <- base::sum(df[[target_col]] == positive_class, na.rm = TRUE)
  } else {
    total_positive_class <- NULL
  }
  
  results <- purrr::map_df(numeric_vars, function(var) {
    stats_vec <- stats::quantile(df[[var]], probs = base::c(0.25, 0.75), na.rm = TRUE)
    q1 <- stats_vec[1]
    q3 <- stats_vec[2]
    iqr <- q3 - q1
    
    upper_threshold <- q3 + (factor * iqr)
    lower_threshold <- q1 - (factor * iqr)
    
    # Count the values that exceed either the upper or the lower limit
    n_outliers <- base::sum(df[[var]] > upper_threshold | df[[var]] < lower_threshold, na.rm = TRUE)
    
    row_data <- base::data.frame(
      variable = var,
      n_outliers = n_outliers,
      pct_over_total = base::round((n_outliers / base::nrow(df)) * 100, 2)
    )
    
    # Add the positive class column only when applicable
    if (!base::is.null(total_positive_class)) {
      row_data$pct_over_positive_class <- base::round((n_outliers / total_positive_class) * 100, 2)
    }
    
    row_data
  })
  
  base::return(results)
}

# ==============================================================================
# UTILS FOR "9. Final recipes generation"
# ==============================================================================

# CUSTOM STEP: step_winsorize_iqr()
# Defines a custom tidymodels step that applies IQR-based winsorization.

#' @title IQR-based winsorization step (constructor)
#'
#' @description Defines the step constructor, called from the recipe, which
#' stores the initial parameters supplied by the user.
#'
#' @param recipe A recipe object to which the step will be added.
#' @param ... One or more selector functions to choose which variables are affected.
#' @param role Not used by this step since no new variables are created.
#' @param trained Logical. Has the step already been trained (internal use by recipes).
#' @param factor Numeric. IQR multiplier used to compute the winsorization limits (default 1.5).
#' @param limits List. Winsorization limits per variable, filled in during `prep()`, not set manually.
#' @param skip Logical. Should the step be skipped when baking the training data itself.
#' @param id String. Unique identifier of the step, autogenerated by `recipes::rand_id()`.
#'
#' @return An updated recipe with the new step appended.
#' @export
step_winsorize_iqr <- function(
    recipe, ...,
    role    = NA,
    trained = FALSE,
    factor  = 1.5,
    limits  = NULL,
    skip    = FALSE,
    id      = recipes::rand_id("winsorize_iqr")
) {
  recipes::add_step(
    recipe,
    base::structure(
      base::list(
        terms   = rlang::enquos(...),
        role    = role,
        trained = trained,
        factor  = factor,
        limits  = limits,
        skip    = skip,
        id      = id
      ),
      class = base::c("step_winsorize_iqr", "step")
    )
  )
}

#' @title Prepares the IQR winsorization step
#'
#' @description Defines the `prep()` method, which calculates the winsorization
#' limits exclusively on the training data, so no information from the test set
#' leaks into the transformation.
#'
#' @param x The untrained step, as created by `step_winsorize_iqr()`.
#' @param training Dataframe with the training data used to compute the limits.
#' @param info Dataframe with metadata about the current set of variables (internal use by recipes).
#' @param ... Additional arguments (not used).
#'
#' @return The trained step, with the winsorization limits stored internally.
#' @export
prep.step_winsorize_iqr <- function(x, training, info = NULL, ...) {
  
  # Resolve which variables were selected by the user
  col_names <- recipes::recipes_eval_select(x$terms, training, info)
  
  # Calculate the limits exclusively on the training data (fold or global train)
  calculated_limits <- purrr::map(col_names, function(var) {
    q    <- stats::quantile(training[[var]], probs = base::c(0.25, 0.75), na.rm = TRUE)
    iqr  <- q[2] - q[1]
    base::list(
      inf = q[1] - x$factor * iqr,
      sup = q[2] + x$factor * iqr
    )
  }) |> purrr::set_names(col_names)
  
  # Return the trained step, with the limits now memorized
  base::structure(
    base::list(
      terms   = x$terms,
      role    = x$role,
      trained = TRUE,
      factor  = x$factor,
      limits  = calculated_limits,
      skip    = x$skip,
      id      = x$id
    ),
    class = base::c("step_winsorize_iqr", "step")
  )
}

#' @title Applies the IQR winsorization step
#'
#' @description Defines the `bake()` method, which applies the limits already
#' memorized during `prep()` to any new dataset (train, validation, or test).
#'
#' @param object The trained step, as returned by `prep.step_winsorize_iqr()`.
#' @param new_data Dataframe to which the winsorization limits should be applied.
#' @param ... Additional arguments (not used).
#'
#' @return The dataframe with the selected variables truncated to their memorized limits.
#' @export
bake.step_winsorize_iqr <- function(object, new_data, ...) {
  
  # Apply the limits learned on train to the new dataset
  for (var in base::names(object$limits)) {
    lim <- object$limits[[var]]
    new_data[[var]] <- base::pmin(base::pmax(new_data[[var]], lim$inf), lim$sup)
  }
  
  base::return(new_data)
}

#' @title Print method for the custom winsorization step
#'
#' @description Allows tidymodels to display the step properly inside the
#' recipe summary, instead of falling back to a generic printout.
#'
#' @param x The step object to print.
#' @param width Integer. Maximum width available for the printed selector text.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the step object, after printing its description.
#' @export
print.step_winsorize_iqr <- function(x, width = base::max(20, base::options()$width - 35), ...) {
  base::cat(
    "IQR Winsorization (factor =", x$factor, ") on",
    recipes::format_selectors(x$terms, width = width), "\n"
  )
  base::return(base::invisible(x))
}

# ------------------------------------------------------------------------------

#' @title Verifies Z-score standardization
#'
#' @description Calculates the mean and standard deviation of the continuous
#' numeric variables in a standardized dataframe, ignoring binary variables
#' (dummies). Means should be close to 0 and deviations close to 1. Useful to
#' confirm that `step_normalize()` was applied correctly.
#'
#' @param df Dataframe. Standardized data (result of baking with step_normalize()).
#' @param dataset_name String. Dataset name used in the report header.
#' @param mean_tolerance Numeric. Tolerance threshold for the mean (default 1e-10).
#' @param sd_tolerance Numeric. Tolerance threshold for the standard deviation (default 0.01).
#'
#' @return A dataframe with the mean, standard deviation, and anomaly flags per variable.
verify_standardization <- function(df, dataset_name = "",
                                   mean_tolerance = 1e-10,
                                   sd_tolerance   = 0.01) {
  
  # Select only the numeric columns
  df_num <- df |> dplyr::select(tidyselect::where(base::is.numeric))
  
  # Filter out binary variables (dummies), which should not be evaluated here.
  # Any variable whose only values are 0, 1 (or NA) is considered binary.
  is_continuous <- base::sapply(df_num, function(x) {
    vals <- base::unique(stats::na.omit(x))
    !(base::length(vals) <= 2 && base::all(vals %in% base::c(0, 1)))
  })
  
  df_evaluate <- df_num[, is_continuous, drop = FALSE]
  
  if (base::ncol(df_evaluate) == 0) {
    base::message("No continuous variables to evaluate (all are binary/dummies).")
    base::return(base::invisible(NULL))
  }
  
  summary_stats <- df_evaluate |>
    dplyr::summarise(dplyr::across(
      dplyr::everything(),
      base::list(mean = ~ base::round(base::mean(.x, na.rm = TRUE), 8),
                 sd   = ~ base::round(stats::sd(.x,   na.rm = TRUE), 6))
    )) |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to  = base::c("variable", ".value"),
      names_sep = "_(?=[^_]+$)"
    ) |>
    dplyr::mutate(
      mean_ok = base::abs(mean) < mean_tolerance,
      sd_ok   = base::abs(sd - 1) < sd_tolerance,
      status  = dplyr::case_when(
        mean_ok & sd_ok ~ "OK",
        !mean_ok        ~ "WARNING: mean != 0",
        !sd_ok          ~ "WARNING: sd != 1"
      )
    )
  
  base::return(summary_stats)
}

# ------------------------------------------------------------------------------

#' @title SVD/PCA diagnostics on standardized data
#'
#' @description Applies `prcomp()` on an already centered and scaled dataframe
#' (Z-score), generates the cumulative explained variance plot with a cutoff line,
#' and creates a clean table with the top variables per component and their loadings.
#'
#' @param df Dataframe. Numeric data already standardized (output of baking with step_normalize()).
#' @param dataset_name String. Dataset name used in the plot titles.
#' @param min_variance Numeric. Cumulative variance threshold (default 0.85).
#' @param n_loadings Integer. Number of top variables to display per component (default 5).
#'
#' @return A list with: n_components (integer), cumulative_variance (numeric vector), 
#' pca (the `prcomp` object), plot (ggplot2 object), and top_loadings (wide dataframe).
diagnose_svd <- function(df, dataset_name = "", min_variance = 0.85, n_loadings = 5) {
  
  # Make sure we are only working with the numeric variables
  df_num <- df |> dplyr::select(tidyselect::where(base::is.numeric))
  
  # Run SVD/PCA: center = FALSE and scale. = FALSE since the data is already Z-scored
  pca <- stats::prcomp(df_num, center = FALSE, scale. = FALSE)
  
  # Explained variance per component
  var_exp   <- pca$sdev^2 / base::sum(pca$sdev^2)
  var_accum <- base::cumsum(var_exp)
  
  # Minimum number of components needed to exceed the variance threshold
  n_comp <- base::which(var_accum >= min_variance)[1]
  
  # --- Cumulative explained variance plot ---
  df_var <- base::data.frame(
    component      = base::seq_along(var_accum),
    cum_variance   = var_accum,
    ind_variance   = var_exp
  )
  
  p <- ggplot2::ggplot(df_var, ggplot2::aes(x = component)) +
    ggplot2::geom_col(ggplot2::aes(y = ind_variance), fill = "#3498db", alpha = 0.6, color = "white") +
    ggplot2::geom_line(ggplot2::aes(y = cum_variance), color = "#2c3e50", linewidth = 1.1) +
    ggplot2::geom_point(ggplot2::aes(y = cum_variance), color = "#2c3e50", size = 2.5) +
    ggplot2::geom_hline(yintercept = min_variance, linetype = "dashed", color = "#e74c3c", linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = n_comp, linetype = "dotted", color = "#e74c3c", linewidth = 0.9) +
    ggplot2::annotate(
      "label", x = n_comp + 0.4, y = min_variance - 0.06,
      label = base::sprintf("k = %d\n%.1f%% var.", n_comp, var_accum[n_comp] * 100),
      color = "#e74c3c", size = 3.5, fontface = "bold", fill = "white", label.size = 0.3
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      sec.axis = ggplot2::sec_axis(~ ., labels = scales::percent_format(accuracy = 1), name = "Individual variance (bars)")
    ) +
    ggplot2::scale_x_continuous(breaks = base::seq_len(base::nrow(df_var))) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title    = base::paste0("Explained variance by SVD component - ", dataset_name),
      subtitle = base::sprintf(
        "Red line: %.0f%% threshold | Bars: individual variance | Dark line: cumulative variance", min_variance * 100
      ),
      x = "Principal component", y = "Cumulative variance"
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey40"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # --- Table with the main loadings per component ---
  loadings_mat <- base::as.data.frame(pca$rotation[, 1:n_comp, drop = FALSE])
  loadings_mat$variable <- base::rownames(loadings_mat)
  
  # Build a wide dataframe where each column corresponds to one principal component
  top_loadings <- purrr::map_dfc(1:n_comp, function(k) {
    col_name <- base::paste0("PC", k)
    var_pct <- base::round(var_exp[k] * 100, 1)
    header <- base::sprintf("%s (%.1f%%)", col_name, var_pct)
    
    top_vars <- loadings_mat |>
      dplyr::select(variable, dplyr::all_of(col_name)) |>
      dplyr::rename(loading = dplyr::all_of(col_name)) |>
      dplyr::mutate(abs_loading = base::abs(loading)) |>
      dplyr::arrange(dplyr::desc(abs_loading)) |>
      utils::head(n_loadings) |>
      dplyr::mutate(label = base::sprintf("%s (%.3f)", variable, loading)) |>
      dplyr::pull(label)
    
    res <- base::data.frame(top_vars, stringsAsFactors = FALSE)
    base::colnames(res) <- header
    base::return(res)
  })
  
  base::return(base::list(
    n_components        = n_comp,
    cumulative_variance = var_accum,
    pca                 = pca,
    plot                = p,
    top_loadings        = top_loadings
  ))
}

# ------------------------------------------------------------------------------

#' @title Visualizes a correlation heatmap with coefficients
#'
#' @description Generates a heatmap of the correlations between numeric variables,
#' displaying the numeric coefficient in each cell to make diagnostics easier.
#'
#' @param df Dataframe with the baked data.
#' @param dataset_name String. Dataset name.
#' @param method String. Correlation method ("spearman" or "pearson").
#' @param lab_size Numeric. Font size for the numeric coefficients (default 2.5).
#' @param digits Integer. Number of decimal places for the coefficients (default 2).
#'
#' @return A ggcorrplot object with numeric labels.
plot_correlation_heatmap <- function(df, dataset_name = "", method = "spearman", lab_size = 2.5, digits = 2) {
  
  # Select only the numeric predictors
  df_num <- df |> dplyr::select(tidyselect::where(base::is.numeric))
  
  if (base::ncol(df_num) < 2) {
    base::return(base::message("Not enough numeric variables to correlate."))
  }
  
  # Calculate the correlation matrix
  cor_mat <- stats::cor(df_num, method = method, use = "pairwise.complete.obs")
  
  # Build the dynamic plot title
  plot_title <- base::paste0("Correlation matrix (", base::toupper(method), ") - ", dataset_name)
  
  # Generate the heatmap with ggcorrplot
  p <- ggcorrplot::ggcorrplot(
    cor_mat,
    hc.order = TRUE,           # Hierarchical ordering, to group similar variables together
    type = "lower",            # Show only the lower triangle
    lab = TRUE,                # Enable the numeric coefficient labels
    lab_size = lab_size,       # Adjust the text size to avoid overlapping
    digits = digits,           # Round the decimals to keep the cells clean
    outline.col = "white",
    colors = base::c("#e74c3c", "white", "#3498db"),
    title = plot_title
  ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(face = "bold", size = 12)
    )
  
  base::return(p)
}

# ------------------------------------------------------------------------------

#' @title Visualizes the linearity assumption for all continuous variables (log-odds)
#'
#' @description Extracts all numeric predictors from a dataset, pivots them to
#' long format, calculates independent quantile bins for each of them, and
#' generates a faceted panel with the empirical log-odds (applying the
#' Haldane-Anscombe correction) against the positive target class.
#'
#' @param data Dataframe or tibble with the baked data.
#' @param target_col String. Name of the target variable.
#' @param positive_class String. Exact value of the "positive" class.
#' @param bins Integer. Number of groups into which each predictor is divided (default 10).
#' @param title String. Global title for the plot panel.
#'
#' @return A ggplot2 object with the faceted panel.
plot_empirical_logodds_multi <- function(data, target_col, positive_class, bins = 10, title = "Linearity audit (empirical log-odds)") {
  
  # Isolate the target and all continuous numeric predictors, ignoring dummies when possible.
  # Since dummies only take values 0 and 1, cut_number() would throw an error on them,
  # so variables with only two unique values are filtered out below.
  df_numeric <- data |>
    dplyr::select(tidyselect::where(base::is.numeric), dplyr::all_of(target_col)) |>
    dplyr::filter(!base::is.na(.data[[target_col]])) |>
    dplyr::mutate(
      .binary_target = base::ifelse(base::as.character(.data[[target_col]]) == positive_class, 1, 0)
    ) |>
    dplyr::select(-dplyr::all_of(target_col))
  
  # Filter out dummy variables (those with at most two unique values)
  cols_to_keep <- base::names(df_numeric)[base::sapply(df_numeric, function(x) base::length(base::unique(x)) > 2)]
  cols_to_keep <- base::c(cols_to_keep, ".binary_target")
  df_numeric <- df_numeric[, cols_to_keep]
  
  # Pivot the predictors to long format
  df_long <- df_numeric |>
    tidyr::pivot_longer(
      cols = -c(.binary_target),
      names_to = "predictor",
      values_to = "value"
    ) |>
    dplyr::filter(!base::is.na(value))
  
  # Calculate the grouped quantiles and the empirical log-odds within each bin
  df_summary <- df_long |>
    dplyr::group_by(predictor) |>
    dplyr::mutate(.bin = ggplot2::cut_number(value, n = bins)) |>
    dplyr::group_by(predictor, .bin) |>
    dplyr::summarise(
      .x_mean = base::mean(value, na.rm = TRUE),
      .events = base::sum(.binary_target, na.rm = TRUE),
      .non_events = dplyr::n() - .events,
      .log_odds = base::log((.events + 0.5) / (.non_events + 0.5)),
      .groups = "drop"
    )
  
  # Build the faceted panel, one facet per predictor
  p <- ggplot2::ggplot(df_summary, ggplot2::aes(x = .x_mean, y = .log_odds)) +
    ggplot2::geom_point(size = 1.5, color = "#2c3e50", alpha = 0.8) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "#e74c3c", se = FALSE, linetype = "dashed", linewidth = 0.8) +
    ggplot2::facet_wrap(~ predictor, scales = "free_x") +
    ggplot2::labs(
      title = title,
      x = "Mean quantile value",
      y = "Empirical log-odds (Haldane-Anscombe)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      strip.text = ggplot2::element_text(face = "bold", size = 8, color = "#34495e"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7)
    )
  
  base::return(p)
}