# This script contains the auxiliary functions necessary for the
# supervised modeling section of the project.

# ==============================================================================
# UTILS FOR "1. Evaluation strategy"
# ==============================================================================

#' @title Builds the clinical metric set
#'
#' @description Defines and returns a `metric_set` from the `yardstick` package,
#' configured with the five core metrics used to evaluate the predictive models
#' in this project: Accuracy, Sensitivity (Recall), F1-Score, ROC-AUC, and
#' PR-AUC. Internally supports both binary classification (Breast Cancer) and
#' multiclass classification (Cardiovascular Risk), automatically applying a
#' macro-average in multiclass contexts.
#'
#' @return A `metric_set` object containing the evaluation functions, ready to
#' be plugged into `tune_grid()` or `fit_resamples()`. Note that the predictions
#' must contain both the hard classes (`.pred_class`) and the continuous
#' probabilities (`.pred_X`), so the areas under the curve can be calculated.
get_clinical_metrics <- function() {
  metrics <- yardstick::metric_set(
    yardstick::accuracy,
    yardstick::sens,
    yardstick::f_meas,
    yardstick::roc_auc,
    yardstick::pr_auc
  )
  
  base::return(metrics)
}

# ------------------------------------------------------------------------------

#' @title Builds the clinical metric set (hard classification)
#'
#' @description Defines and returns a `metric_set` from the `yardstick` package
#' restricted to metrics that only require the hard predicted class, for models
#' that do not return class probabilities (e.g. C5.0 with a cost matrix).
#'
#' @return A `metric_set` object containing Accuracy, Sensitivity, and F1-Score,
#' ready to be plugged into `tune_grid()` or `fit_resamples()`.
get_hard_metrics <- function() {
  metrics <- yardstick::metric_set(
    yardstick::accuracy,
    yardstick::sens,
    yardstick::f_meas
  )
  
  base::return(metrics)
}

# ==============================================================================
# UTILS FOR "2, 3, and 4. Model tuning and evaluation"
# ==============================================================================

#' @title Hierarchical selection of hyperparameters
#'
#' @description Extracts the best combination of hyperparameters from a set of
#' tuning results, based on a strict clinical hierarchy: Sensitivity > PR-AUC >
#' F1-Score > ROC-AUC > Accuracy. Ties in the first metric are broken by the
#' next one in the hierarchy, and so on.
#'
#' @param tune_results Object with tuning results, as returned by `tune::tune_grid()`.
#'
#' @return A one-row dataframe with the winning combination of hyperparameters
#' and its associated metrics.
get_best_hierarchical <- function(tune_results) {
  # Note: id_cols requires tidyselect's own c() to combine .config (a bare
  # column name) with any_of(...). Using base::c() here would force eager
  # evaluation and R would look for a *variable* called .config instead of
  # a *column* called .config, causing an "object not found" error.
  best_params <- tune_results |>
    tune::collect_metrics() |>
    # Pivot dynamically, keeping whichever hyperparameters are actually present
    tidyr::pivot_wider(
      id_cols = c(.config, tidyselect::any_of(base::c("penalty", "mixture", "min_n", "trees", "mtry"))),
      names_from = .metric,
      values_from = mean
    ) |>
    dplyr::arrange(
      dplyr::desc(sens),
      dplyr::desc(pr_auc),
      dplyr::desc(f_meas),
      dplyr::desc(roc_auc),
      dplyr::desc(accuracy)
    ) |>
    dplyr::slice(1)
  
  base::return(best_params)
}

# ------------------------------------------------------------------------------

#' @title Hierarchical selection of hyperparameters (hard classification)
#'
#' @description Extracts the best combination of hyperparameters from a set of
#' tuning results, based on a strict clinical hierarchy restricted to the
#' metrics available for hard classification models: Sensitivity > F1-Score >
#' Accuracy.
#'
#' @param tune_results Object with tuning results, as returned by `tune::tune_grid()`.
#'
#' @return A one-row dataframe with the winning combination of hyperparameters
#' and its associated metrics.
get_best_hard_hierarchical <- function(tune_results) {
  # Note: see get_best_hierarchical() above for why id_cols must use
  # tidyselect's own c() instead of base::c() to combine .config with any_of(...).
  best_params <- tune_results |>
    tune::collect_metrics() |>
    tidyr::pivot_wider(
      id_cols = c(.config, tidyselect::any_of(base::c("penalty", "mixture", "min_n", "trees", "mtry"))),
      names_from = .metric,
      values_from = mean
    ) |>
    dplyr::arrange(
      dplyr::desc(sens),
      dplyr::desc(f_meas),
      dplyr::desc(accuracy)
    ) |>
    dplyr::slice(1)
  
  base::return(best_params)
}

# ------------------------------------------------------------------------------

#' @title Generates a visual confusion matrix
#'
#' @description Takes a dataframe of predictions, calculates the confusion
#' matrix, and generates a professional heatmap-style plot using ggplot2.
#'
#' @param data Dataframe containing the truth column and the `.pred_class` column.
#' @param truth Unquoted name of the column with the actual target variable.
#' @param title String. Title to display on the plot.
#' @param subtitle String. Subtitle to display (defaults to internal training context).
#'
#' @return A ggplot2 object ready to be rendered.
plot_confusion_matrix <- function(data, truth, title = "Confusion Matrix", subtitle = "Internal evaluation on the training set") {
  truth_col <- rlang::enquo(truth)
  
  p <- data |>
    # Synchronize the predicted factor levels with the actual truth levels
    dplyr::mutate(.pred_class = base::factor(.pred_class, levels = base::levels(!!truth_col))) |>
    yardstick::conf_mat(truth = !!truth_col, estimate = .pred_class) |>
    ggplot2::autoplot(type = "heatmap") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Clinical truth",
      y = "Model prediction"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  
  base::return(p)
}

# ------------------------------------------------------------------------------

# Explicit parsnip engine dependency declarations for renv detection
# nocov start
if (FALSE) {
  library(C50)
  library(ranger)
}
# nocov end

# ==============================================================================
# UTILS FOR "4. Random Forest"
# ==============================================================================

#' @title Plots Random Forest variable importance (Gini)
#'
#' @description Extracts the native engine from a trained Random Forest model,
#' pulls the full variable importance ranking as a table, and generates a
#' custom bar chart. Guarantees that all requested variables are displayed.
#'
#' @param rf_fit A finalized and trained Random Forest model (workflow object).
#' @param num_features Integer. Number of predictor variables to display.
#' @param plot_title String. Main title of the plot.
#' @param plot_subtitle String. Descriptive subtitle of the plot.
#' @param fill_color String. HEX code or color name used for the bars.
#'
#' @return A ggplot2 object ready to be rendered.
plot_rf_importance <- function(rf_fit, num_features, plot_title, plot_subtitle, fill_color = "#2c3e50") {
  
  # Extract the native engine (ranger) from the trained workflow
  engine_fit <- workflows::extract_fit_engine(rf_fit)
  
  # Extract the raw importance values as a tibble instead of relying on an auto-plot
  importance_data <- vip::vi(engine_fit) |>
    utils::head(num_features) |>
    # Reorder the factor levels so the bar chart is sorted from highest to lowest
    dplyr::mutate(Variable = stats::reorder(Variable, Importance))
  
  # Build the bar chart manually with ggplot2, for full control over the styling
  p <- ggplot2::ggplot(importance_data, ggplot2::aes(x = Importance, y = Variable)) +
    ggplot2::geom_col(fill = fill_color, color = "black", alpha = 0.8) +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Importance (cumulative reduction in Gini impurity)",
      y = "Predictors"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text.y = ggplot2::element_text(face = "bold", size = 9) # Slightly smaller font, so every label fits
    )
  
  base::return(p)
}