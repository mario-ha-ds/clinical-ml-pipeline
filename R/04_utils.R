# This script contains the auxiliary functions necessary to evaluate all
# project models (unsupervised and supervised) against the test set.

# ==============================================================================
# UTILS FOR "2.1. Centroid-based models: K-Means and K-Medians"
# ==============================================================================

#' @title Cluster labeling by majority vote (using exclusively Train)
#'
#' @description Translates the unnamed geometric labels of a clustering model
#' (Cluster 1, Cluster 2...) into the actual clinical classes of the target,
#' based exclusively on the most frequent actual class among the Train
#' patients that fell into each cluster. This mapping is calculated only once
#' on the Train set and remains fixed to be applied later on the Test set,
#' avoiding any information leakage from the Test set.
#'
#' @param cluster_labels_train Vector with the cluster labels assigned by the
#' model on the Train set (e.g. model$cluster or flexclust::clusters(model)).
#' @param target_train Factor vector with the actual target of those same Train patients.
#'
#' @return A dataframe with two columns: cluster (the numeric cluster
#' identifier) and translated_class (the actual clinical label assigned by majority).
create_cluster_class_mapping <- function(cluster_labels_train, target_train) {
  
  cross_table <- base::table(cluster_labels_train, target_train)
  
  # Identify the most frequent actual class (column) for each cluster (row)
  mapping <- base::data.frame(
    cluster          = base::rownames(cross_table),
    translated_class = base::colnames(cross_table)[base::apply(cross_table, 1, base::which.max)],
    stringsAsFactors = FALSE
  )
  
  base::return(mapping)
}

# ------------------------------------------------------------------------------

#' @title New patient assignment to the nearest centroid
#'
#' @description Calculates, for each patient in a new dataset (typically the
#' Test set), the distance to each of the already fixed centroids from a
#' trained clustering model, and assigns the patient to the nearest one. The
#' model is never retrained, nor are the centroids recalculated: they remain
#' exactly as they were fixed during training.
#'
#' @param new_data Dataframe or purely numeric matrix (same variables and
#' order as those used to train the original model).
#' @param centers Matrix of already fixed centroids (rows = clusters, columns = variables).
#' @param method String. Distance to use: "euclidean" (K-Means) or "manhattan" (K-Medians).
#'
#' @return An integer vector with the identifier of the cluster assigned to each row of new_data.
assign_nearest_cluster <- function(new_data, centers, method = "euclidean") {
  
  new_data_mat <- base::as.matrix(new_data)
  
  # Distance matrix: each row is a patient, each column a centroid
  dist_to_centers <- base::matrix(
    NA_real_,
    nrow = base::nrow(new_data_mat),
    ncol = base::nrow(centers)
  )
  
  for (k in 1:base::nrow(centers)) {
    diffs <- base::scale(new_data_mat, center = centers[k, ], scale = FALSE)
    
    if (method == "euclidean") {
      dist_to_centers[, k] <- base::sqrt(base::rowSums(diffs^2))
    } else if (method == "manhattan") {
      dist_to_centers[, k] <- base::rowSums(base::abs(diffs))
    } else {
      base::stop("The method must be 'euclidean' or 'manhattan'.")
    }
  }
  
  # Assign each patient to the column index (centroid) with the minimum distance
  assigned_cluster <- base::apply(dist_to_centers, 1, base::which.min)
  
  base::return(assigned_cluster)
}

# ------------------------------------------------------------------------------

#' @title Full external evaluation of a centroid-based model
#'
#' @description Encapsulates the entire external evaluation flow for a
#' clustering model (K-Means or K-Medians): builds the cluster-to-class
#' mapping on Train, assigns Test patients to the nearest centroid, translates
#' those assignments into clinical classes, and calculates the hard metrics
#' (Accuracy, Sensitivity, F1-Score) along with the confusion matrix.
#'
#' @param train_data Train dataframe, already projected into the PCA space (numeric predictors only).
#' @param train_target Factor vector with the actual Train target.
#' @param train_cluster_labels Vector with the cluster labels assigned on Train by the model.
#' @param test_data Test dataframe, already projected into the same PCA space (numeric predictors only).
#' @param test_target Factor vector with the actual Test target.
#' @param centers Matrix of already fixed model centroids.
#' @param method String. Distance to use: "euclidean" (K-Means) or "manhattan" (K-Medians).
#'
#' @return A list with three elements: predictions (dataframe with truth,
#' assigned cluster, and translated predicted class), mapping (the
#' cluster-to-class mapping used), and metrics (a yardstick tibble).
evaluate_centroid_model <- function(train_data, train_target, train_cluster_labels,
                                    test_data, test_target, centers, method = "euclidean") {
  
  # Cluster-to-class mapping, calculated exclusively on Train
  mapping <- create_cluster_class_mapping(train_cluster_labels, train_target)
  
  # Assignment of Test patients to the nearest centroid (without retraining)
  test_cluster <- assign_nearest_cluster(test_data, centers, method = method)
  
  # Translation of the Test assignments to clinical classes, using the fixed Train mapping
  predicted_class <- mapping$translated_class[base::match(test_cluster, mapping$cluster)]
  
  # Assemble the predictions dataframe, aligning the factor levels with the actual target
  predictions <- base::data.frame(
    truth            = test_target,
    assigned_cluster = test_cluster,
    .pred_class      = base::factor(predicted_class, levels = base::levels(test_target))
  )
  
  # Calculate the hard metrics (Accuracy, Sensitivity, F1-Score)
  metrics <- get_hard_metrics()(
    data     = predictions,
    truth    = truth,
    estimate = .pred_class
  )
  
  base::return(base::list(
    predictions = predictions,
    mapping     = mapping,
    metrics     = metrics
  ))
}

# ------------------------------------------------------------------------------

#' @title Binary collapse of a multiclass target
#'
#' @description Collapses a target with three or more classes into two
#' categories, grouping every class other than the reference class into a
#' single "rest" category. Used to allow a direct comparison between models
#' with a different number of clusters (e.g. K-Means k=2 vs K-Medians k=3 in
#' the CVD scenario).
#'
#' @param target_vector Factor vector with the original target (3 or more classes).
#' @param positive_class String with the exact name of the class to keep isolated.
#' @param rest_label String with the name of the collapsed category (default "OTHER").
#'
#' @return A 2-level factor: positive_class and rest_label.
collapse_binary_target <- function(target_vector, positive_class, rest_label = "OTHER") {
  
  collapsed <- base::ifelse(
    base::as.character(target_vector) == positive_class,
    positive_class,
    rest_label
  )
  
  base::return(base::factor(collapsed, levels = base::c(positive_class, rest_label)))
}

# ------------------------------------------------------------------------------

#' @title PCA visualization of the external evaluation of a centroid-based model
#'
#' @description Projects the Test set patients onto their first two principal
#' components (PC1 and PC2), coloring each point according to the geometric
#' cluster assigned by the model and using a different shape according to its
#' actual clinical class. This visually assesses how well the geometric
#' partition (color) matches the actual clinical label (shape): a perfect
#' match would mean each color is always associated with the same shape.
#' A high-contrast, colorblind-safe palette (deuteranopia/protanopia) is used
#' instead of the default ggplot2 palette.
#'
#' @param test_data Test dataframe, already projected into the PCA space (numeric predictors only).
#' @param assigned_cluster Vector with the cluster assigned to each test patient.
#' @param actual_class Factor vector with the actual clinical class of each test patient.
#' @param title String. Plot title.
#'
#' @return A ggplot2 object with the PC1 vs PC2 projection.
plot_pca_evaluation <- function(test_data, assigned_cluster, actual_class, title = "External evaluation (PC1 vs PC2)") {
  
  # High-contrast, colorblind-safe palette (based on Okabe-Ito / Wong)
  colorblind_palette <- base::c(
    "#000000", # black
    "#F0E442", # yellow
    "#0072B2", # blue
    "#D55E00", # orange/brick red
    "#009E73", # teal/blue-green
    "#CC79A7"  # pink/magenta
  )
  
  df_plot <- base::data.frame(
    PC1          = test_data[[1]],
    PC2          = test_data[[2]],
    cluster      = base::factor(assigned_cluster),
    actual_class = actual_class
  )
  
  n_clusters <- base::length(base::unique(df_plot$cluster))
  
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = PC1, y = PC2, color = cluster, shape = actual_class)) +
    ggplot2::geom_point(size = 2.2, alpha = 0.75, stroke = 1) +
    ggplot2::scale_color_manual(values = colorblind_palette[1:n_clusters]) +
    ggplot2::scale_shape_manual(values = base::c(16, 17, 15, 3, 8, 4)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      legend.position = "bottom"
    ) +
    ggplot2::labs(
      title = title,
      x     = "Principal component 1 (PC1)",
      y     = "Principal component 2 (PC2)",
      color = "Assigned cluster",
      shape = "Actual class"
    )
  
  base::return(p)
}

# ------------------------------------------------------------------------------

#' @title Binary external evaluation of a centroid-based model
#'
#' @description Evaluates a centroid-based model against a multiclass target
#' collapsed to a binary scheme (e.g. K-Means k=2 on the CVD scenario).
#' The target is collapsed to binary prior to building the majority-vote
#' mapping on Train, ensuring a direct 2-class alignment on Test without
#' intermediate redundant mappings.
#'
#' @param train_data Train dataframe, already projected into the PCA space.
#' @param train_target Factor vector with the actual Train target (original scale).
#' @param train_cluster_labels Vector with the cluster labels assigned on Train by the model.
#' @param test_data Test dataframe, already projected into the same PCA space.
#' @param test_target Factor vector with the actual Test target (original scale).
#' @param centers Matrix of already fixed model centroids.
#' @param method String. Distance to use: "euclidean" (K-Means) or "manhattan" (K-Medians).
#' @param positive_class String with the exact name of the class to keep isolated (e.g. "HIGH").
#' @param rest_label String with the name of the collapsed category (default "OTHER").
#'
#' @return A list with three elements: predictions, mapping, and metrics.
evaluate_centroid_model_binary <- function(train_data, train_target, train_cluster_labels,
                                           test_data, test_target, centers, method = "euclidean",
                                           positive_class, rest_label = "OTHER") {
  
  # Collapse both the Train and Test targets directly into the binary scheme
  binary_train_target <- collapse_binary_target(train_target, positive_class, rest_label)
  binary_test_target  <- collapse_binary_target(test_target, positive_class, rest_label)
  
  # Cluster-to-class mapping calculated directly over the binary labels on Train
  mapping <- create_cluster_class_mapping(train_cluster_labels, binary_train_target)
  
  # Assignment of Test patients to the nearest centroid (without retraining)
  test_cluster <- assign_nearest_cluster(test_data, centers, method = method)
  
  # Translation of Test assignments to binary classes using the fixed Train mapping
  predicted_class <- mapping$translated_class[base::match(test_cluster, mapping$cluster)]
  
  # Assemble the predictions dataframe
  predictions <- base::data.frame(
    truth            = binary_test_target,
    assigned_cluster = test_cluster,
    .pred_class      = base::factor(predicted_class, levels = base::levels(binary_test_target))
  )
  
  # Calculate the hard metrics (Accuracy, Sensitivity, F1-Score)
  metrics <- get_hard_metrics()(
    data        = predictions,
    truth       = truth,
    estimate    = .pred_class,
    event_level = "first"
  )
  
  base::return(base::list(
    predictions = predictions,
    mapping     = mapping,
    metrics     = metrics
  ))
}

# ==============================================================================
# UTILS FOR "2.2. Penalized logistic regression"
# ==============================================================================

#' @title External evaluation of a binary model with a configurable threshold
#'
#' @description Applies an already trained binary classification workflow on
#' a raw Test set, converts the continuous probability of the positive class
#' to a hard class using the specified decision threshold, and calculates the
#' full set of clinical metrics (Accuracy, Sensitivity, F1-Score, ROC-AUC, PR-AUC).
#'
#' @param workflow_fit Already trained workflow (parsnip/workflows).
#' @param test_data_raw Raw Test dataframe (unbaked).
#' @param target_col String with the name of the actual target column.
#' @param positive_class String with the exact name of the positive class.
#' @param threshold Numeric. Decision threshold (default 0.5).
#'
#' @return A list with two elements: predictions (dataframe) and metrics (tibble).
evaluate_logistic_binary <- function(workflow_fit, test_data_raw, target_col,
                                     positive_class, threshold = 0.5) {
  
  target_vec <- test_data_raw[[target_col]]
  col_prob_positive <- base::paste0(".pred_", positive_class)
  negative_class <- base::setdiff(base::levels(target_vec), positive_class)
  
  # Detect whether the positive class is the first or second level of the target
  event_level <- base::ifelse(base::levels(target_vec)[1] == positive_class, "first", "second")
  
  # --- Tidymodels bypass ---
  prepped_recipe <- workflows::extract_recipe(workflow_fit, estimated = TRUE)
  parsnip_model  <- workflows::extract_fit_parsnip(workflow_fit)
  
  test_baked <- recipes::bake(prepped_recipe, new_data = test_data_raw)
  
  # Continuous probabilities from the native engine
  probs <- stats::predict(parsnip_model, new_data = test_baked, type = "prob")
  
  # Assemble the predictions, preserving the original factor levels
  predictions <- probs |>
    dplyr::bind_cols(truth = target_vec) |>
    dplyr::mutate(
      .pred_class = dplyr::if_else(.data[[col_prob_positive]] > threshold, positive_class, negative_class),
      .pred_class = base::factor(.pred_class, levels = base::levels(target_vec))
    )
  
  # Calculate the full metric set, using the native probability column
  metrics <- get_clinical_metrics()(
    data        = predictions,
    truth       = truth,
    estimate    = .pred_class,
    !!rlang::sym(col_prob_positive),
    event_level = event_level
  )
  
  base::return(base::list(predictions = predictions, metrics = metrics))
}

# ------------------------------------------------------------------------------

#' @title External evaluation of a multiclass model (multinomial logistic regression)
#'
#' @description Applies an already trained multiclass classification workflow
#' on a raw Test set. Manually extracts the recipe and the engine to predict
#' without hiding the target. Extracts the continuous probabilities of all
#' classes and assigns the class with the highest probability.
#'
#' @param workflow_fit Already trained workflow (parsnip/workflows).
#' @param test_data_raw Raw Test dataframe (unbaked).
#' @param target_col String with the name of the actual target column.
#'
#' @return A list with two elements: predictions (dataframe with truth,
#' probabilities, and predicted class) and metrics.
evaluate_logistic_multiclass <- function(workflow_fit, test_data_raw, target_col) {
  
  target_vec <- test_data_raw[[target_col]]
  
  # --- Tidymodels bypass ---
  # Uses the same workaround strategy as in the binary scenario
  prepped_recipe <- workflows::extract_recipe(workflow_fit, estimated = TRUE)
  parsnip_model  <- workflows::extract_fit_parsnip(workflow_fit)
  
  test_baked <- recipes::bake(prepped_recipe, new_data = test_data_raw)
  
  # Continuous probabilities for every class, straight from the native engine
  probs <- stats::predict(parsnip_model, new_data = test_baked, type = "prob")
  
  # Hard class, obtained directly from predict(type = "class") on the native engine
  classes <- stats::predict(parsnip_model, new_data = test_baked, type = "class")
  
  predictions <- probs |>
    dplyr::bind_cols(classes) |>
    dplyr::bind_cols(truth = target_vec) |>
    dplyr::mutate(.pred_class = base::factor(.pred_class, levels = base::levels(target_vec)))
  
  # Probability column names, passed dynamically to yardstick
  prob_cols <- base::paste0(".pred_", base::levels(target_vec))
  
  # Calculate the full set of clinical metrics (one-vs-rest, macro-average)
  metrics <- get_clinical_metrics()(
    data     = predictions,
    truth    = truth,
    estimate = .pred_class,
    !!!rlang::syms(prob_cols)
  )
  
  base::return(base::list(predictions = predictions, metrics = metrics))
}

# ==============================================================================
# UTILS FOR "2.3. C5.0 with cost matrix"
# ==============================================================================

#' @title External evaluation of a C5.0 model with a cost matrix (binary)
#'
#' @description Applies an already trained C5.0 workflow (with an integrated
#' asymmetric cost matrix) on a raw Test set, and calculates only the hard
#' metric set (Accuracy, Sensitivity, F1-Score). ROC-AUC and PR-AUC are not
#' calculated, since the cost matrix distorts the tree's internal probabilities
#' and invalidates any analysis based on a continuous spectrum of thresholds.
#' Manually extracts the trained recipe and the C5.0 engine from the workflow
#' (instead of calling predict(workflow_fit, ...) directly), since the recipe
#' includes a step_select(all_outcomes()) that requires seeing the target
#' column in new_data when baking. predict.workflow() strips the target before
#' passing the data to bake(), which causes an error. Baking manually on
#' test_data_raw (which still has the target) and then predicting with the
#' extracted parsnip engine avoids that conflict.
#'
#' @param workflow_fit Already trained workflow (parsnip/workflows), with an integrated recipe and model.
#' @param test_data_raw Raw Test dataframe (unbaked), with the target column present.
#' @param target_col String with the name of the actual target column.
#' @param positive_class String with the exact name of the positive class (e.g. "M"),
#' used exclusively to set the correct event_level in the hard metrics.
#'
#' @return A list with two elements: predictions (dataframe with truth and predicted class) and metrics.
evaluate_c50_binary <- function(workflow_fit, test_data_raw, target_col, positive_class) {
  
  target_vec  <- test_data_raw[[target_col]]
  event_level <- base::ifelse(base::levels(target_vec)[1] == positive_class, "first", "second")
  
  # --- Tidymodels bypass ---
  # Bake manually, preserving the target, then predict with the extracted parsnip engine
  prepped_recipe <- workflows::extract_recipe(workflow_fit, estimated = TRUE)
  parsnip_model  <- workflows::extract_fit_parsnip(workflow_fit)
  
  test_baked <- recipes::bake(prepped_recipe, new_data = test_data_raw)
  
  classes <- stats::predict(parsnip_model, new_data = test_baked, type = "class")
  
  predictions <- classes |>
    dplyr::bind_cols(truth = target_vec) |>
    dplyr::mutate(.pred_class = base::factor(.pred_class, levels = base::levels(target_vec)))
  
  metrics <- get_hard_metrics()(
    data        = predictions,
    truth       = truth,
    estimate    = .pred_class,
    event_level = event_level
  )
  
  base::return(base::list(predictions = predictions, metrics = metrics))
}

# ------------------------------------------------------------------------------

#' @title External evaluation of a C5.0 model with a cost matrix (multiclass)
#'
#' @description Multiclass variant of evaluate_c50_binary(), designed for the
#' CVD scenario. Calculates only the hard metrics, without ROC-AUC or PR-AUC.
#' Uses the same manual bypass (extract_recipe + bake + predict on the
#' parsnip engine) to avoid the conflict between step_select(all_outcomes())
#' and predict.workflow().
#'
#' @param workflow_fit Already trained workflow (parsnip/workflows), with an integrated recipe and model.
#' @param test_data_raw Raw Test dataframe (unbaked), with the target column present.
#' @param target_col String with the name of the actual target column.
#'
#' @return A list with two elements: predictions (dataframe with truth and predicted class) and metrics.
evaluate_c50_multiclass <- function(workflow_fit, test_data_raw, target_col) {
  
  target_vec <- test_data_raw[[target_col]]
  
  # --- Tidymodels bypass ---
  prepped_recipe <- workflows::extract_recipe(workflow_fit, estimated = TRUE)
  parsnip_model  <- workflows::extract_fit_parsnip(workflow_fit)
  
  test_baked <- recipes::bake(prepped_recipe, new_data = test_data_raw)
  
  classes <- stats::predict(parsnip_model, new_data = test_baked, type = "class")
  
  predictions <- classes |>
    dplyr::bind_cols(truth = target_vec) |>
    dplyr::mutate(.pred_class = base::factor(.pred_class, levels = base::levels(target_vec)))
  
  metrics <- get_hard_metrics()(
    data     = predictions,
    truth    = truth,
    estimate = .pred_class
  )
  
  base::return(base::list(predictions = predictions, metrics = metrics))
}

# ==============================================================================
# UTILS FOR "2.4. Random Forest"
# ==============================================================================

#' @title External evaluation of a Random Forest model (binary)
#'
#' @description Applies an already trained Random Forest workflow on a raw
#' Test set. Manually extracts the recipe and the predictive engine to avoid
#' the step_select(all_outcomes()) conflict. Converts the continuous
#' probabilities to a hard class using the specified threshold, and calculates
#' the full set of metrics (Accuracy, Sensitivity, F1, ROC-AUC, PR-AUC), since
#' Random Forest does not rely on a cost matrix.
#'
#' @param workflow_fit Already trained workflow (parsnip/workflows).
#' @param test_data_raw Raw Test dataframe (unbaked).
#' @param target_col String with the name of the actual target column.
#' @param positive_class String with the exact name of the positive class (e.g. "M").
#' @param threshold Numeric. Decision threshold (default 0.5).
#'
#' @return A list with two elements: predictions (dataframe) and metrics (tibble).
evaluate_rf_binary <- function(workflow_fit, test_data_raw, target_col,
                               positive_class, threshold = 0.5) {
  
  target_vec <- test_data_raw[[target_col]]
  col_prob_positive <- base::paste0(".pred_", positive_class)
  negative_class <- base::setdiff(base::levels(target_vec), positive_class)
  
  event_level <- base::ifelse(base::levels(target_vec)[1] == positive_class, "first", "second")
  
  # --- Tidymodels bypass ---
  prepped_recipe <- workflows::extract_recipe(workflow_fit, estimated = TRUE)
  parsnip_model  <- workflows::extract_fit_parsnip(workflow_fit)
  
  test_baked <- recipes::bake(prepped_recipe, new_data = test_data_raw)
  
  # Continuous probabilities from the native engine
  probs <- stats::predict(parsnip_model, new_data = test_baked, type = "prob")
  
  # Assemble the hard predictions, preserving the original factor levels
  predictions <- probs |>
    dplyr::bind_cols(truth = target_vec) |>
    dplyr::mutate(
      .pred_class = dplyr::if_else(.data[[col_prob_positive]] > threshold, positive_class, negative_class),
      .pred_class = base::factor(.pred_class, levels = base::levels(target_vec))
    )
  
  # Calculate the full set of clinical metrics
  metrics <- get_clinical_metrics()(
    data        = predictions,
    truth       = truth,
    estimate    = .pred_class,
    !!rlang::sym(col_prob_positive),
    event_level = event_level
  )
  
  base::return(base::list(predictions = predictions, metrics = metrics))
}

# ------------------------------------------------------------------------------

#' @title External evaluation of a Random Forest model (multiclass)
#'
#' @description Multiclass variant of evaluate_rf_binary() for the CVD
#' scenario. Uses the tidymodels bypass and calculates the full set of
#' clinical metrics (Accuracy, macro Sensitivity, macro F1, macro ROC-AUC,
#' and macro PR-AUC).
#'
#' @param workflow_fit Already trained workflow (parsnip/workflows).
#' @param test_data_raw Raw Test dataframe (unbaked).
#' @param target_col String with the name of the actual target column.
#'
#' @return A list with two elements: predictions (dataframe) and metrics (tibble).
evaluate_rf_multiclass <- function(workflow_fit, test_data_raw, target_col) {
  
  target_vec <- test_data_raw[[target_col]]
  
  # --- Tidymodels bypass ---
  prepped_recipe <- workflows::extract_recipe(workflow_fit, estimated = TRUE)
  parsnip_model  <- workflows::extract_fit_parsnip(workflow_fit)
  
  test_baked <- recipes::bake(prepped_recipe, new_data = test_data_raw)
  
  # Probabilities and hard classes, straight from the native engine
  probs   <- stats::predict(parsnip_model, new_data = test_baked, type = "prob")
  classes <- stats::predict(parsnip_model, new_data = test_baked, type = "class")
  
  predictions <- probs |>
    dplyr::bind_cols(classes) |>
    dplyr::bind_cols(truth = target_vec) |>
    dplyr::mutate(.pred_class = base::factor(.pred_class, levels = base::levels(target_vec)))
  
  prob_cols <- base::paste0(".pred_", base::levels(target_vec))
  
  # Calculate the full set of clinical metrics
  metrics <- get_clinical_metrics()(
    data     = predictions,
    truth    = truth,
    estimate = .pred_class,
    !!!rlang::syms(prob_cols)
  )
  
  base::return(base::list(predictions = predictions, metrics = metrics))
}