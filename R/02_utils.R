# This script contains the auxiliary functions necessary for the
# unsupervised modeling section of the project.

# ==============================================================================
# UTILS FOR "1. Centroid-based models: K-Means and K-Medians"
# ==============================================================================

#' @title Evaluates metrics for hyperparameter k (clustering)
#'
#' @description Iterates over a range of k values, fitting a K-Means or
#' K-Medians model for each of them. Calculates the cost function (WCSS for L2,
#' WCAE for L1) and the mean Silhouette coefficient, strictly respecting the
#' topological metric of the chosen algorithm (Euclidean or Manhattan), and
#' generates the corresponding evaluation plots.
#'
#' @param data Dataframe or purely numeric matrix (without the target variable).
#' @param k_range Numeric vector with the k values to explore (e.g. 2:10).
#' @param algorithm String. Algorithm to use: "kmeans" or "kmedians".
#' @param dataset_name String. Dataset name used in the plot title.
#'
#' @return A list with two elements: metrics (a dataframe with the columns k,
#' cost, and silhouette) and plot (a ggplot2 object ready to be rendered).
evaluate_k_search <- function(data, k_range = 2:10, algorithm = "kmeans", dataset_name = "Dataset") {
  
  # Validate that the algorithm argument is one of the supported options
  if (!algorithm %in% base::c("kmeans", "kmedians")) {
    base::stop("The algorithm must be 'kmeans' or 'kmedians'.")
  }
  
  # Initialize vectors to store the results of every k value
  cost_vals <- base::numeric(base::length(k_range))
  sil_vals  <- base::numeric(base::length(k_range))
  
  # Pre-calculate the distance matrix used later for the Silhouette coefficient
  if (algorithm == "kmeans") {
    dist_matrix <- stats::dist(data, method = "euclidean")
  } else {
    dist_matrix <- stats::dist(data, method = "manhattan")
  }
  
  # Iterate over every candidate value of k
  for (i in base::seq_along(k_range)) {
    k <- k_range[i]
    
    if (algorithm == "kmeans") {
      # Fit K-Means (L2), which seeks to minimize the WCSS
      base::set.seed(42)
      model <- stats::kmeans(data, centers = k, nstart = 25)
      
      cost_vals[i]   <- model$tot.withinss
      cluster_labels <- model$cluster
      
    } else if (algorithm == "kmedians") {
      # Fit K-Medians (L1) via kcca, which seeks to minimize the WCAE
      base::set.seed(42)
      model <- flexclust::kcca(base::as.matrix(data), k = k,
                               family = flexclust::kccaFamily("kmedians"),
                               control = base::list(initcent = "kmeanspp"))
      
      # Calculate the L1 cost function (WCAE): cluster size times average distance to the median
      cost_vals[i]   <- base::sum(model@clusinfo$size * model@clusinfo$av_dist)
      cluster_labels <- flexclust::clusters(model)
    }
    
    # Calculate the mean Silhouette coefficient, only when more than one cluster exists
    if (base::length(base::unique(cluster_labels)) > 1) {
      sil_obj <- cluster::silhouette(cluster_labels, dist_matrix)
      sil_vals[i] <- base::mean(sil_obj[, "sil_width"])
    } else {
      sil_vals[i] <- 0
    }
  }
  
  # Structure the results into a single dataframe
  results_df <- base::data.frame(
    k = k_range,
    cost = cost_vals,
    silhouette = sil_vals
  )
  
  # --- Build the evaluation plot ---
  
  # Pivot to long format so both metrics can be faceted with facet_wrap()
  df_long <- tidyr::pivot_longer(results_df,
                                 cols = base::c("cost", "silhouette"),
                                 names_to = "metric",
                                 values_to = "value")
  
  # Refactor the metric names into readable facet titles, keeping the insertion order
  df_long$metric <- base::factor(df_long$metric,
                                 levels = base::c("cost", "silhouette"),
                                 labels = base::c("Cost function (Elbow method)", "Mean Silhouette coefficient"))
  
  # Build the dynamic plot title
  plot_title <- base::paste("Optimal k search:", base::toupper(algorithm), "-", dataset_name)
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = k, y = value, color = metric)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_x_continuous(breaks = k_range) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5)
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Number of clusters (k)",
      y = "Evaluation value"
    )
  
  base::return(base::list(metrics = results_df, plot = p))
}

# ------------------------------------------------------------------------------

#' @title Trains candidate centroid-based models
#'
#' @description Iteratively trains K-Means and K-Medians models for a vector
#' of specific k values, so they can later be compared and evaluated.
#'
#' @param data Dataframe or purely numeric matrix (without the target variable).
#' @param k_candidates Numeric vector with the k values to train (e.g. c(2, 3, 4)).
#' @param algorithm String. Algorithm to use: "kmeans" or "kmedians".
#'
#' @return A named list with the trained models, accessible by their corresponding k.
train_centroid_candidates <- function(data, k_candidates = base::c(2, 3, 4), algorithm = "kmeans") {
  
  # Validate that the algorithm argument is one of the supported options
  if (!algorithm %in% base::c("kmeans", "kmedians")) {
    base::stop("Invalid algorithm.")
  }
  
  models_list <- base::list()
  
  # Train one model per candidate value of k
  for (k in k_candidates) {
    k_name <- base::paste0("k_", k)
    
    if (algorithm == "kmeans") {
      base::set.seed(42)
      models_list[[k_name]] <- stats::kmeans(data, centers = k, nstart = 25)
      
    } else {
      base::set.seed(42)
      models_list[[k_name]] <- flexclust::kcca(base::as.matrix(data), k = k,
                                               family = flexclust::kccaFamily("kmedians"),
                                               control = base::list(initcent = "kmeanspp"))
    }
  }
  
  base::return(models_list)
}

# ------------------------------------------------------------------------------

#' @title Visual PCA projection of clusters
#'
#' @description Extracts the first two principal components (PC1 and PC2) from
#' the already preprocessed dataset, to visualize the resulting clusters in 2D.
#'
#' @param data Dataframe or training matrix, already transformed by SVD/PCA.
#' @param models_list List of models, as returned by train_centroid_candidates().
#' @param algorithm String. Algorithm used to train the models: "kmeans" or "kmedians".
#' @param dataset_name String. Dataset name used in the plot title.
#'
#' @return A ggplot2 object with one facet per k value.
plot_candidates_pca <- function(data, models_list, algorithm, dataset_name) {
  
  # Extract the first two dimensions (the data is already expressed in the PCA space)
  df_pca <- base::as.data.frame(data[, 1:2])
  base::colnames(df_pca) <- base::c("PC1", "PC2")
  
  # Extract the cluster labels from every model and consolidate them into a long dataframe
  plot_data <- base::data.frame()
  
  for (k_name in base::names(models_list)) {
    model <- models_list[[k_name]]
    
    # Extract the labels depending on which package produced the model
    if (algorithm == "kmeans") {
      labels <- base::as.factor(model$cluster)
    } else {
      labels <- base::as.factor(flexclust::clusters(model))
    }
    
    # Build a temporary dataframe for this specific k
    temp_df <- df_pca
    temp_df$cluster <- labels
    temp_df$k_value <- base::toupper(base::gsub("_", "=", k_name)) # e.g. "k_2" becomes "K=2"
    
    plot_data <- base::rbind(plot_data, temp_df)
  }
  
  # Build the dynamic plot title
  plot_title <- base::paste("Visual projection (PC1 vs PC2):", base::toupper(algorithm), "-", dataset_name)
  
  # Map both color and shape to the cluster variable, one facet per k value
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC1, y = PC2, color = cluster, shape = cluster)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::facet_wrap(~ k_value, ncol = base::length(models_list)) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "bottom"
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Principal component 1 (PC1)",
      y = "Principal component 2 (PC2)",
      color = "Cluster",
      shape = "Cluster"
    )
  
  base::return(p)
}

# ------------------------------------------------------------------------------

#' @title Evaluates internal metrics for centroid-based models
#'
#' @description Calculates the Silhouette, Davies-Bouldin, and Calinski-Harabasz
#' indices for a list of candidate models. The metrics are computed via distance
#' matrices, strictly respecting the topological metric of the algorithm: L2
#' (Euclidean) for K-Means and L1 (Manhattan) for K-Medians. For K-Medians, the
#' centroid is manually recalculated as the exact median per coordinate, to
#' avoid relying on flexclust's internal approximation.
#'
#' @param data Dataframe or purely numeric training matrix.
#' @param models_list List of models, as returned by train_centroid_candidates().
#' @param algorithm String. Algorithm used to train the models: "kmeans" or "kmedians".
#'
#' @return A dataframe with the metrics summary for each k value.
evaluate_centroid_models <- function(data, models_list, algorithm = "kmeans") {
  
  # Validate that the algorithm argument is one of the supported options
  if (!algorithm %in% base::c("kmeans", "kmedians")) {
    base::stop("Invalid algorithm.")
  }
  
  # Prepare the global variables shared across all candidate models
  n_rows <- base::nrow(data)
  dist_method <- base::ifelse(algorithm == "kmeans", "euclidean", "manhattan")
  dist_mat <- stats::dist(data, method = dist_method)
  
  results_df <- base::data.frame()
  
  for (k_name in base::names(models_list)) {
    model <- models_list[[k_name]]
    
    # Extract the cluster assignments and centers depending on the package used
    if (algorithm == "kmeans") {
      clusters <- model$cluster
      centers  <- model$centers
      k        <- base::nrow(centers)
      global_center <- base::colMeans(data) # L2 center of gravity (mean)
    } else {
      clusters <- flexclust::clusters(model)
      k        <- base::length(base::unique(clusters))
      
      # Recalculate the centroid as the exact median per coordinate, instead of
      # relying on flexclust's internal iterative approximation
      centers <- base::matrix(0, nrow = k, ncol = base::ncol(data))
      for (i in 1:k) {
        cluster_pts <- base::as.matrix(data[clusters == i, , drop = FALSE])
        centers[i, ] <- base::apply(cluster_pts, 2, stats::median)
      }
      
      global_center <- base::apply(data, 2, stats::median) # L1 center of gravity (median)
    }
    
    # --- Silhouette ---
    sil <- base::mean(cluster::silhouette(clusters, dist_mat)[, "sil_width"])
    
    # Initialize the variables needed for the manual Davies-Bouldin and Calinski-Harabasz calculation
    s_dispersion <- base::numeric(k) # Mean intra-cluster dispersion
    w_dispersion <- 0                # Total intra-cluster dispersion
    b_dispersion <- 0                # Total inter-cluster dispersion
    
    # --- Dispersions (L1 or L2, depending on the algorithm) ---
    for (i in 1:k) {
      cluster_pts <- base::as.matrix(data[clusters == i, , drop = FALSE])
      n_i <- base::nrow(cluster_pts)
      
      if (n_i == 0) next
      
      # Calculate the distance from every point in the cluster to its centroid
      diffs <- base::scale(cluster_pts, center = centers[i, ], scale = FALSE)
      
      if (algorithm == "kmeans") {
        dists_to_center <- base::sqrt(base::rowSums(diffs^2)) # Euclidean distance
        s_dispersion[i] <- base::mean(dists_to_center)        # Used for Davies-Bouldin
        w_dispersion    <- w_dispersion + base::sum(dists_to_center^2) # Used for Calinski-Harabasz (sum of squares)
        b_dispersion    <- b_dispersion + n_i * base::sum((centers[i, ] - global_center)^2)
      } else {
        dists_to_center <- base::rowSums(base::abs(diffs))    # Manhattan distance
        s_dispersion[i] <- base::mean(dists_to_center)        # Used for Davies-Bouldin
        w_dispersion    <- w_dispersion + base::sum(dists_to_center)   # Used for Calinski-Harabasz (sum of absolutes)
        b_dispersion    <- b_dispersion + n_i * base::sum(base::abs(centers[i, ] - global_center))
      }
    }
    
    # --- Davies-Bouldin index ---
    db_index <- 0
    for (i in 1:k) {
      max_r <- 0
      for (j in 1:k) {
        if (i != j) {
          # Calculate the distance between centroids i and j
          if (algorithm == "kmeans") {
            d_ij <- base::sqrt(base::sum((centers[i, ] - centers[j, ])^2))
          } else {
            d_ij <- base::sum(base::abs(centers[i, ] - centers[j, ]))
          }
          
          # Avoid a division by zero if two centroids collapse onto the same point
          if (d_ij == 0) d_ij <- 1e-10
          
          r_ij <- (s_dispersion[i] + s_dispersion[j]) / d_ij
          if (r_ij > max_r) max_r <- r_ij
        }
      }
      db_index <- db_index + max_r
    }
    db_index <- db_index / k
    
    # --- Calinski-Harabasz index ---
    # Apply mathematical protections in case the within-cluster dispersion collapses to 0
    # (which can happen with perfectly separated clusters on duplicated data)
    if (w_dispersion == 0) w_dispersion <- 1e-10
    ch_index <- (b_dispersion / (k - 1)) / (w_dispersion / (n_rows - k))
    
    # Append the results of this model to the results dataframe
    results_df <- base::rbind(results_df, base::data.frame(
      k_candidate       = k,
      silhouette        = sil,
      davies_bouldin    = db_index,
      calinski_harabasz = ch_index
    ))
  }
  
  base::return(results_df)
}

# ==============================================================================
# UTILS FOR "2. Density-based models: DBSCAN and OPTICS"
# ==============================================================================

#' @title Generates OPTICS reachability plots
#'
#' @description Trains OPTICS models for a list of MinPts values, using an
#' infinite Epsilon to map the full density structure of the space. Generates
#' a topographic plot (reachability plot) to visually compare how the
#' definition of "valleys" (potential clusters) changes across different MinPts.
#'
#' @param data Dataframe or purely numeric training matrix.
#' @param minpts_values Numeric vector with the MinPts values to evaluate.
#' @param dataset_name String. Dataset name used in the plot title.
#'
#' @return A ggplot2 object with one facet per MinPts value.
plot_optics_reachability <- function(data, minpts_values, dataset_name = "Dataset") {
  
  plot_data <- base::data.frame()
  
  for (m_val in minpts_values) {
    
    # Train OPTICS: eps = Inf allows mapping the entire space in a single run.
    # Note: this is computationally heavy for n > 100000, but works fine for our
    # biomedical dataset sizes.
    opt_model <- dbscan::optics(data, minPts = m_val, eps = Inf)
    
    # Extract the reachability distance, following the OPTICS ordering
    reach_dist <- opt_model$reachdist[opt_model$order]
    
    # Impute the first point (infinite reachability by definition) with the
    # maximum finite value plus a small margin, purely for plotting purposes
    max_reach <- base::max(reach_dist[!base::is.infinite(reach_dist)], na.rm = TRUE)
    reach_dist[base::is.infinite(reach_dist)] <- max_reach * 1.05
    reach_dist[base::is.na(reach_dist)] <- max_reach * 1.05 # Extra protection
    
    # Assemble a temporary dataframe for this MinPts value
    temp_df <- base::data.frame(
      order_idx    = 1:base::length(reach_dist),
      reachability = reach_dist,
      minpts_label = base::paste("minPts =", m_val)
    )
    
    plot_data <- base::rbind(plot_data, temp_df)
  }
  
  # Refactor to keep the facets in the original insertion order
  plot_data$minpts_label <- base::factor(plot_data$minpts_label,
                                         levels = base::paste("minPts =", minpts_values))
  
  # Build the reachability (horizon-style) plot
  plot_title <- base::paste("OPTICS topographic analysis (reachability plot) -", dataset_name)
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = order_idx, y = reachability)) +
    # geom_segment() simulates the classic dense OPTICS bar plot
    ggplot2::geom_segment(ggplot2::aes(xend = order_idx, yend = 0), color = "darkblue", alpha = 0.7) +
    ggplot2::facet_wrap(~ minpts_label, ncol = 1, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Topological ordering of points (OPTICS sequence)",
      y = "Reachability distance"
    )
  
  base::return(p)
}

# ------------------------------------------------------------------------------

#' @title Extracts and visualizes DBSCAN clusters from OPTICS
#'
#' @description Performs horizontal cuts (Epsilon) on a base OPTICS model to
#' simulate DBSCAN extractions. Calculates the resulting cluster distribution
#' and the noise rate (cluster 0) for each cut, and returns a unified
#' reachability plot colored by extraction value.
#'
#' @param optics_model Object generated by dbscan::optics().
#' @param eps_values Numeric vector with the Epsilon values to test.
#' @param dataset_name String. Dataset name used in the plot title (e.g. "Breast Cancer").
#'
#' @return A list with three elements: summary (a dataframe with the extraction
#' statistics), plot (a ggplot2 object with the visual panel), and models (a
#' list with the extraction objects returned by extractDBSCAN()).
extract_and_plot_dbscan <- function(optics_model, eps_values, dataset_name = "Dataset") {
  
  # Prepare the base topological data space shared across all extractions
  n_points <- base::length(optics_model$order)
  reach_dist <- optics_model$reachdist[optics_model$order]
  
  # Impute the first point, which has no prior reachability by definition
  max_r <- base::max(reach_dist[!base::is.infinite(reach_dist)], na.rm = TRUE)
  reach_dist[base::is.infinite(reach_dist) | base::is.na(reach_dist)] <- max_r * 1.05
  
  # Initialize the containers used to accumulate results across cuts
  plot_data <- base::data.frame()
  summary_list <- base::list()
  extracted_models <- base::list()
  
  for (eps in eps_values) {
    
    # Perform the simulated DBSCAN extraction at this Epsilon value
    ext_model <- dbscan::extractDBSCAN(optics_model, eps_cl = eps)
    extracted_models[[base::paste0("eps_", eps)]] <- ext_model
    
    # Audit the resulting clusters and the amount of noise.
    # extractDBSCAN() returns cluster labels in the original input order, so we
    # reorder them here since the plot follows the topological OPTICS order.
    labels_ordered <- base::as.factor(ext_model$cluster[optics_model$order])
    
    counts <- base::table(ext_model$cluster)
    n_noise <- if ("0" %in% base::names(counts)) counts["0"] else 0
    pct_noise <- base::round((base::as.integer(n_noise) / n_points) * 100, 2)
    n_clusters <- base::length(counts[base::names(counts) != "0"])
    
    # Save the summary statistics for this cut
    summary_list[[base::paste0("eps_", eps)]] <- base::data.frame(
      epsilon          = eps,
      num_clusters     = n_clusters,
      noise_patients   = base::as.integer(n_noise),
      noise_percentage = pct_noise
    )
    
    # Prepare the plotting data for this cut
    temp_df <- base::data.frame(
      order_idx    = 1:n_points,
      reachability = reach_dist,
      cluster      = labels_ordered,
      cut_label    = base::paste0("Epsilon = ", eps, " | Noise = ", pct_noise, "%")
    )
    
    plot_data <- base::rbind(plot_data, temp_df)
  }
  
  # Consolidate the summary dataframe across all cuts
  summary_df <- base::do.call(base::rbind, summary_list)
  base::rownames(summary_df) <- NULL
  
  # Adjust the plot levels to respect the original Epsilon insertion order
  cut_levels <- base::unique(plot_data$cut_label)
  plot_data$cut_label <- base::factor(plot_data$cut_label, levels = cut_levels)
  
  # Build a custom color palette: cluster "0" (noise) is always dark gray
  my_colors <- base::c(
    "0" = "#555555",
    "1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A",
    "4" = "#984EA3", "5" = "#FF7F00", "6" = "#FFFF33"
  )
  
  # Build the continuous reachability bar plot
  plot_title <- base::paste("Dynamic DBSCAN extraction -", dataset_name)
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = order_idx, y = reachability, color = cluster)) +
    ggplot2::geom_segment(ggplot2::aes(xend = order_idx, yend = 0), alpha = 0.8) +
    ggplot2::scale_color_manual(values = my_colors) +
    # Add the theoretical Epsilon cut as a horizontal reference line
    ggplot2::geom_hline(
      data = base::data.frame(cut_label = cut_levels, eps = eps_values),
      ggplot2::aes(yintercept = eps),
      color = "red", linetype = "dashed", linewidth = 0.8
    ) +
    ggplot2::facet_wrap(~ cut_label, ncol = 1, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      legend.position = "bottom"
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Point order (OPTICS)",
      y = "Reachability distance",
      color = "Assigned cluster (0 = noise)"
    )
  
  base::return(base::list(summary = summary_df, plot = p, models = extracted_models))
}