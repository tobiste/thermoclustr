#' Locate elbow in elbow curve of wss (within-cluster sum of squares)
#'
#' @param wss_values numeric vector of wss values.
#'
#' @returns integer index of the elbow point in the wss curve.
#' @noRd
find_elbow <- function(wss_values) {
  n_points <- length(wss_values)
  all_coords <- cbind(1:n_points, wss_values)

  # Line from first to last point
  line_vec <- all_coords[n_points, ] - all_coords[1, ]
  line_vec_norm <- line_vec / sqrt(sum(line_vec^2))

  # Compute the distance from each point to the line
  vec_from_first <- sweep(all_coords, 2, all_coords[1, ])
  scalar_product <- rowSums(vec_from_first * matrix(rep(line_vec_norm, n_points), ncol = 2, byrow = TRUE))
  proj <- outer(scalar_product, line_vec_norm)
  vec_to_line <- vec_from_first - proj
  dist_to_line <- sqrt(rowSums(vec_to_line^2))

  which.max(dist_to_line)
}


#' Optimal number of clusters
#'
#' Convenience function to determine the optimal number of clusters in a
#' dissimilarity matrix using the three methods
#' *within cluster sums of squares*, *average silhouette*, and *gap statistic*.
#'
#' @param x either an object of class `"HeFTy"` (output of [read_hefty()]),
#' `"tTdiss` (output of [path_diss()]), or
#' a `data.frame` containing the `time`, `temperature`, and `segment` columns of the modeled paths.
#' @param FUNcluster cluster function
#' @param wss,silhouette,gap.statistic logical. Select if this statistic should be included.
#' @param k.max integer. the maximum number of clusters to consider, must be at least two.
# #' @param nboot integer. integer, number of Monte Carlo ("bootstrap") samples.
# #' Used only for determining the number of clusters using gap statistic.
#' @param dim integer. Number of dimensions for [stats::cmdscale()] used for the gap statistic, default is 4.
#' @param n.threshold integer. If the number of paths is greater than this value,
#' a random sample of size `n.threshold` of paths will be used to determine the optimal number of clusters.
#' Must be greater than `k.max`. Default is `Inf`, i.e., no sampling.
#' Consider to specify this parameter if processing takes too much time.
#' @param linecolor line color
#' @param ... optionally further arguments for [factoextra::fviz_nbclust()]
#'
#' @returns list containing the results of the three methods, the optimal number
#' of clusters (median of the results), and the plots for each method.
#' @export
#'
#' @importFrom factoextra fviz_nbclust
#' @importFrom ggplot2 labs geom_vline
#'
#' @examples
#' # example data
#' data(tT_paths1)
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.4)
#' res <- path_nbclust(tT_paths_subset)
#' res$optimal
path_nbclust <- function(x, 
                         FUNcluster = factoextra::hcut,
                         k.max = 10,
                         wss = FALSE, silhouette = TRUE, gap.statistic = FALSE,
                         # nboot = 50,
                         dim = 4,
                         n.threshold = Inf,
                         linecolor = "#1D1147FF",
                         ...) {
  if (inherits(x, "HeFTy")) {
    x <- path_diss(x) # calculate dissimilarities
  } else if (inherits(x, "tTdiss")) {
  } else if (inherits(x, "data.frame")) {
    stopifnot(c("time", "temperature", "segment") %in% colnames(x))
    x <- path_diss(x)
  }

  median <- temp_q <- time_min <- time_median <- time_max <- temp_sd <- temp_IQR <- temp_median <- temp_5 <- temp_95 <- temp_max <- temp_min <- NULL

  diss <- x$diss
  diss_mat <- as.matrix(diss)

  n_paths <- nrow(diss)
  stopifnot(n_paths > k.max)
  stopifnot(is.numeric(n.threshold))
  if (n_paths > n.threshold) {
    stopifnot(n.threshold > k.max)
    rnd <- sample(1:n_paths, size = n.threshold)
    diss <- diss_mat <- diss_mat[rnd, rnd]
  }

  #### within cluster sums of squares
  if(wss){
  plot_wss <- factoextra::fviz_nbclust(diss_mat,
    FUNcluster = FUNcluster,
    method = "wss",
    k.max = k.max,
    verbose = FALSE,
    # print.summary = TRUE,
    linecolor = linecolor,
    ...
  ) +
    ggplot2::labs(title = NULL)

  nb1 <- find_elbow(plot_wss$data$y)

  plot_wss <- plot_wss +
    ggplot2::geom_vline(xintercept = nb1, lty = 2, color = linecolor)
  } else {
    plot_wss <- nb1 <- NULL
  }
  

  #### Silhouette method
  if(silhouette){
  plot_silh <- factoextra::fviz_nbclust(diss_mat,
    FUNcluster = FUNcluster,
    method = "silhouette",
    k.max = k.max,
    verbose = FALSE,
    # print.summary = FALSE,
    linecolor = linecolor,
    ...
  ) +
    ggplot2::labs(title = NULL)

  nb2 <- plot_silh$layers[[3]]$data |>
    unname() |>
    as.numeric()
  } else {
    plot_silh <- nb2 <- NULL
  }
  

  #### Gap statistic method
  if(gap.statistic){
  plot_gapst <- factoextra::fviz_nbclust(stats::cmdscale(diss, k = dim),
    FUNcluster = FUNcluster,
    method = "gap_stat",
    # nboot = nboot,
    k.max = k.max,
    verbose = FALSE,
    # print.summary = FALSE,
    linecolor = linecolor,
    ...
  ) +
    ggplot2::labs(title = NULL)

  nb3 <- plot_gapst$layers[[4]]$data |>
    unname() |>
    as.numeric()
  } else {
    plot_gapst <- snb3 <- NULL
  }


  ### Combine results
  results <- c("wss" = nb1, "silhouette" = nb2, "gap_stat" = nb3)
  opt_nb <- median(results)
  list(results = results, optimal = opt_nb, wss = plot_wss, silhouette = plot_silh, gap_stat = plot_gapst)
}
