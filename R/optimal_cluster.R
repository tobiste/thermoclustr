#' Locate elbow in elbow curve of wss (within-cluster sum of squares)
#'
#' @param wss_values numeric vector of wss values.
#'
#' @returns integer index of the elbow point in the wss curve.
#' @noRd
.find_elbow <- function(wss_values) {
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



#' Convenience function to find average silhouette width in cluster
#'
#' @param d distance matrix
#' @param integer vector with *k* different integer cluster codes
#'
#' @returns numeric. The average silhouette width
#' @noRd
#'
#' @importFrom cluster silhouette
#'
#' @examples
#' # example data
#' data(tT_paths1)
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.4) |> path_diss()
#' .get_ave_sil_width(tT_paths_subset$diss, cluster = kmeans(dmat, 3)$cluster)
#' .get_ave_sil_width(tT_paths_subset$diss, cluster = path_hcut(dmat, 3)$cluster)
.get_ave_sil_width <- function(d, cluster) {
  ss <- cluster::silhouette(x = cluster, dist = d)
  mean(ss[, 3])
}


#' Optimal number of clusters
#'
#' Determines the optimal number of clusters in a
#' dissimilarity matrix using *average silhouette width*.
#'
#' @param x either an object of class `"HeFTy"` (output of [read_hefty()]),
#' `"tTdiss` (output of [path_diss()]), or
#' a `data.frame` containing the `time`, `temperature`, and `segment` columns
#' of the modeled paths. If `x` is not a `"tTdiss"` object, it will be
#' calculated first by [path_diss()] using its default settings.
#' @param FUNcluster cluster function. Default is [path_hcut()].
#' @param k.max integer. the maximum number of clusters to consider, must be at least two.
#' @param n.threshold integer. If the number of paths is greater than this value,
#' a random sample of size `n.threshold` of paths will be selected from the
#' dissimilarity matrix to determine the optimal number of clusters.
#' Must be greater than `k.max`. Default is `Inf`, i.e., no sampling.
#' Consider to specify this parameter if processing takes too much time.
#' @param linecolor line color
#' @param ... optionally further arguments for `FUNcluster`
#'
#' @returns list. `$optimal` (integer) is the optimal number
#' of clusters, `$plot` is a ggplot showing the average silhouette width for each number of clusters.
#' @export
#'
#' @importFrom ggplot2 ggplot geom_line geom_point aes labs geom_vline
#'
#' @examples
#' # example data
#' data(tT_paths1)
#' set.seed(20250411)
#' res1 <- path_nbclust(tT_paths1, n.threshold = 100)
#' res1
#'
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.4)
#' res2 <- path_nbclust(tT_paths_subset)
#' res2
path_nbclust <- function(x,
                         FUNcluster = path_hcut,
                         k.max = 10,
                         n.threshold = Inf,
                         linecolor = "#B63679FF",
                         ...) {
  if (inherits(x, "HeFTy") || inherits(x, "data.frame")) {
    if (inherits(x, "data.frame")) {
      stopifnot(all(c("time", "temperature", "segment") %in% colnames(x)))
    }
    x <- path_diss(x)
  } else if (!inherits(x, "tTdiss")) {
    stop("Input `x` must be of class 'HeFTy', 'tTdiss', or 'data.frame' with columns 'time', 'temperature', and 'segment'.")
  }

  y <- clusters <- cluster <- NULL

  diss <- x$diss

  n_paths <- nrow(diss)
  stopifnot(n_paths > k.max)
  stopifnot(k.max > 2)
  stopifnot(is.numeric(n.threshold))

  # cluster only a subset:
  if (n_paths > n.threshold) {
    stopifnot(n.threshold > k.max)
    diss_mat <- as.matrix(diss)
    rnd <- sample(1:n_paths, size = n.threshold)
    diss_mat <- diss_mat[rnd, rnd]
    diss <- as.dist(diss_mat)
  }

  v <- sapply(2:k.max, function(i) {
    clust <- FUNcluster(diss, i, ...)
    .get_ave_sil_width(diss, clust$cluster)
  })
  v <- c(0, v)

  optimal_nbc <- which.max(v)

  silh_plot <- data.frame(
    clusters = as.factor(1:k.max),
    y = v,
    stringsAsFactors = TRUE
  ) |>
    ggplot(aes(x = clusters, y = y)) +
    geom_line(group = 1) +
    geom_point() +
    labs(
      title = "Optimal number of clusters",
      x = "Number of clusters (k)",
      y = "Overall average silhouette width"
    ) +
    geom_vline(xintercept = which.max(v), linetype = 2, color = linecolor)

  return(list(optimal = optimal_nbc, plot = silh_plot))
  invisible(silh_plot)
}
