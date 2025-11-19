#' Clustering thermal histories models
#'
#' Groups t-T paths into "path families" based on the *Hausdorff* or
#' *Fr&#233;chet distance* between paths.
#'
#' @param x either an object of class `"HeFTy"` (output of [read_hefty()]), an
#' object of `"tTdiss"` (output of [path_diss()]),
#' or a `data.frame` containing the `time`, `temperature`, and `segment` columns of the
#' modeled paths.
#' @param k an integer scalar or vector with the desired number of groups.
#' Ignored when `dist` equal to `dbscan` or `hdbscan` (**see Note**).
#' @param dist character. Algorithm to calculate a dissimilarity matrix
#' (distance) for lines; either `Hausdorff` (the default), or
#' `Frechet`.
#' @param method character. Clustering method to be applied. Currently implemented are
#' \describe{
#'  \item{`"hclust"`}{Hierarchical Clustering using [stats::hclust()], the default)}
#'  \item{`"kmeans"`}{K-Means Clustering using [stats::kmeans()])}
#'  \item{`"pam"`}{Partitioning Around Medoids using [cluster::pam()]}
#'  \item{`"dbscan"`}{Density-based Spatial Clustering of Applications with
#'  Noise (DBSCAN) using [dbscan::dbscan()] (**see Note**)}
#'  \item{`"hdbscan"`}{Hierarchical DBSCAN using [dbscan::hdbscan()]
#'  (**see Note**)}
#'  \item{`"specc"`}{Spectral Clustering using [kernlab::specc()]}
#'  \item{`"agnes"`}{Agglomerative hierarchical clustering using [cluster::agnes()]}
#'  \item{`"diana"`}{Divisive hierarchical clustering using [cluster::diana()]}
#'  \item{`"clara"`}{Clustering Large Applications using [cluster::clara()]}
#'  \item{`"fanny"`}{Fuzzy Analysis Clustering using [cluster::fanny()]}
#' }
#' @param naming character. Naming scheme for clusters. One of
#' `"asis"` (output of underlying cluster algorithm),
#' `"GOF"` (ranks of the mean GOF values within clusters), and
#' `"size"` (ranks of the size of clusters).
#' Outliers detected by `dbscan` of `hdbscan` will always named as `0`.
#' @param warn logical. Should there be a warning message if at least one cluster
#' contains less than (`threshold` * 100)% of the total paths?
#' @param threshold numeric. The significance threshold as a fraction of the
#' total amount of paths (`0.01` by default). If one cluster contains less
#' paths per total paths than this value, a warning message is created. Ignored
#' if `warn` is `FALSE`.
#' @param ... additional arguments passed to cluster method.
#'
#' @note that `dbscan` and `hdbscan` methods require `eps` and `minPts` arguments.
#' Optimal `eps` values can be visually estimated from the "knee" in a k-nearest
#' neighbor distance plot using [dbscan::kNNdistplot()]'.
#'
#' @details If you want to use a different clustering method that is not built
#' in the current version of `thermochron`, you can use the distance matrix
#' produced by [path_diss()] and feed your cluster algorithm.
#'
#' @return a data.frame with the path `segment` (integer) and `cluster` (factor)
#' @export
#'
#' @importFrom sf st_as_sf st_distance
#' @importFrom dplyr summarise group_by left_join right_join join_by count transmute min_rank row_number mutate
#' @importFrom stats hclust cutree as.dist kmeans
#' @importFrom cluster pam agnes diana clara fanny
#' @importFrom kernlab specc
#' @importFrom dbscan dbscan hdbscan
#'
#' @seealso [path_diss()] for calculating dissimilarities and [path_nbclust()]
#' for determining the optimal number of clusters.
#'
#' @examples
#' data(tT_paths1)
#' tT_paths1$paths <- subset(tT_paths1$paths, Comp_GOF >= 0.4)
#'
#' # cluster the paths
#' cluster_paths(tT_paths1, k = 3)
cluster_paths <- function(
    x, k,
    dist = c("Hausdorff", "Frechet"),
    method = c("hclust", "kmeans", "pam", "dbscan", "hdbscan", "specc", "diana", "agnes", "clara", "fanny"),
    naming = c("asis", "GOF", "size"),
    warn = TRUE,
    threshold = 0.01,
    ...) {
  segment <- Comp_GOF <- mean_GOF <- cluster_sort <- n <- segment <- cluster <- NULL

  naming <- match.arg(naming)
  method <- match.arg(method)

  stopifnot(inherits(x, "HeFTy") || inherits(x, "tTdiss") || inherits(x, "data.frame"))

  has_GOF <- inherits(x, "HeFTy") || inherits(x, "data.frame")
  if (has_GOF) {
    # dat <- x
    paths_GOF <- x$paths
    x <- path_diss(x, dist)
  }

  dmat <- x$diss
  paths <- x$paths

  cl <- switch(method,
    hclust = path_hcut(dmat, k = k, FUN = stats::hclust, ...)$cluster,
    diana = path_hcut(dmat, k = k, FUN = cluster::diana, ...)$cluster,
    agnes = path_hcut(dmat, k = k, FUN = cluster::agnes, ...)$cluster,
    kmeans = stats::kmeans(dmat, centers = k, ...)$cluster,
    pam = cluster::pam(dmat, k = k, ...)$clustering,
    dbscan = dbscan::dbscan(dmat, ...)$cluster,
    hdbscan = dbscan::hdbscan(dmat, ...)$cluster,
    specc = as.integer(kernlab::specc(as.matrix(dmat), centers = k, ...)),
    clara = cluster::clara(as.matrix(dmat), k = k, ...)$clustering,
    fanny = cluster::fanny(dmat, k = k, ...)$clustering
  )


  res <- data.frame(segment = paths$segment, cluster = as.factor(cl))

  if (naming == "GOF" && has_GOF) {
    outliers <- res |> dplyr::filter(as.integer(cluster) == 0)
    res <- res |>
      dplyr::filter(as.integer(cluster) > 0) |>
      dplyr::left_join(paths_GOF, dplyr::join_by(segment)) |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(
        mean_GOF = mean(Comp_GOF), .groups = "drop"
      ) |>
      dplyr::mutate(cluster_sort = as.factor(dplyr::min_rank(mean_GOF))) |>
      # dplyr::select(-mean_GOF) |>
      dplyr::right_join(res, dplyr::join_by(cluster)) |>
      dplyr::transmute(segment, cluster = cluster_sort) |>
      dplyr::bind_rows(outliers)
  } else if (naming == "size") {
    outliers <- res |> dplyr::filter(as.integer(cluster) == 0)
    res <- res |>
      dplyr::filter(as.integer(cluster) > 0) |>
      dplyr::count(cluster) |>
      dplyr::mutate(cluster_sort = as.factor(dplyr::min_rank(n))) |>
      dplyr::right_join(res, dplyr::join_by(cluster)) |>
      dplyr::transmute(segment, cluster = cluster_sort) |>
      dplyr::bind_rows(outliers)
  }

  if (isTRUE(warn)) {
    count.cluster <- count_cluster(res)

    n.cluster <- length(count.cluster)
    percentage.cluster <- count.cluster / sum(count.cluster)

    if (any(percentage.cluster < threshold)) {
      warning(paste0("Cluster with less than ", threshold * 100, "% of total paths detected \U1F622"))
    }
  }

  # if(has_GOF){
  #   dat$cluster <- dplyr::left_join(res, dat$paths, by = 'segment')
  #   res <- dat
  # }
  #
  res
}

#' Count the number of paths in each cluster
#'
#' Helper function to calculate the paths in each cluster and give a warning
#' message if one cluster contains not enough paths to be considered significant
#'
#' @param x data.frame with a `$cluster` column
#'
#' @export
#' @return named array of integers, the number of paths in each cluster
#' @examples
#' data(tT_paths1)
#' tT_paths1$paths <- subset(tT_paths1$paths, Comp_GOF >= 0.4)
#'
#' # cluster the paths
#' res <- cluster_paths(tT_paths1, k = 3)
#' count_cluster(res)
count_cluster <- function(x) split(x, x$cluster) |> sapply(function(i){nrow(unique(i))})


#' Hierarchical Clustering of Cooling Paths
#'
#' Convenience function to compute hierarchical clustering and cut the tree into k clusters
#'
#' @param x an object of class `"HeFTy"`, `"tTdiss`" or `"dist"` (dissimilarity matrix).
#' @param k integer. number of clusters to be generated
#' @param FUN hierarchical clustering function to be used, i.e. one of
#' [stats::hclust()] (the default), [cluster::agnes()], [cluster::diana()]).
#' @param ... optional arguments past to `hc_func`
#'
#' @returns an object of class `"hcut"` containing the result of the standard
#' function used (read the documentation of [stats::hclust()], [cluster::agnes()], [cluster::diana()]).
#'
#' It includes also
#' \describe{
#' \item{cluster}{the cluster assignment of observations after cutting the tree}
#' \item{nbclust}{the number of clusters}
#' \item{size}{the size of clusters}
#' }
#'
#' @importFrom stats hclust cutree
#' @importFrom cluster agnes diana
#'
#' @export
#'
#' @examples
#' data(tT_paths1)
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.4)
#'
#' # calculate the dissimilarities of the paths:
#' tT_diss <- path_diss(tT_paths_subset)
#' path_hcut(tT_diss$diss, 3)
path_hcut <- function(x, k, FUN = stats::hclust, ...) {
  if (inherits(x, "HeFTy")) x <- path_diss(x)
  if (inherits(x, "tTdiss")) x <- x$diss
  stopifnot(inherits(x, "dist") && is.numeric(k))

  hc <- FUN(x, ...)

  hc.cut <- stats::cutree(hc, k = k)
  hc$cluster <- hc.cut
  hc$nbclust <- k

  hc$size <- tabulate(hc.cut, nbins = k)
  class(hc) <- c(class(hc), "hcut")
  hc
}


CPL_geos_dist <- utils::getFromNamespace("CPL_geos_dist", "sf")

#' Helper function to calculate distances between paths
#'
#' @param x an `sf` linestring object with geometries of paths.
#' @param which character. Either `"Hausdorff"` (the default) or `"Frechet"`.
#' @param par numeric. tolerance parameter for distance calculation
#'
#' @returns a `dist` object with the pairwise distances between paths.
#' @noRd
path_distances <- function(x, which = c("Hausdorff", "Frechet"), par = 0) {
  x <- sf::st_geometry(x)
  CPL_geos_dist(x, x, which, par) |> as.dist()
}



#' Dissimilarity of thermochronology cooling paths
#'
#' Calculates the dissimilarity matrix and Hopkins statistic for
#' thermochronology cooling paths using the *Hausdorff* or the *Fr&#233;chet distance*.
#'
#' @param x either an object of class `"HeFTy"` (output of [read_hefty()]) or
#' a `data.frame` containing the `time`, `temperature`, and `segment` columns of the modeled paths.
#' @param dist character. Algorithm to calculate a dissimilarity matrix
#' (distance) for lines; either `Hausdorff` (the default) or `Frechet`.
#' @param densify numeric. optionally use a value
#' between 0 and 1 to densify the geometry description. Default is 0.
#' @param simplify numeric. Optional tolerance parameter to specify the
#' amount of simplification of path geometries in order to save processing time.
#' The tolerance is in the same unit as temperature and should relate to the
#' expected temperature spread.
#' The simplification is applied before the densification.
#' Some good values are 0.5, 2.5, 5, and 10. Default is 0 (no simplification).
#' @param ... (optional) arguments passed to [cluster_tendency()].
#'
#' @returns `tTdiss` object, i.e. a list containing
#' \describe{
#' \item{`paths`}{t-T paths as `sf` object in time-temperature space}
#' \item{`diss`}{the \eqn{n \times n} dissimilarity matrix (\eqn{n} is number of paths)}
#' \item{`method`}{the dissimilarity algorithm used}
#' \item{`hopkins`}{the Hopkins statistic and its p-value.}
#' \item{`mds`}{the multidimensional scaling coordinates of the dissimilarity matrix}
#' }
#'
#' @details The Hausdorff distance is the greatest of all the distances from a
#' point in one set to the closest point in the other set.
#' The Fr&#233;chet distance additionally takes into account the location and
#' ordering of the points along the curves (the "flow").
#'
#' @note The algorithm calculates the pairwise dissimilarities between the n
#' paths and creates a n x n matrix.
#' A large number of paths may end up in very long processing times and can
#' cause memory issues.
#' Thus, it is recommended to filter the paths beforehand using either
#' [crop_paths()] or [subset()] to, for example,
#' include only paths with large GOF (see examples),
#' or create a random subset using, e.g., [sample()] (see examples).
#' Additionally,
#' the `simplify` parameter can be adjusted to reduce the number of vertices
#' along the paths to speed up the processing time.
#'
#' @importFrom sf st_simplify st_distance st_as_sf st_cast sf_use_s2
#' @importFrom hopkins hopkins
#'
#' @export
#'
#' @examples
#' data(tT_paths1)
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.4)
#'
#' # calculate the dissimilarities of the paths:
#' tT_diss <- path_diss(tT_paths_subset, densify = 1)
#'
#' # the `diss` object of the output can be used for clustering using any
#' # available cluster algorithm, e.g.:
#' stats::kmeans(tT_diss$diss, centers = 3)
#'
#' # using a random subset of paths (in case of too large data)
#' set.seed(20250411)
#'
#' ## select 100 random path segments:
#' random_segments <- sample(unique(tT_paths1$paths$segment), size = 100)
#' tT_paths_rnd <- subset(tT_paths1$paths, segment %in% random_segments)
#'
#' tT_diss_rnd <- path_diss(tT_paths_rnd)
path_diss <- function(x, dist = c("Hausdorff", "Frechet"), densify = 0, simplify = 0, ...) {
  if (inherits(x, "HeFTy")) x <- x$paths
  stopifnot(inherits(x, "data.frame") & c("time", "temperature", "segment") %in% colnames(x))

  segment <- NULL
  dist <- match.arg(dist)

  suppressMessages({
    sf::sf_use_s2(FALSE)
    paths <- x |>
      sf::st_as_sf(coords = c("time", "temperature")) |>
      dplyr::group_by(segment) |>
      dplyr::summarise(do_union = FALSE) |>
      sf::st_cast("LINESTRING") |>
      sf::st_simplify(dTolerance = simplify)

    dmat <- path_distances(paths, which = dist, par = densify)
    sf::sf_use_s2(TRUE)
  })

  h <- cluster_tendency(dmat, ...)

  paths_mds <- stats::cmdscale(dmat)
  # paths_mds <- cbind(MDS1 = paths_mds[, 1], MDS2 = paths_mds[, 2])
  rownames(paths_mds) <- paths$segment



  res <- list(paths = paths, diss = dmat, hopkins = h, dist = dist, mds = paths_mds)
  class(res) <- append(class(res), "tTdiss")
  return(res)
}

# path_diss_fast <- function(x, dist = c("Hausdorff", "Frechet"), densify = 0, simplify = 0, ...) {
#   if (inherits(x, "HeFTy")) x <- x$paths
#   stopifnot(inherits(x, "data.frame"))
#
#   segment <- NULL
#   dist <- match.arg(dist)
#
#   suppressMessages({
#     sf::sf_use_s2(FALSE)
#     paths <- x |>
#       sf::st_as_sf(coords = c("time", "temperature")) |>
#       dplyr::group_by(segment) |>
#       dplyr::summarise(do_union = FALSE) |>
#       sf::st_cast("LINESTRING") |>
#       sf::st_simplify(dTolerance = simplify)
#     sf::sf_use_s2(TRUE)
#   })
#
#   paths_mat <- split(paths, paths$segment) |>
#     lapply(sf::st_coordinates)
#
#   dmat <- hausdorff_dmat(paths_mat)
#
#   h <- cluster_tendency(dmat, ...)
#
#   res <- list(paths = paths, diss = dmat, hopkins = h, dist = dist)
#   class(res) <- append(class(res), "tTdiss")
#   return(res)
# }

#' Cluster tendency of thermochronology cooling paths
#'
#' Calculate the Hopkins statistic of Hausdorff or Fr&#233;chet distance matrices
#' to check clusterability of thermochronologic cooling paths.
#' Calculated values 0-0.3 indicate regularly-spaced data. Values around 0.5
#' indicate random data. Values 0.7-1 indicate clustered data.
#'
#' @param x either an object of class `"tTdiss"` (output of [path_diss()]) or a
#' distance matrix.
#' @param m integer. Number of rows to sample from `x`. If `NULL`, all
#' rows are considered.
#' @param method character. Either `"simple"` (default) or `"torus"`.
#'
#' @returns numeric. The value of Hopkins statistic in the range of 0 and 1 and
#' its p-value (under the null hypothesis of spatial randomness).
#'
#' @importFrom hopkins hopkins hopkins.pval
#' @importFrom stats cmdscale
#'
#' @details Calculates the Hopkins statistic on the transformed
#' Hausdorff or Fr&#233;chet distance matrix using metric multidimensional scaling.
#'
#'
#' @export
#'
#' @examples
#' data(tT_paths1)
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.4)
#'
#' # calculate the dissimilarities of the paths:
#' tT_diss <- path_diss(tT_paths_subset, densify = 1)
#'
#' set.seed(20250411)
#' cluster_tendency(tT_diss$diss) # H=0.68, p=8.6e-7
cluster_tendency <- function(x, m = NULL, method = c("simple", "torus")) {
  if (inherits(x, "tTdiss")) x <- x$diss

  mds <- stats::cmdscale(x)
  if (is.null(m)) {
    m <- nrow(mds) - 1
    # m <- ceiling(nrow(mds)/10)
  }
  method <- match.arg(method)

  h <- hopkins::hopkins(mds, m = m, method = method)
  p <- hopkins::hopkins.pval(h, m)
  stats::setNames(c(h, p), c("statistic", "p-value"))
}
