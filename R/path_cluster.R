#' Cluster thermal histories
#'
#' Groups t-T paths into "path families" based on the *Hausdorff distance*
#' between paths.
#'
#' @param x either an object of class `"HeFTy"` (output of [read_hefty()]), an 
#' object of `"tTdiss"` (output of [path_diss()]), 
#' or a `data.frame` containing the `time`, `temperature` columns of the 
#' modeled paths.
#' @param cluster an integer scalar or vector with the desired number of groups.
#' Ignored when `dist` equal to `dbscan` or `hdbscan` (**see Note**).
#' @param dist character. Algorithm to calculate a dissimilarity matrix
#' (distance) for lines; one of `Hausdorff` (the default),
#' `Frechet`, or `Euclidean`.
#' @param method character. Clustering method to use. Currently implemented are
#' \describe{
#'  \item{`"hclust"`}{for Hierarchical Clustering using [stats::hclust()], the default)}
#'  \item{`"kmeans"`}{for K-Means Clustering using [stats::kmeans()])}
#'  \item{`"pam"`}{for Partitioning Around Medoids using [cluster::pam()]}
#'  \item{`"dbscan"`}{for Density-based Spatial Clustering of Applications with
#'  Noise (DBSCAN) using [dbscan::dbscan()] (**see Note**)}
#'  \item{`"hdbscan"`}{for Hierarchical DBSCAN using [dbscan::hdbscan()]
#'  (**see Note**)}
#' }
#' @param ... additional arguments passed to cluster method.
#'
#' @note that `dbscan` and `hdbscan` methods require `eps` and `minPts` arguments.
#' Optimal `eps` values can be visually estimated from the "knee" in a k-nearest
#' neighbor distance plot using:
#' [dbscan::kNNdistplot()]
#'
#' @return a tibble with the path `segment` (integer) and `cluster` (factor)
#' @export
#'
#' @importFrom sf st_as_sf st_distance
#' @importFrom dplyr summarise group_by tibble
#' @importFrom stats hclust cutree as.dist kmeans
#' @importFrom cluster pam
#' @importFrom forcats as_factor
#' @importFrom dbscan dbscan hdbscan
#'
#' @seealso [path_diss()]
#'
#' @examples
#' \dontrun{
#' # example data
#' data(tT_paths1)
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.5)
#'
#' # cluster the paths
#' cluster_paths(tT_paths_subset, cluster = 3)
#' }
cluster_paths <- function(
    x, cluster,
    dist = c("Hausdorff", "Frechet", "Euclidean"),
    method = c("hclust", "kmeans", "pam", "dbscan", "hdbscan"),
    ...) {
  if (inherits(x, "HeFTy")) x <- x$paths
  if (!inherits(x, "tTdiss")) x <- path_diss(x, dist)
  dmat <- x$diss
  paths <- x$paths

  method <- match.arg(method)

  if (method == "hclust") {
    cl <- stats::hclust(stats::as.dist(dmat), ...) |>
      stats::cutree(k = cluster)
  } else if (method == "kmeans") {
    cl <- stats::kmeans(dmat, centers = cluster, ...)$cluster
  } else if (method == "pam") {
    cl <- cluster::pam(dmat, k = cluster, ...)$clustering
  } else if (method == "dbscan") {
    cl <- dbscan::dbscan(dmat, ...)$cluster
  } else if (method == "hdbscan") {
    cl <- dbscan::hdbscan(dmat, ...)$cluster
  }

  dplyr::tibble(segment = paths$segment, cluster = forcats::as_factor(cl))
}




#' Dissimilarity matrix for t-T paths
#'
#' @param x either an object of class `"HeFTy"` (output of [read_hefty()]) or 
#' a `data.frame` containing the `time`, `temperature` columns of the modeled paths.
#' @param method character. Algorithm to calculate a dissimilarity matrix
#' (distance) for lines; one of `Hausdorff` (the default),
#' `Frechet`, or `Euclidean`.
#' @param densify for `dist` equal to `Hausdorff` or `Frechet`, optionally use a value
#' between 0 and 1 to densify the geometry description.
#'
#' @returns `tTdiss` object, i.e. a list containing
#' \describe{
#' \item{`paths`}{t-T paths as `sf` object in time-temperature space}
#' \item{`diss`}{the \eqn{n \times n} dissimilarity matrix (\eqn{n} is number of paths)}
#' }
#'
#' @details The Hausdorff distance is the greatest of all the distances from a
#' point in one set to the closest point in the other set.
#' The Fréchet distance additionally takes into account the location and
#' ordering of the points along the curves (the "flow").
#'
#' @export
#'
#' @examples
#' # example data
#' data(tT_paths1)
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.5)
#'
#' # calculate the dissimilarities of the paths:
#' tT_diss <- path_diss(tT_paths_subset, densify = 1)
#' print(tT_diss)
#'
#' # the `diss` object of the output can be used for clustering, e.g.:
#' stats::kmeans(tT_diss$diss, centers = 3)
path_diss <- function(x, method = c("Hausdorff", "Frechet", "Euclidean"), densify = 0) {
  if(inherits(x, 'HeFTy')) x <- x$paths
  stopifnot(inherits(x, "data.frame"))
  
  segment <- NULL
  method <- match.arg(method)

  paths <- x |>
    sf::st_as_sf(coords = c("time", "temperature")) |>
    dplyr::group_by(segment) |>
    dplyr::summarise(do_union = FALSE) |>
    sf::st_cast("LINESTRING")

  dmat <- sf::st_distance(paths, which = method, par = densify)

  res <- list(paths = paths, diss = dmat, diss = method)
  class(res) <- append(class(res), "tTdiss")
  return(res)
}
