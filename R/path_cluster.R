#' Clustering thermal histories models
#'
#' Groups t-T paths into "path families" based on the *Hausdorff* or
#' *Fréchet distance* between paths.
#'
#' @param x either an object of class `"HeFTy"` (output of [read_hefty()]), an
#' object of `"tTdiss"` (output of [path_diss()]),
#' or a `data.frame` containing the `time`, `temperature` columns of the
#' modeled paths.
#' @param cluster an integer scalar or vector with the desired number of groups.
#' Ignored when `dist` equal to `dbscan` or `hdbscan` (**see Note**).
#' @param dist character. Algorithm to calculate a dissimilarity matrix
#' (distance) for lines; either `Hausdorff` (the default), or
#' `Frechet`.
#' @param method character. Clustering method to use. Currently implemented are
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
#' @importFrom dplyr summarise group_by tibble left_join right_join join_by count select min_rank row_number mutate
#' @importFrom stats hclust cutree as.dist kmeans
#' @importFrom cluster pam agnes diana clara fanny
#' @importFrom kernlab specc
#' @importFrom forcats as_factor
#' @importFrom dbscan dbscan hdbscan
#'
#' @seealso [path_diss()]
#'
#' @examples
#' # example data
#' data(tT_paths1)
#' tT_paths1$paths <- subset(tT_paths1$paths, Comp_GOF >= 0.2)
#'
#' # cluster the paths
#' cluster_paths(tT_paths1, cluster = 3)
cluster_paths <- function(
    x, cluster,
    dist = c("Hausdorff", "Frechet"),
    method = c("hclust", "kmeans", "pam", "specc", "dbscan", "hdbscan", "diana", "agnes", "clara", "fanny"),
    naming = c("asis", "GOF", "size"),
    ...) {
  naming <- match.arg(naming)
  method <- match.arg(method)
  stopifnot(inherits(x, "HeFTy") | inherits(x, "tTdiss"))

  has_GOF <- FALSE
  if (inherits(x, "HeFTy")) {
    has_GOF <- TRUE
    paths_GOF <- x$paths
    x <- path_diss(x, dist) # calculate dissimilarities
  } else if (inherits(x, "tTdiss")) {
    has_GOF <- FALSE
  }

  dmat <- x$diss
  paths <- x$paths

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
  } else if (method == "specc") {
    cl <- kernlab::specc(dmat, centers = cluster, ...) |>
      as.integer()
  } else if (method == "diana") {
    cl <- cluster::diana(dmat, ...) |>
      stats::cutree(k = cluster)
  } else if (method == "agnes") {
    cl <- cluster::agnes(dmat, ...) |>
      stats::cutree(k = cluster)
  } else if (method == "clara") {
    cl <- cluster::clara(dmat, k = cluster, ...)$clustering
  } else if (method == "fanny") {
    cl <- cluster::fanny(dmat, k = cluster, ...)$clustering
  }

  res <- dplyr::tibble(segment = paths$segment, cluster = forcats::as_factor(cl))

  if (naming == "GOF" & has_GOF) {
    outliers <- res |> dplyr::filter(as.numeric(as.character(cluster)) == 0)
    res <- res |>
      dplyr::filter(as.numeric(as.character(cluster)) > 0) |>
      dplyr::left_join(paths_GOF, dplyr::join_by(segment)) |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(
        mean_GOF = mean(Comp_GOF)
      ) |>
      dplyr::mutate(cluster_sort = forcats::as_factor(dplyr::min_rank(mean_GOF))) |>
      # dplyr::select(-mean_GOF) |>
      dplyr::right_join(res, dplyr::join_by(cluster)) |>
      dplyr::select(segment, cluster = cluster_sort) |>
      dplyr::bind_rows(outliers)
  } else if (naming == "size") {
    outliers <- res |> dplyr::filter(as.numeric(as.character(cluster)) == 0)
    res <- res |>
      dplyr::filter(as.numeric(as.character(cluster)) > 0) |>
      dplyr::count(cluster) |>
      dplyr::mutate(cluster_sort = forcats::as_factor(dplyr::min_rank(n))) |>
      dplyr::right_join(res, dplyr::join_by(cluster)) |>
      dplyr::select(segment, cluster = cluster_sort) |>
      dplyr::bind_rows(outliers)
  }

  res
}




#' Dissimilarity of thermochronologic cooling paths
#'
#' Calculates the dissimilarity matrix and Hopkins statistic for
#' thermochronologic cooling paths using the *Hausdorff* or the *Fréchet distance*.
#'
#' @param x either an object of class `"HeFTy"` (output of [read_hefty()]) or
#' a `data.frame` containing the `time`, `temperature` columns of the modeled paths.
#' @param dist character. Algorithm to calculate a dissimilarity matrix
#' (distance) for lines; either `Hausdorff` (the default) or `Frechet`.
#' @param densify numeric. optionally use a value
#' between 0 and 1 to densify the geometry description.
#' @param ... (optional) arguments passed to [cluster_tendency()].
#'
#' @returns `tTdiss` object, i.e. a list containing
#' \describe{
#' \item{`paths`}{t-T paths as `sf` object in time-temperature space}
#' \item{`diss`}{the \eqn{n \times n} dissimilarity matrix (\eqn{n} is number of paths)}
#' \item{`method`}{the dissimilarity algorithm used}
#' \item{`hopkins`}{the Hopkins statistic and its p-value.}
#' }
#'
#' @details The Hausdorff distance is the greatest of all the distances from a
#' point in one set to the closest point in the other set.
#' The Fréchet distance additionally takes into account the location and
#' ordering of the points along the curves (the "flow").
#'
#' @importFrom hopkins hopkins
#'
#' @export
#'
#' @examples
#' # example data
#' data(tT_paths1)
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.2)
#'
#' # calculate the dissimilarities of the paths:
#' tT_diss <- path_diss(tT_paths_subset, densify = 1)
#' print(tT_diss)
#'
#' # the `diss` object of the output can be used for clustering, e.g.:
#' stats::kmeans(tT_diss$diss, centers = 3)
path_diss <- function(x, dist = c("Hausdorff", "Frechet"), densify = 0, ...) {
  if (inherits(x, "HeFTy")) x <- x$paths
  stopifnot(inherits(x, "data.frame"))

  segment <- NULL
  dist <- match.arg(dist)

  paths <- x |>
    sf::st_as_sf(coords = c("time", "temperature")) |>
    dplyr::group_by(segment) |>
    dplyr::summarise(do_union = FALSE) |>
    sf::st_cast("LINESTRING")

  dmat <- sf::st_distance(paths, which = dist, par = densify) |> 
    as.dist()

  h <- cluster_tendency(dmat, ...)

  res <- list(paths = paths, diss = dmat, hopkins = h, dist = dist)
  class(res) <- append(class(res), "tTdiss")
  return(res)
}


#' Cluster tendency of thermochronologic cooling paths
#'
#' Calculate the Hopkins statistic of Hausdorff or Fréchet distance matrices
#' to check clusterability of thermochronologic cooling paths.
#' If the value is close to 1 (far above 0.5), then the
#' dataset is significantly clusterable.
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
#' @details The algorithm calculates the Hopkins statistic on a transformed
#' Hausdorff or Fréchet distance matrix using metric multidimensional scaling.
#'
#'
#' @export
#'
#' @examples
#' # example data
#' data(tT_paths1)
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.4)
#'
#' # calculate the dissimilarities of the paths:
#' tT_diss <- path_diss(tT_paths_subset, densify = 1)
#'
#' set.seed(1)
#' cluster_tendency(tT_diss$diss) # H=0.867, p=0.003
cluster_tendency <- function(x, m = NULL, method = c("simple", "torus")) {
  if (inherits(x, "tTdiss")) x <- x$diss

  mds <- stats::cmdscale(x)
  if (is.null(m)) {
    m <- nrow(mds)-1
    # m <- ceiling(nrow(mds)/10)
  }
  method <- match.arg(method)

  h <- hopkins::hopkins(mds, m = m, method = method)
  p <- hopkins::hopkins.pval(h, m)
  setNames(c(h, p), c("statistic", "p-value"))
}
