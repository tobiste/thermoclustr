#' Cluster thermal histories
#'
#' Groups t-T paths into "path families" based on the *Hausdorff distance*
#' between paths.
#'
#' @param x t-T and GOF data of the modeled paths. Output of [read_hefty()].
#' @param cluster an integer scalar or vector with the desired number of groups
#' @param dist character; algorithm to calculate a dissimilarity matrix
#' (distance) for lines; one of `Euclidean`, `Hausdorff` (the default) or
#' `Frechet`.
#' @param method character. Clustering method to use; one of `"hclust"`
#' (for [stats::hclust()], the default), `"kmeans"` (for [stats::kmeans()]),
#' or `"pam"` (for [cluster::pam()].
#' @param ... additional arguments passed to cluster method.
#'
#' @return a tibble with the path `segment` (integer) and `cluster` (factor)
#' @export
#'
#' @importFrom sf st_as_sf st_distance
#' @importFrom dplyr summarise group_by tibble
#' @importFrom stats hclust cutree as.dist kmeans
#' @importFrom cluster pam
#' @importFrom forcats as_factor
#'
#' @examples
#' \dontrun{
#' data(tT_paths1)
#' tT_paths_subset <- subset(tT_paths1, Comp_GOF >= 0.5)
#' cluster_paths(tT_paths_subset, cluster = 3)
#' }
cluster_paths <- function(x, cluster, dist = c("Hausdorff", "Frechet", "Euclidean"), method = c("hclust", "kmeans", 'pam'), ...) {
  segment <- NULL
  x_lines <- x |>
    sf::st_as_sf(coords = c("time", "temperature")) |>
    dplyr::group_by(segment) |>
    dplyr::summarise(do_union = FALSE) |>
    sf::st_cast("LINESTRING")

  dist <- match.arg(dist)
  dmat <- sf::st_distance(x_lines, which = dist)

  method <- match.arg(method)

  if (method == "hclust") {
    cl <- stats::hclust(stats::as.dist(dmat), ...) |>
      stats::cutree(k = cluster)
  } else if (method == "kmeans") {
    cl <- stats::kmeans(dmat, centers = cluster, ...)$cluster
  } else if(method == 'pam') {
    cl <- cluster::pam(dmat, k = cluster, ...)$clustering
  }

  dplyr::tibble(segment = x_lines$segment, cluster = forcats::as_factor(cl))
}
