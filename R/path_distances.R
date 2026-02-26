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
#' @name path_diss
#'
#' @export
#'
#' @examples
#' data(tT_paths)
#' tT_paths_subset <- subset(tT_paths$paths, Comp_GOF >= 0.4)
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
#' random_segments <- sample(unique(tT_paths$paths$segment), size = 100)
#' tT_paths_rnd <- subset(tT_paths$paths, segment %in% random_segments)
#'
#' path_diss(tT_paths_rnd)
#' 
#' \dontrun{
#' use_parallel(TRUE)
#' path_diss_parallel(tT_paths_rnd)
#' use_parallel(FALSE)
#' }
NULL

#' @rdname path_diss
path_diss <- function(x, dist = c("Hausdorff", "Frechet"), densify = 0, simplify = 0, ...) {
  if (inherits(x, "HeFTy")) x <- x$paths
  stopifnot(inherits(x, "data.frame") & c("time", "temperature", "segment") %in% colnames(x))
  dist <- match.arg(dist)
  
  x_dist <- path_distances_R(x, dist, densify, simplify)
  
  paths <- x_dist$paths
  dmat <- x_dist$dmat
  
  # Hopkins statistics
  h <- cluster_tendency(dmat, ...)
  
  # MDS
  paths_mds <- stats::cmdscale(dmat)
  rownames(paths_mds) <- paths$segment
  
  res <- list(paths = paths, diss = dmat, hopkins = h, dist = dist, mds = paths_mds)
  class(res) <- append(class(res), "tTdiss")
  return(res)
}

#' @rdname path_diss
path_diss_parallel <- function(x, simplify = 0, ...) {
  if (inherits(x, "HeFTy")) x <- x$paths
  stopifnot(inherits(x, "data.frame") & c("time", "temperature", "segment") %in% colnames(x))
  # dist <- match.arg(dist)
  
  x_dist <- path_distances_R_parallel(x, simplify)
  
  paths <- x_dist$paths
  dmat <- x_dist$dmat
  
  # Hopkins statistics
  h <- cluster_tendency(dmat, ...)
  
  # MDS
  paths_mds <- stats::cmdscale(dmat)
  rownames(paths_mds) <- paths$segment
  
  res <- list(paths = paths, diss = dmat, hopkins = h, dist = dist, mds = paths_mds)
  class(res) <- append(class(res), "tTdiss")
  return(res)
}


# path_dissC <- function(x, dist = c("Hausdorff", "Frechet"), simplify = 0, ...) {
#   if (inherits(x, "HeFTy")) x <- x$paths
#   stopifnot(inherits(x, "data.frame"))
#   
#   x_dist <- path_distances_C(x, dist, simplify)
#   
#   paths <- x_dist$paths
#   dmat <- x_dist$dmat
#   
#   # Hopkins statistics
#   h <- cluster_tendency(dmat, ...)
#   
#   # MDS
#   paths_mds <- stats::cmdscale(dmat)
#   rownames(paths_mds) <- paths$segment
#   
#   res <- list(paths = paths, diss = dmat, hopkins = h, dist = dist, mds = paths_mds)
#   class(res) <- append(class(res), "tTdiss")
#   return(res)
# }


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

path_distances_R <- function(x, dist = c("Hausdorff", "Frechet"), densify = 0, simplify = 0){
  segment <- NULL
  dist <- match.arg(dist)
  
  suppressMessages({
    sf::sf_use_s2(FALSE)
    paths <- paths2lines(x, simplify)
    
    dmat <- path_distances(paths, which = dist, par = densify)
    sf::sf_use_s2(TRUE)
  })
  
  return(list(paths = paths, dmat = dmat))
}

#' @importFrom future.apply future_lapply
path_distances_R_parallel <- function(x, simplify = 0, ...){
  segment <- NULL
  # dist <- match.arg(dist)
  
  suppressMessages({
    sf::sf_use_s2(FALSE)
    paths <- paths2lines(x, simplify)
  })
  
  path_list <- 
    split(x, x$segment) |> 
    future.apply::future_lapply(function(x) {
       as.matrix(select(x, time, temperature))
    })
  
  dmat <- hausdorff_list_parallel(path_list)
  
  return(list(paths = paths, dmat = dmat))
}


# path_distances_C <- function(x, dist = c("Hausdorff", "Frechet"), simplify = 0){
#   segment <- NULL
#   
#   suppressMessages({
#     sf::sf_use_s2(FALSE)
#     paths <- paths2lines(x, simplify)
#     sf::sf_use_s2(TRUE)
#   })
#   paths_mat <- split(paths, paths$segment) |>
#     lapply(sf::st_coordinates) |> 
#     lapply(function(i){
#       as.matrix(select(x, time, temperature))
#     })
#   
#   dist <- match.arg(dist)
#   dmat <- switch(dist,
#     Hausdorff = hausdorff_list_C(paths_mat),
#     Frechet = disc_frechet_list_C(paths_mat)
#   )
# 
#   return(list(paths = paths, dmat = dmat))
# }


paths2lines <- function(x, simplify = 0){
  sf::st_as_sf(x, coords = c("time", "temperature")) |>
    dplyr::group_by(segment) |>
    dplyr::summarise(do_union = FALSE) |>
    sf::st_cast("LINESTRING") |>
    sf::st_simplify(dTolerance = simplify)
}



#' Hausdorff Distance
#'
#' Hausdorff distance (aka Hausdorff dimension)
#'
#' @param P,Q numerical matrices, representing points in an m-dim. space.
#'
#' @returns A single scalar, the Hausdorff distance (dimension).
#'
#' @name hausdorff
#'
#' @examples
#' P <- matrix(c(1, 1, 2, 2, 5, 4, 5, 4), 4, 2)
#' Q <- matrix(c(4, 4, 5, 5, 2, 1, 2, 1), 4, 2)
#' hausdorff(P, Q)
NULL

#' @rdname hausdorff
#' @export
hausdorff <- function(P, Q) {
  max(
    directed_hausdorff(P, Q),
    directed_hausdorff(Q, P)
  )
}

#' @rdname hausdorff
#' @export
#' @importFrom RANN nn2
directed_hausdorff <- function(P, Q) {
  nn <- RANN::nn2(Q, query = P, k = 1)
  max(nn$nn.dists)
}


#' Hausdorff distance matrix
#'
#' @param x list of matrices or vectors, each representing a segment of a path.
#'
#' @returns distance matrix
#' 
#' @noRd 
#' @importFrom future.apply future_sapply
#'
#' @examples
#' tT_paths_subset <- subset(tT_paths$paths, Comp_GOF >= 0.2)
#' tT_paths_list <- split(tT_paths_subset, tT_paths_subset$segment) |> lapply(function(x) {
#'   as.matrix(select(x, time, temperature))
#' })
#' hausdorff_list_parallel(tT_paths_list)
hausdorff_list_parallel <- function(x) {
  keep <- vapply(x, function(m)
    !is.null(m) && length(m) > 0L,
    logical(1))
  
  mats <- x[keep]
  n <- length(mats)
  
  combs <- which(lower.tri(matrix(0, n, n)), arr.ind = TRUE)
  
  d <- future_sapply(
    seq_len(nrow(combs)),
    function(k) {
      i <- combs[k,1]
      j <- combs[k,2]
      hausdorff(mats[[i]], mats[[j]])
    }
  )
  
  out <- matrix(0, n, n)
  out[lower.tri(out)] <- d
  
  as.dist(out)
}

