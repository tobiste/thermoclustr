
#' Cluster tendency of thermochronology cooling paths
#'
#' Calculate the Hopkins statistic of Hausdorff or Fr&#233;chet distance matrices
#' to check clusterability of thermochronologic cooling paths.
#' Calculated values 0-0.3 indicate regularly-spaced data. Values around 0.5
#' indicate random data. Values 0.7-1 indicate clustered data.
#'
#' @param x either an object of class `"tTdiss"` (output of [path_diss()]), a
#' distance matrix (object `"dist"`), or a `data.frame` containing coordinates from a dimension-reducing algorithm.
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
#' Hausdorff or Fr&#233;chet distance matrix using metric multidimensional scaling by default.
#' 
#'
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
#' set.seed(20250411)
#' cluster_tendency(tT_diss) # H=0.77, p=2.88e-11
#' 
#' # Alternative dimension reducing algorithms can be used:
#' @examplesIf require("MASS")
#' # Kriskal's non-metric MDS
#' library(MASS)
#' tT_nmds <- isoMDS(tT_diss$diss)$points
#' cluster_tendency(tT_nmds)
#'
#' @examplesIf require("uwot")
#' # Uniform Manifold Approximation and Projection (UMAP)  
#' library(uwot)
#' tT_umap <- umap2(tT_diss$diss)
#' cluster_tendency(tT_umap)
cluster_tendency <- function(x, m = NULL, method = c("simple", "torus")) {
  if (inherits(x, "tTdiss")) x <- x$diss
  
  if(inherits(x, 'dist')) {
    mds <- stats::cmdscale(x)
  } else {
    mds <- x
  }
  
  if (is.null(m)) {
    m <- nrow(mds) - 1
    # m <- ceiling(nrow(mds)/10)
  }
  method <- match.arg(method)
  
  h <- hopkins::hopkins(mds, m = m, method = method)
  p <- hopkins::hopkins.pval(h, m)
  stats::setNames(c(h, p), c("statistic", "p-value"))
}