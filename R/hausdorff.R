#' Hausdorff Distance
#' 
#' Hausdorff distance (aka Hausdorff dimension)
#'
#' @param P,Q numerical matrices, representing points in an m-dim. space.
#'
#' @returns A single scalar, the Hausdorff distance (dimension).
#' @export
#' @noRd
#' 
#' @importFrom matrixStats rowMins colMins
#'
#' @examples
#' P <- matrix(c(1,1,2,2, 5,4,5,4), 4, 2)
#' Q <- matrix(c(4,4,5,5, 2,1,2,1), 4, 2)
#' hausdorff(P, Q)    # 4.242641 = sqrt(sum((c(4,2)-c(1,5))^2))
hausdorff <- function(P, Q) {
  stopifnot(is.numeric(P), is.numeric(Q))
  
  if (is.vector(P)) 
    P <- matrix(P, ncol = 1)
  if (is.vector(Q)) 
    Q <- matrix(Q, ncol = 1)
  if (ncol(P) != ncol(Q)) 
    stop("'P' and 'Q' must have the same number of columns.")
  
  D <- distmat(P, Q)
  
  dhd_PQ <- max(matrixStats::rowMins(D))
  dhd_QP <- max(matrixStats::colMins(D))
  
  max(dhd_PQ, dhd_QP)
}


#' Distance Matrix
#' 
#' Computes the Euclidean distance between rows of two matrices.
#'
#' @param X matrix of some size m x k; vector will be taken as row matrix.
#' @param Y matrix of some size m x k; vector will be taken as row matrix.
#'
#' @returns matrix
#' @noRd
#'
#' @examples
#' A <- c(0.0, 0.0)
#' B <- matrix(c(0,0, 1,0, 0,1, 1,1), nrow=4, ncol = 2, byrow = TRUE)
#' distmat(A, B)  #=> 0 1 1 sqrt(2)
#' distmat_fast(A, B)  #=> 0 1 1 sqrt(2)
distmat <- function(X, Y) {
  XX <- rowSums(X^2)
  YY <- rowSums(Y^2)
  XY <- X %*% t(Y)
  
  # Outer sum of XX and YY + subtract 2*XY
  d2 <- outer(XX, YY, `+`) - 2 * XY
  sqrt(pmax(d2, 0))
}


#' Hausdorff distance matrix
#'
#' @param x list of matrices or vectors, each representing a segment of a path.
#' @param workers integer. Number of parallel workers to use for computation.
#'
#' @returns matrix
#' @export
#' 
#' @importFrom future.apply future_lapply
#' @import future
#' 
#' @examples
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.2)
#' tT_paths_list <- split(tT_paths_subset, tT_paths_subset$segment) |> lapply(function(x) {as.matrix(select(x, time, temperature))})
#' hausdorff_dmat(tT_paths_list)
#' hausdorff_dmat2(tT_paths_list)
hausdorff_dmat <- function(x){
x1 <- x[lapply(x,length)>0]
n <- length(x1)
seg_names <- names(x1)
x2 <- unname(x1)

combs <- which(lower.tri(matrix(NA, n, n)), arr.ind = TRUE)

hdists <- lapply(seq_len(nrow(combs)), function(k) {
  # for(k in seq_len(nrow(combs))){
  i <- combs[k, 1]
  j <- combs[k, 2]
  d <- hausdorff_rcpp(x2[[i]], x2[[j]])
  c(i = i, j = j, d = d)
}
)

hdist_mat  <- do.call(rbind, hdists)
mat <- matrix(0, n, n, dimnames = list(seg_names, seg_names))
mat[cbind(hdist_mat[, "i.row"], hdist_mat[, "j.col"])] <- hdist_mat[, "d"]
mat[cbind(hdist_mat[, "j.col"], hdist_mat[, "i.row"])] <- hdist_mat[, "d"]

as.dist(mat)
}

hausdorff_dmat2 <- function(x){
  x1 <- x[lapply(x,length)>0]
  n <- length(x1)
  seg_names <- names(x1)
  x2 <- unname(x1)
  
  hausdorff_matrix_rcpp(x2) |>  as.dist()
}
