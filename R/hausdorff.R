#' Hausdorff Distance
#'
#' Hausdorff distance (aka Hausdorff dimension)
#'
#' @param P,Q numerical matrices, representing points in an m-dim. space.
#' @param check logical. Checks if input arguments are valid. 
#'
#' @returns A single scalar, the Hausdorff distance (dimension).
#' @export
#'
#' @importFrom matrixStats rowMins colMins
#'
#' @examples
#' P <- matrix(c(1, 1, 2, 2, 5, 4, 5, 4), 4, 2)
#' Q <- matrix(c(4, 4, 5, 5, 2, 1, 2, 1), 4, 2)
#' hausdorff(P, Q) # 4.242641 = sqrt(sum((c(4,2)-c(1,5))^2))
hausdorff <- function(P, Q, check = TRUE) {
  if(isTRUE(check)){
  if (is.vector(P)) P <- matrix(P, ncol = 1)
  if (is.vector(Q)) Q <- matrix(Q, ncol = 1)

  if (ncol(P) != ncol(Q)) {
    stop("'P' and 'Q' must have the same number of columns.")
  }

  if (nrow(P) == 0L || nrow(Q) == 0L)
    return(NA_real_)

  stopifnot(is.numeric(P), is.numeric(Q))
  }

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
#' B <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, ncol = 2, byrow = TRUE)
#' distmat(A, B) # => 0 1 1 sqrt(2)
#' distmat_fast(A, B) # => 0 1 1 sqrt(2)
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
#'
#' @returns matrix
#' @name hausdorff_dmat
#'
#' @examples
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.2)
#' tT_paths_list <- split(tT_paths_subset, tT_paths_subset$segment) |> lapply(function(x) {
#'   as.matrix(select(x, time, temperature))
#' })
#' hausdorff_list_R(tT_paths_list)
#' hausdorff_list_C(tT_paths_list)
NULL

# uses R functions
#' @rdname hausdorff_dmat
#' @export
hausdorff_list_R <- function(x) {
  # drop empty list elements
  keep <- vapply(x, function(m) !is.null(m) && length(m) > 0L, logical(1))
  mats <- x[keep]

  n <- length(mats)
  if (n == 0)
    return(as.dist(matrix(numeric(0), 0, 0)))

  out <- matrix(NA_real_, n, n)

  for (i in seq.int(n)) {
    Pi <- mats[[i]]
    for (j in seq_len(i - 1)) {
      out[i, j] <- hausdorff(Pi, mats[[j]], check = FALSE)
    }
  }

  as.dist(out)
}

# uses Cpp functions
#' @rdname hausdorff_dmat
#' @export
hausdorff_list_C <- function(x) {
  keep <- vapply(x, function(m) !is.null(m) && length(m) > 0L, logical(1))
  mats <- x[keep]
  
  D_cpp <- pairwise_hausdorff_cpp(mats)
  
  as.dist(D_cpp)
}


#' Frechet distance matrix
#' 
#' Calculates the canonical discrete Fréchet algorithm (standard for sampled paths) after Eiter & Mannila (1994) or the 
#' continuous Fréchet distance (free-space diagram algorithm) after Alt & Godau (1995).
#'
#' @inheritParams hausdorff_dmat_C
#'
#' @returns matrix
#' @name frechet_dmat
#'
#' @examples
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.2)
#' tT_paths_list <- split(tT_paths_subset, tT_paths_subset$segment) |> lapply(function(x) {
#'   as.matrix(select(x, time, temperature))
#' })
#' cont_frechet_list_C(tT_paths_list)
#' disc_frechet_list_C(tT_paths_list)
NULL

#' @rdname frechet_dmat
#' @export
cont_frechet_list_C <- function(x) {
  keep <- vapply(x, function(m) !is.null(m) && length(m) > 0L, logical(1))
  mats <- x[keep]
  
  D_cpp <- pairwise_cont_frechet_cpp(mats)
  
  as.dist(D_cpp)
}

#' @rdname frechet_dmat
#' @export
disc_frechet_list_C <- function(x) {
  keep <- vapply(x, function(m) !is.null(m) && length(m) > 0L, logical(1))
  mats <- x[keep]
  
  D_cpp <- pairwise_discr_frechet_cpp(mats)
  
  as.dist(D_cpp)
}
