#' Hausdorff Distance
#'
#' Hausdorff distance (aka Hausdorff dimension)
#'
#' @param P,Q numerical matrices, representing points in an m-dim. space.
#' @param check logical. Checks if input arguments are valid. 
#'
#' @returns A single scalar, the Hausdorff distance (dimension).
#'
#' @name hausdorff
#'
#' @examples
#' P <- matrix(c(1, 1, 2, 2, 5, 4, 5, 4), 4, 2)
#' Q <- matrix(c(4, 4, 5, 5, 2, 1, 2, 1), 4, 2)
#' hausdorff(P, Q) # 4.242641 = sqrt(sum((c(4,2)-c(1,5))^2))
#' hausdorff_fast(P, Q)

#' @rdname hausdorff
#' @importFrom matrixStats rowMins colMins
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
#' hausdorff_list_fast(tT_paths_list)
#' hausdorff_list_C(tT_paths_list)
#' hausdorff_list_parallel(tT_paths_list)
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

#' @rdname hausdorff_dmat
#' @export
hausdorff_list_fast <- function(x) {
  keep <- vapply(x, function(m)
    !is.null(m) && length(m) > 0L,
    logical(1))
  
  mats <- x[keep]
  n <- length(mats)
  
  out <- matrix(0, n, n)
  
  for (i in seq_len(n)) {
    Pi <- mats[[i]]
    
    for (j in seq_len(i - 1)) {
      out[i, j] <- hausdorff_fast(Pi, mats[[j]])
    }
  }
  
  as.dist(out)
}


# uses Cpp functions
#' @rdname hausdorff_dmat
#' @export
hausdorff_list_C <- function(x) {
  # This code is faster than the R version. However, the code is still overall 
  # way slower than the sf version. The reason is that it still loops through 
  # 1000's of matrices with each 1000's of rows...
  keep <- vapply(x, function(m) !is.null(m) && length(m) > 0L, logical(1))
  mats <- x[keep]
  
  D_cpp <- pairwise_hausdorff_cpp(mats)
  
  as.dist(D_cpp)
}

# My algorithm builds an n × m dense matrix. 
# If paths have 300 points each 1000 paths you compute ~500,000 Hausdorff calls × 90,000 distances each. That explodes quickly.
# But Hausdorff distance only needs max(min(d(p,q))) Thus you never need the full matrix.

# uses parallel calculations, no distance matrix, O(n log n) instead of O(nm), uses RANN (very fast, C++ kd-tree)
#' @rdname hausdorff_dmat
#' @export future.apply future_sapply
hausdorff_list_parallel2 <- function(x, ...) {
  # setup_parallel(...)
  
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
      hausdorff_fast(mats[[i]], mats[[j]])
    }
  )
  
  out <- matrix(0, n, n)
  out[lower.tri(out)] <- d
  
  as.dist(out)
  # plan(sequential)
}

#' Frechet Distance
#'
#' Frechet distance (aka Frechet dimension)
#'
#' @param P,Q numerical matrices, representing points in an m-dim. space.
#' @param check logical. Checks if input arguments are valid. 
#'
#' @returns A single scalar, the Hausdorff distance (dimension).
#'
#' @name frechet
#'
#' @examples
#' P <- matrix(c(1, 1, 2, 2, 5, 4, 5, 4), 4, 2)
#' Q <- matrix(c(4, 4, 5, 5, 2, 1, 2, 1), 4, 2)
#' frechet(P, Q) # 4.242641 = sqrt(sum((c(4,2)-c(1,5))^2))
#' frechet_discrete_fast(P, Q)

#' @rdname frechet
#' @importFrom SimilarityMeasures Frechet
frechet <- function(P, Q, check = TRUE){
  # This code is implements the Frechet algorithm from the SimilarityMeasures package. 
  
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
  
  SimilarityMeasures::Frechet(P, Q)
}

#' @rdname frechet
frechet_discrete_fast <- function(P, Q) {
  n <- nrow(P)
  m <- nrow(Q)
  
  ca <- matrix(-1, n, m)
  
  d <- function(i, j) {
    sqrt(sum((P[i,] - Q[j,])^2))
  }
  
  rec <- function(i, j) {
    
    if (ca[i,j] > -1)
      return(ca[i,j])
    
    if (i == 1 && j == 1)
      ca[i,j] <<- d(1,1)
    
    else if (i > 1 && j == 1)
      ca[i,j] <<- max(rec(i-1,1), d(i,1))
    
    else if (i == 1 && j > 1)
      ca[i,j] <<- max(rec(1,j-1), d(1,j))
    
    else
      ca[i,j] <<- max(
        min(
          rec(i-1,j),
          rec(i-1,j-1),
          rec(i,j-1)
        ),
        d(i,j)
      )
    
    ca[i,j]
  }
  
  rec(n,m)
}



#' Fréchet distance matrix
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
#' frechet_list_R(tT_paths_list)
#' cont_frechet_list_C(tT_paths_list)
#' disc_frechet_list_C(tT_paths_list)
#' 
#' setup_parallel()
#' frechet_list_parallel(tT_paths_list)
#' future::plan(sequential)
NULL

#' @rdname frechet_dmat
#' @export
frechet_list_R <- function(x) {
  # This code implements the discrete Frechet algorithm. 
  # The code is 10x faster than the normal version
  
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
      out[i, j] <- frechet_discrete_fast(Pi, mats[[j]])
    }
  }
  
  as.dist(out)
}

#' @rdname frechet_dmat
#' @export
cont_frechet_list_C <- function(x) {
  # The code is still overall 
  # way slower than the sf version. The reason is that it still loops through 
  # 1000's of matrices with each 1000's of rows...
  keep <- vapply(x, function(m) !is.null(m) && length(m) > 0L, logical(1))
  mats <- x[keep]
  
  D_cpp <- pairwise_cont_frechet_cpp(mats)
  
  as.dist(D_cpp)
}

#' @rdname frechet_dmat
#' @export
disc_frechet_list_C <- function(x) {
  # The code is still overall 
  # way slower than the sf version. The reason is that it still loops through 
  # 1000's of matrices with each 1000's of rows...
  keep <- vapply(x, function(m) !is.null(m) && length(m) > 0L, logical(1))
  mats <- x[keep]
  
  D_cpp <- pairwise_discr_frechet_cpp(mats)
  
  as.dist(D_cpp)
}


# uses parallel calculations 
#' @rdname hausdorff_dmat
#' @export future.apply future_sapply
#' @importFrom package function
frechet_list_parallel <- function(x, ...) {
  # setup_parallel(...)
  
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
      frechet_discrete_fast(mats[[i]], mats[[j]])
    }
  )
  
  out <- matrix(0, n, n)
  out[lower.tri(out)] <- d
  
  as.dist(out)
  # plan(sequential)
}
