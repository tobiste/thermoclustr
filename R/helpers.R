#' Fit of beta distribution to data
#' 
#' Estimates the shape parameters for a beta distribution using 
#'
#' @param x numeric vector.
#' @param range numeric. two-element vector specifying the interval of the beta distribution
#' @param na.rm a logical value indicating whether `NA` values should be stripped before the computation proceeds.
#'
#' @returns mean, mode, and median (sorted by their values), and the variance 
#' of the beta distribution fitted to `x`.
#' 
#' @export
#' 
#' @references https://en.wikipedia.org/wiki/Beta_distribution
#' 
#' @seealso [beta-stats], [stats::Beta]
#'
#' @examples
#' set.seed(20250411)
#' x <- stats::rbeta(1e5, 3.43, 5.55)
#' fitbeta(x)
fitbeta <- function(x, range = c(0, 1), na.rm = FALSE){
  xm <- (mean(x, na.rm = na.rm) - range[1]) / (range[2] - range[1])
  xv <- var(x, na.rm = na.rm) / (range[2] - range[1])^2
  
  temp1 <- xm * (1 - xm)
  
  if(xv < temp1){
    temp2 <- (temp1 / xv - 1)
    alpha <- xm * temp2
    beta <- (1 - xm) * temp2
  } else {
    alpha <- NA
    beta <- NA
  }
  c(shape1 = alpha, shape2 = beta)
  
}


#' Summary statistics of a beta distribution
#'
#' @param x numeric vector
#' @param shape1,shape2 non-negative parameters of the Beta distribution.
#' @param ... optional arguments passed to [fitbeta()]
#'
#' @returns numeric.
#' @name beta-stats
#' 
#' @seealso [fitbeta()], [stats::Beta]
#' 
#' @references https://en.wikipedia.org/wiki/Beta_distribution
#' 
#' @examples
#' set.seed(20250411)
#' x <- stats::rbeta(1e5, 3.43, 5.55)
#' beta_summary(x)
NULL

#' @rdname beta-stats
#' @export
beta_mean <- function(x, shape1 = NULL, shape2 = NULL){
  if(is.null(shape1) | is.null(shape2)) {
    ab <- fitbeta(x) |>  unname ()
    shape1 <- ab[1]
    shape2 <- ab[2]
  }
  shape1 / (shape1 + shape2)
}

#' @rdname beta-stats
#' @export
beta_var <- function(x, shape1 = NULL, shape2 = NULL){
  if(is.null(shape1) | is.null(shape2)) {
    ab <- fitbeta(x) |> unname()
    shape1 <- ab[1]
    shape2 <- ab[2]
  }
  
  shape1 * shape2 /(
    (shape1 + shape2)^2 * (shape1 + shape2 + 1)
  )
  
}

#' @rdname beta-stats
#' @export
beta_mode <- function(x, shape1 = NULL, shape2 = NULL){
  if(is.null(shape1) | is.null(shape2)) {
    ab <- fitbeta(x) |> unname()
    shape1 <- ab[1]
    shape2 <- ab[2]
  }
  if(shape1>1 & shape2 > 1)  (shape1 - 1) / (shape1 + shape2 - 2) else NA
}

beta_median <- function(x, shape1 = NULL, shape2 = NULL){
  if(is.null(shape1) | is.null(shape2)) {
    ab <- fitbeta(x) |> unname()
    shape1 <- ab[1]
    shape2 <- ab[2]
  }
  if(shape1>1 & shape2 > 1)  (shape1 - 1/3) / (shape1 + shape2 - 2/3) else NA
}

#' @rdname beta-stats
#' @export
beta_summary <- function(x, ...){
  ab <- fitbeta(x, ...) |> unname()
  c(
    Mean = beta_mean(NULL, ab[1], ab[2]),
    Mode = beta_mode(NULL, ab[1], ab[2]),
    Median = beta_median(NULL, ab[1], ab[2])
  ) |> sort() |> 
    append(c(Variance = beta_var(NULL, ab[1], ab[2])))
}


#' Locate elbow in a curve
#' 
#' Useful to locate elbow in elbow curve of within-cluster sum of squares
#'
#' @param x numeric vector
#'
#' @returns integer index of the elbow point
#' @noRd
find_elbow <- function(x) {
  n_points <- length(x)
  all_coords <- cbind(1:n_points, x)
  
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



#' Helper function to set up parallel workflows depending on the OS
#'
#' @param workers The number of parallel processes to use.
#' @export
#' @importFrom future plan
#' @importFrom parallel detectCores
#'
#' @returns NULL
setup_parallel <- function(workers = NULL) {
  if (is.null(workers))
    workers <- max(1, parallel::detectCores() - 1)
  
  if (.Platform$OS.type == "windows") {
    
    future::plan(future::multisession, workers = workers)
    
  } else {
    
    # macOS + Linux
    future::plan(future::multicore, workers = workers)
  }
}

