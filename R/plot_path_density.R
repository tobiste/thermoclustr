#' Path Density Plot
#'
#' Creates a 2d kernel density estimate of the t-T paths and plots it using
#' ggplot
#'
#' @param x either an object of class `"HeFTy"` (output of [read_hefty()]),
#' a `data.frame` containing the `time`, `temperature` columns of the modeled
#' paths, or output of [densify_paths()].
#' @param bins integer. Amount of bins used for the kernel density estimate. 50 by default.
#' @param densify logical. Whether the paths in `x` should be densified first?
#' Default is `TRUE`.
#' @param show.legend logical. Should this layer be included in the legends? `NA`,
#' the default, includes if any aesthetics are mapped. `FALSE` never includes,
#' and `TRUE` always includes. It can also be a named logical vector to finely
#' select the aesthetics to display.
#' @param geom Use to override the default connection between
#' [ggplot2::geom_density_2d()] and [ggplot2::stat_density_2d()].
#' For more information at overriding these connections, see how the stat and geom arguments work.
#' @param n integer. Number of grid points in each direction.
#' @param weights numeric vector. Weights for each path segment. 
#' @param ... Arguments passed on to [densify_paths()] (only if `densify=TRUE`).
#'
#' @return ggplot
#' @import ggplot2
#'
#' @name plt_density
#'
#' @examples
#' data(tT_paths1)
#' plot_path_density(tT_paths1)
#' plot_path_density_filled(tT_paths1)
#' plot_path_density_filled(tT_paths1, geom = "raster")
#' plot_path_density_filled_weighted(tT_paths1, weights = tT_paths1$paths$Comp_GOF)
NULL

#' @rdname plt_density
#' @export
#' @importFrom ggplot2 aes geom_density2d_filled ggplot
plot_path_density_filled <- function(x, bins = 50L, densify = TRUE,
                                     show.legend = NA, geom = "density_2d_filled",
                                     n = 100L, ...) {
  time <- temperature <- ndensity <- NULL
  if (densify) x <- densify_paths(x, ...)

  if (geom == "density_2d_filled") {
    ggplot(data = x, aes(x = time, y = temperature)) +
      geom_density_2d_filled(
        contour_var = "ndensity",
        n = n,
        bins = bins, show.legend = show.legend
      )
  } else {
    ggplot(data = x, aes(x = time, y = temperature)) +
      stat_density_2d_filled(
        geom = geom,
        # contour_var = "ndensity",
        aes(fill = ggplot2::after_stat(ndensity)),
        contour = FALSE,
        n = n,
        bins = bins, show.legend = show.legend
      )
  }
}

#' @rdname plt_density
#' @export
#' @importFrom ggplot2 aes geom_density2d ggplot
plot_path_density <- function(x, bins = 50L, densify = TRUE, show.legend = NA, n = 100L, ...) {
  time <- temperature <- NULL
  if (densify) x <- densify_paths(x, ...)

  ggplot(data = x, aes(x = time, y = temperature)) +
    geom_density2d(contour_var = "ndensity", bins = bins, show.legend = show.legend, n = n)
}


kde2d.weighted <- function(x, y, w, h, n = 25, lims = c(range(x), range(y))) {
  nx <- length(x)
  if (length(y) != nx) stop("data vectors must be the same length")
  if (length(w) != nx & length(w) != 1) stop("weight vectors must be 1 or length of data")
  gx <- seq(lims[1], lims[2], length = n) # gridpoints x
  gy <- seq(lims[3], lims[4], length = n) # gridpoints y
  if (missing(h)) h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))
  if (missing(w)) w <- numeric(nx) + 1
  h <- h / 4
  ax <- outer(gx, x, "-") / h[1] # distance of each point to each grid point in x-direction
  ay <- outer(gy, y, "-") / h[2] # distance of each point to each grid point in y-direction
  z <- (matrix(rep(w, n), nrow = n, ncol = nx, byrow = TRUE) * matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx)) / (sum(w) * h[1] * h[2]) # z is the density
  return(list(x = gx, y = gy, z = z))
}


#' @rdname plt_density
#' @export
#' @importFrom ggplot2 aes geom_tile ggplot
#' @importFrom dplyr left_join bind_cols
plot_path_density_filled_weighted <- function(x, bins = 50L, densify = TRUE,
                                              show.legend = NA, #geom = "density_2d_filled",
                                              n = 100L, weights = NULL, ...) {
  time <- temperature <- ndensity <- NULL

  segments <- x$paths$segment
  if (is.null(weights)) {
    w <- rep(1, length(segments))
  } else {
    w <- x$paths$Comp_GOF
  }

  weights <- data.frame(segment = segments, w = w) |>
    unique()

  if (isTRUE(densify)) {
    x <- densify_paths(x, ...) |>
      left_join(weights, by = "segment")
  }

  kde_results <- kde2d.weighted(x$time, x$temperature, w = x$w, n = n, h = rep(bins, 2))

  kde_data <- expand.grid(x = kde_results$x, y = kde_results$y)
  kde_data$z <- as.vector(kde_results$z)

  #if (geom == "density_2d_filled") {
  ggplot(kde_data, aes(x = x, y = y, fill = z)) +
    geom_tile(show.legend = show.legend)
  # } else {
  #   interpdf <- interp::interp2xyz(interp::interp(x=kde_data$x, y=kde_data$y, z=kde_data$z, duplicate="mean"), data.frame=TRUE)
  #   ggplot(interpdf, aes(x = x, y = y, z = z, fill = z)) +
  #    # geom_tile() +
  #    geom_contour() 
  # }
}
