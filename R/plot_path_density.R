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
NULL

#' @rdname plt_density
#' @export
#' @importFrom ggplot2 aes geom_density2d_filled ggplot
plot_path_density_filled <- function(x, bins = 50L, densify = TRUE, show.legend = NA, ...) {
  time <- temperature <- NULL
  if (densify) x <- densify_paths(x, ...)

  ggplot(data = x, aes(x = time, y = temperature)) +
    geom_density2d_filled(contour_var = "ndensity", bins = bins, show.legend = show.legend)
}

#' @rdname plt_density
#' @export
#' @importFrom ggplot2 aes geom_density2d ggplot
plot_path_density <- function(x, bins = 50L, densify = TRUE, show.legend = NA, ...) {
  time <- temperature <- NULL
  if (densify) x <- densify_paths(x, ...)

  ggplot(data = x, aes(x = time, y = temperature)) +
    geom_density2d(contour_var = "ndensity", bins = bins, show.legend = show.legend)
}
