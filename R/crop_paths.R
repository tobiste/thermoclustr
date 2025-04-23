#' Crop cooling paths by time and temperature
#'
#' @param x an object of class `"HeFTy"` (output of [read_hefty()])
#' @param time numeric. two column vector of the desired time range.
#' @param temperature numeric. two column vector of the desired temperature range.
#'
#' @returns object of class `"HeFTy"` (output of [read_hefty()])
#' @export
#'
#' @note Only the `path` slot of `x` will be modified.
#'
#' @examples
#' # example data
#' data(tT_paths1)
#' crop_paths(tT_paths1, time = c(0, 300), temperature = c(0, 200))
crop_paths <- function(x, time = c(0, Inf), temperature = c(0, Inf)) {
  stopifnot(inherits(x, "HeFTy"))
  paths <- dplyr::distinct(x$paths)

  time_range <- sort(time)
  temperature_range <- sort(temperature)

  time_min <- max(time_range[1], min(paths$time))
  time_max <- min(time_range[2], max(paths$time))
  temperature_min <- max(temperature_range[1], min(paths$temperature)) # min temp not smaller than min temp in data
  temperature_max <- min(temperature_range[2], max(paths$temperature)) # max temp not greater than max temp in data

  meta <- paths |> dplyr::select(-time, -temperature)

  suppressWarnings(
    paths_sf <- paths |> 
      sf::st_as_sf(coords = c("time", "temperature")) |>
      dplyr::group_by(segment) |>
      dplyr::summarise(do_union = FALSE) |>
      sf::st_cast("LINESTRING") |>
      dplyr::left_join(meta, dplyr::join_by(segment)) |>
      sf::st_crop(c(xmin = time_min, xmax = time_max, ymin = temperature_min, ymax = temperature_max)) |>
      sf::st_cast("MULTIPOINT") |> 
      sf::st_cast("POINT")
  )

  coords <- sf::st_coordinates(paths_sf) |>
    dplyr::as_tibble()

  x$paths <- sf::st_drop_geometry(paths_sf) |>
    dplyr::bind_cols(coords) |>
    dplyr::select(segment, time = X, temperature = Y, Fit, Comp_GOF) |>
    dplyr::distinct()

  return(x)
}