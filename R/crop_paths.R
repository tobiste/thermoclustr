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
  stopifnot(inherits(x, "HeFTy") || inherits(x, "data.frame"))
  segment <- geometry <- Fit <- Comp_GOF <- NULL
  # if(inherits(x, "HeFTy")){
  #   paths <- dplyr::distinct(x$paths)
  # } else {
  #   paths <- x
  # }
  # Extract paths
  paths <- if (inherits(x, "HeFTy")) x$paths else x

  time_range <- sort(time)
  temperature_range <- sort(temperature)

  time_range[1] <- max(time_range[1], min(paths$time))
  time_range[2] <- min(time_range[2], max(paths$time))
  temperature_range[1] <- max(temperature_range[1], min(paths$temperature)) # min temp not smaller than min temp in data
  temperature_range[2] <- min(temperature_range[2], max(paths$temperature)) # max temp not greater than max temp in data
  #
  # meta <- paths |> dplyr::select(-time, -temperature)
  #
  # suppressWarnings(
  #   paths_sf <- paths |>
  #     sf::st_as_sf(coords = c("time", "temperature")) |>
  #     dplyr::group_by(segment) |>
  #     dplyr::summarise(do_union = FALSE) |>
  #     sf::st_cast("LINESTRING") |>
  #     dplyr::left_join(meta, dplyr::join_by(segment)) |>
  #     sf::st_crop(c(xmin = time_min, xmax = time_max, ymin = temperature_min, ymax = temperature_max)) |>
  #     sf::st_cast("MULTIPOINT") |>
  #     sf::st_cast("POINT")
  # )
  #
  # coords <- sf::st_coordinates(paths_sf) |>
  #   dplyr::as_tibble()
  #
  # paths_cropped <- sf::st_drop_geometry(paths_sf) |>
  #   dplyr::bind_cols(coords) |>
  #   dplyr::select(segment, time = X, temperature = Y, Fit, Comp_GOF) |>
  #   dplyr::distinct()
  #
  # if(inherits(x, "HeFTy")) {
  #   res <- x
  #   res$paths <- paths_cropped
  # } else {
  #   res <- paths_cropped
  # }
  # return(res)
  crop_box <- c(
    xmin = time_range[1], xmax = time_range[2],
    ymin = temperature_range[1], ymax = temperature_range[2]
  )

  suppressWarnings(suppressMessages({
    sf::sf_use_s2(FALSE)

    # Convert to LINESTRING per segment
    paths_sf <- sf::st_as_sf(paths, coords = c("time", "temperature")) |>
      dplyr::group_by(segment) |>
      dplyr::summarise(geometry = sf::st_combine(geometry), .groups = "drop") |>
      sf::st_cast("LINESTRING")

    # Crop lines (this splits them at box boundaries!)
    cropped <- sf::st_crop(paths_sf, crop_box)

    # Cast to individual points again (adds intersection vertices!)
    cropped_points <- cropped |>
      sf::st_cast("MULTIPOINT", warn = FALSE) |>
      sf::st_cast("POINT", warn = FALSE)

    # Add attributes back (Fit and Comp_GOF from original)
    meta <- dplyr::distinct(paths[, c("segment", "Fit", "Comp_GOF")])

    result <- cropped_points |>
      dplyr::left_join(meta, by = "segment")

    # Extract coordinates back to columns
    coords <- sf::st_coordinates(result)
    out <- sf::st_drop_geometry(result)
    out$time <- coords[, "X"]
    out$temperature <- coords[, "Y"]

    sf::sf_use_s2(TRUE)
  }))


  out <- dplyr::select(out, segment, time, temperature, Fit, Comp_GOF) |>
    dplyr::distinct()

  if (inherits(x, "HeFTy")) {
    x$paths <- out
    return(x)
  } else {
    return(out)
  }
}
