#' Path statistics
#'
#' Provides binned statistics (mean, median, IQR, quantiles, etc.) of modeled t-T paths.
#'
#' @param x either an object of class `"HeFTy"` (output of [read_hefty()]) or
#' a `data.frame` containing the `time` and `temperature` columns of the modeled
#' paths.
#' @param breaks either a numeric vector of two or more unique cut points or a
#' single number (greater than or equal to 2) giving the number of intervals
#' into which `x` is to be cut.
#'
#' @importFrom dplyr group_by summarise select
#' @importFrom purrr map_dbl
#'
#' @return tibble.
#' @export
#'
#' @examples
#' data(tT_paths1)
#' path_statistics(tT_paths1)
path_statistics <- function(x, breaks = 50) {
  if (inherits(x, "HeFTy")) x <- x$paths

  bins <- time <- temperature <- time_min <- time_median <- time_max <- temp_q <- temp_min <- temp_5 <- temp_median <- temp_95 <- temp_max <- temp_sd <- temp_IQR <-  NULL
  x$bins <- cut(x$time, breaks = breaks)

  dplyr::group_by(x, bins, .add = TRUE) |>
    dplyr::summarise(
      time_min = min(time, na.rm = TRUE),
      time_median = stats::median(time, na.rm = TRUE),
      time_max = max(time, na.rm = TRUE),
      temp_sd = stats::sd(temperature, na.rm = TRUE),
      temp_IQR = stats::IQR(temperature, na.rm = TRUE),
      temp_median = stats::median(temperature, na.rm = TRUE),
      temp_q = list(stats::quantile(temperature, probs = c(0.05, 0.95), na.rm = TRUE)),
      # temp_5 = stats::quantile(temperature, probs = .05),
      # temp_95 = stats::quantile(temperature, probs = .95),
      temp_max = min(temperature, na.rm = TRUE),
      temp_min = max(temperature, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      temp_5 = purrr::map_dbl(temp_q, 1),
      temp_95 = purrr::map_dbl(temp_q, 2)
    ) |>
    dplyr::select(-temp_q, bins, time_min, time_median, time_max, temp_min, temp_5, temp_median, temp_95, temp_max, temp_sd, temp_IQR)
}
