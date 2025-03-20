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
#' @return tibble.
#' @export
#'
#' @examples
#' data(tT_paths1)
#' path_statistics(tT_paths1)
path_statistics <- function(x, breaks = 50){
  if(inherits(x, 'HeFTy')) x <- x$paths
  
  bins <- time <- temperature <- NULL
  dplyr::mutate(x,
         bins = cut(time, breaks = breaks)) |>
    dplyr::group_by(bins, .add = TRUE) |>
    dplyr::summarise(
      time_min = min(time),
      time_median = stats::median(time, na.rm = TRUE),
      time_max = max(time),
      temp_sd = stats::sd(temperature, na.rm = TRUE),
      temp_IQR = stats::IQR(temperature, na.rm = TRUE),
      temp_median = stats::median(temperature),
      temp_5 = stats::quantile(temperature, probs = .05),
      temp_95 = stats::quantile(temperature, probs = .95),
      temp_max = min(temperature),
      temp_min = max(temperature),
      .groups = 'drop'
    )
}
