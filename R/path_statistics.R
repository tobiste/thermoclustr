#' Path statistics
#'
#' Provides binned statistics (mean, median, IQR, quantiles, etc.) of modeled t-T paths.
#'
#' @param x either an object of class `"HeFTy"` (output of [read_hefty()]) or
#' a `data.frame` containing the `time` and `temperature` columns of the modeled
#' paths.
#' @param w numeric vector. Weights for each path segment. 
#' @param breaks either a numeric vector of two or more unique cut points or a
#' single number (greater than or equal to 2) giving the number of intervals
#' into which `x` is to be cut.
#' 
#' @seealso [gof_weighting()] for rescaling goodness-of-fit values to weights.
#'
#' @importFrom dplyr group_by summarise select filter
#' @importFrom purrr map_dbl
#' @importFrom Hmisc wtd.quantile wtd.var
#'
#' @return tibble.
#' @export
#'
#' @examples
#' data(tT_paths)
#' path_statistics(tT_paths, w = gof_weighting(tT_paths$paths$Comp_GOF))
path_statistics <- function(x, w = NULL, breaks = 50) {
  if (inherits(x, "HeFTy")) x <- x$paths

  bins <- time <- temperature <- time_min <- time_median <- time_max <- temp_q <- temp_min <- temp_5 <- temp_median <- temp_95 <- temp_max <- temp_sd <- temp_IQR <- NULL
  x$bins <- cut(x$time, breaks = breaks)
  if(is.null(w)) w <- rep(1, nrow(x))
  x$w <- w

  dplyr::group_by(x, bins, .add = TRUE) |>
    dplyr::filter(!is.na(time), !is.na(temperature)) |> 
    dplyr::summarise(
      time_min = min(time, na.rm = TRUE),
      time_median = Hmisc::wtd.quantile(time, w, probs = 0.5, na.rm = TRUE),
      time_max = max(time, na.rm = TRUE),
      temp_mean = Hmisc::wtd.mean(temperature, w, na.rm = TRUE),
      temp_sd = weighted_sd(temperature, w, na.rm = TRUE),
      temp_IQR = stats::IQR(temperature, na.rm = TRUE),
      temp_median =  Hmisc::wtd.quantile(temperature, w, probs = 0.5, na.rm = TRUE),
      temp_q = list(Hmisc::wtd.quantile(temperature, w, probs = c(0.05, 0.95), na.rm = TRUE)),
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
    dplyr::select(-temp_q, bins, time_min, time_median, time_max, temp_min, temp_5, temp_mean, temp_median, temp_95, temp_max, temp_sd, temp_IQR)
}




weighted_sd <- function(x, w, ...){
  sqrt(Hmisc::wtd.var(x, w, ...))
}

