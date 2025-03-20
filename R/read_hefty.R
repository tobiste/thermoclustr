combine_rows <- function(x) {
  xm <- as.matrix(x)

  odd_rows <- xm[seq_along(xm[, 1]) %% 2 != 0, ] # all odd rows
  even_rows <- xm[seq_along(xm[, 1]) %% 2 == 0, ] # all even rows

  time.col <- as.vector(t(odd_rows))
  temp.col <- as.vector(t(even_rows))


  data.frame("time" = as.numeric(time.col), "temperature" = as.numeric(temp.col))
}

assign_segment <- function(x) { # create new function called assign_segment
  n <- length(x)
  segment <- integer(n)

  segment[1] <- 1L # create a vector with first value of 1
  for (i in 2:n) { # start with index 2 of vector (second position in vector)
    if (x[i] < x[i - 1]) { # if the value of x (time) in position i is less than the value of x in the previous position (i -1),
      segment[i] <- segment[i - 1] # then give that position the same segment value as position before it
    } else { # otherwise
      segment[i] <- segment[i - 1] + 1L # take the segment value from the position before it, add one, and give that position the new segment value
    }
  }
  return(forcats::as_factor(segment)) # return the output of segment as a factor (this is easier for ggplot according to Tobi)
}




#' Read time, temperature and GOF data from HeFTy output from excel file
#'
#'
#' @description
#' `r lifecycle::badge("superseded")`
#'
#'  `read_hefty_xlsx()` has been superseded in favor of [read_hefty()] to allow
#' direct import from HeFTy output files with out additional format
#' manipulation. The function `read_hefty_xlsx()`will be removed in a future
#' release.
#'
#' @param fname path to the excel spreadsheet that contains the HeFTy outputs, i.e. the t-T-paths in sheet 1 and the GOF values in sheet 2
#'
#' @return `data.frame` of the combined data
#'
#' @importFrom readxl read_xlsx
#' @importFrom forcats fct fct_reorder as_factor
#' @importFrom dplyr mutate as_tibble between case_when
#'
#' @seealso [read_hefty()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' path2myfile <- system.file("s14MM_v1.xlsx", package = "HeFTy.SmoothR")
#' read_hefty_xlsx(path2myfile)
#' }
read_hefty_xlsx <- function(fname) {
  Fit <- time <- NULL

  x <- readxl::read_xlsx(fname, sheet = 1, col_names = FALSE) # load t-T data
  GOF <- readxl::read_xlsx(fname, sheet = 2, col_names = TRUE) # load GOF sheet

  hs.input <- dplyr::select(x, -1) |>  # select and remove the '...1' from the column names
    combine_rows() |>  # run combine_rows function
    dplyr::mutate(segment = assign_segment(time)) # run Assign_segment function on the time column of hs.input; then mutate function takes output of assign_segment

  merge(hs.input, GOF, by = "segment") |>
    dplyr::as_tibble() |>
    dplyr::mutate(
      Fit = dplyr::case_when(
        Comp_GOF >= 0.9 ~ "Best",
        dplyr::between(Comp_GOF, 0.5, 0.9) ~ "Good",
        dplyr::between(Comp_GOF, 0.05, .5) ~ "Acceptable",
        .default = NA
      ),
      Fit = forcats::fct(Fit, levels = c(NA, "Acceptable", "Good", "Best"))
    ) |>
    dplyr::arrange(Fit, Comp_GOF) |>
    dplyr::mutate(segment = forcats::fct_reorder(segment, Comp_GOF))
}

#' Import inverse thermal history model
#'
#' Read time, temperature and GOF data of modeled thermal history paths out of
#' HeFTy output .txt file
#'
#' @param fname path to the .txt file that contains the HeFTy outputs
#'
#' @return object of class `"HeFTy"`, i.e. `list` of the individual paths, 
#' the constraints, the weighted mean path and the grain summary statistics.
#'
#' @importFrom dplyr across as_tibble everything matches rename
#'
#' @export
#' @examples
#' path2myfile <- system.file("112-9-30-zr-inv.txt", package = "thermochron")
#' read_hefty(path2myfile)
read_hefty <- function(fname) {
  Fit <- value <- Comp_GOF <- time <- segment <- temperature <- V1 <- V2 <- V3 <- V4 <- V5 <- constraint <- NULL
  file <- readLines(fname) |>
    strsplit("\t", useBytes = TRUE)

  # extract individual paths
  individual_paths_loc <- grep("Individual paths", file)
  individual_paths_mat <-  do.call(rbind, file[(individual_paths_loc + 4):length(file)])

  odd <- seq_along(individual_paths_mat[, 1]) %% 2 == 0

  GOF <- individual_paths_mat[odd, 2] |>
    as_tibble() |>
    rename(Comp_GOF = value) |>
    mutate(
      Comp_GOF = as.numeric(Comp_GOF),
      segment = forcats::as_factor(dplyr::row_number()),
      Fit = dplyr::case_when(
        Comp_GOF >= 0.9 ~ "Best",
        dplyr::between(Comp_GOF, 0.5, 0.9) ~ "Good",
        dplyr::between(Comp_GOF, 0.05, .5) ~ "Acceptable",
        .default = NA
      ),
      Fit = forcats::fct(Fit, levels = c(NA, "Acceptable", "Good", "Best"))
    )

  time_loc <- grep('Time', individual_paths_mat[1, ]) + 1
  paths <- individual_paths_mat[, time_loc:ncol(individual_paths_mat)] |>
    as_tibble() |>
    combine_rows() |>
    # dplyr::filter(!is.na(time) & !is.na(temperature)) |>
    mutate(segment = assign_segment(time)) |>
    dplyr::right_join(GOF, dplyr::join_by(segment)) |>
    dplyr::select(segment, time, temperature, Fit, Comp_GOF) |>
    dplyr::arrange(Fit, Comp_GOF) |>
    dplyr::mutate(segment = forcats::fct_reorder(segment, Comp_GOF)) |>
    as_tibble()

  # extract weighted mean path
  wm_path_loc <- grep("Weighted mean path", file)
  wm <- do.call(rbind, file[wm_path_loc + c(1, 2)] ) |>
    t() |>
    utils::tail(-1) |>
    as_tibble() |>
    rename(time = V1, temperature = V2) |>
    mutate(across(everything(), as.numeric))

  # extract constraints
  inversion_terminated_loc <- grep("Inversion", file)
  constr1 <- do.call(rbind, file[3:inversion_terminated_loc - 1])

  constraints <- constr1[, 1:5] |>
    utils::tail(-1) |>
    as_tibble() |>
    rename("constraint" = V1, "max_time" = V2, "min_time" = V3, "max_temp" = V4, "min_temp" = V5) |>
    mutate(across(everything(), as.numeric),
      constraint = forcats::as_factor(constraint)
    )

  summaries_loc <- grep("Summaries", file)
  grain_summary <-  do.call(rbind, file[summaries_loc:(wm_path_loc - 1)]) |>
    t() |>
    utils::tail(-1) |>
    as_tibble() |>
    rename("grain" = V1, "mean" = V2, sd = V3, min = V4, max = V5) |>
    mutate(across(!matches("grain"), as.numeric))

  res <- list(
    paths = paths,
    constraints = constraints,
    weighted_mean_path = wm,
    summary = grain_summary
  )
  class(res) <- append(class(res), 'HeFTy')
  return(res)
}
