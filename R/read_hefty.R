combine_rows <- function(x) {
  xm <- as.matrix(x)

  odd_rows <- xm[seq_along(xm[, 1]) %% 2 != 0, ] # all odd rows
  even_rows <- xm[seq_along(xm[, 1]) %% 2 == 0, ] # all even rows

  time.col <- as.vector(t(odd_rows))
  temp.col <- as.vector(t(even_rows))
  
  # n_pairs <- nrow(xm) / 2
  # 
  # time.col <- matrix(xm[seq(1, nrow(xm), 2), ], ncol = ncol(xm), byrow = FALSE)
  # temp.col <- matrix(xm[seq(2, nrow(xm), 2), ], ncol = ncol(xm), byrow = FALSE)
  # 
  data.frame(
    time = as.numeric(time.col),
    temperature = as.numeric(temp.col),
    stringsAsFactors = FALSE
  )
}

assign_segment <- function(x) { # create new function called assign_segment
  # n <- length(x)
  # segment <- integer(n)
  # 
  # segment[1] <- 1L # create a vector with first value of 1
  # for (i in 2:n) { # start with index 2 of vector (second position in vector)
  #   if (x[i] < x[i - 1]) { # if the value of x (time) in position i is less than the value of x in the previous position (i -1),
  #     segment[i] <- segment[i - 1] # then give that position the same segment value as position before it
  #   } else { # otherwise
  #     segment[i] <- segment[i - 1] + 1L # take the segment value from the position before it, add one, and give that position the new segment value
  #   }
  # }
  # return(forcats::as_factor(segment)) # return the output of segment as a factor (this is easier for ggplot according to Tobi)
  
  # Find where x increases or stays the same (mark with 1), where it decreases mark 0
  increases <- c(TRUE, x[-1] >= x[-length(x)])
  
  # Cumulative sum of increases gives segment number
  segment <- cumsum(increases)
  
  # Return factor
  as.factor(segment)
}

#' Parse the hefty txt files
#'
#' @param fname path to the .txt file that contains the HeFTy outputs
#' @importFrom data.table fread
#' 
#' @returns list of character vectors, where each vector corresponds to a line in the file.
#' @noRd
parse_hefty <- function(fname){
  f <- data.table::fread(fname, sep = "\n", header = FALSE, skip = 0, colClasses = 'character', blank.lines.skip	= TRUE, col.names = "") |>
    lapply(function(x) {
      strsplit(x, "\t")
    })
  f[[1]]
}



#' Read time, temperature and GOF data from HeFTy output from excel file
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
#' @examples
#' \dontrun{
#' path2myfile <- system.file("s14MM_v1.xlsx", package = "HeFTy.SmoothR")
#' read_hefty_xlsx(path2myfile)
#' }
read_hefty_xlsx <- function(fname) {
  Fit <- time <- NULL
  
  x <- readxl::read_xlsx(fname, sheet = 1, col_names = FALSE) # load t-T data
  GOF <- readxl::read_xlsx(fname, sheet = 2, col_names = TRUE) # load GOF sheet
  
  hs.input <- dplyr::select(x, -1) |> # select and remove the '...1' from the column names
    combine_rows() |> # run combine_rows function
    dplyr::mutate(segment = assign_segment(time)) # run Assign_segment function on the time column of hs.input; then mutate function takes output of assign_segment
  
  merge(hs.input, GOF, by = "segment") |>
    # dplyr::as_tibble() |>
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
  # file <- readLines(fname) |>
  #   strsplit("\t", useBytes = TRUE)
  
  file <- parse_hefty(fname)
  
  # extract individual paths
  individual_paths_loc <- grep("Individual paths", file)
  # individual_paths_mat <- do.call(rbind, file[(individual_paths_loc + 4):length(file)])
  
  # repeated rbind is slow:
  individual_data <- file[(individual_paths_loc + 4):length(file)]
  individual_paths_mat <- matrix(
    unlist(individual_data, use.names = FALSE),
    ncol = length(individual_data[[1]]),
    byrow = TRUE
  )
  
  # odd <- seq_along(individual_paths_mat[, 1]) %% 2 == 0
  # GOF <- individual_paths_mat[odd, 2] |>
  #   as_tibble() |>
  #   rename(Comp_GOF = value) |>
  #   mutate(
  #     Comp_GOF = as.numeric(Comp_GOF),
  #     segment = as.factor(dplyr::row_number()),
  #     Fit = dplyr::case_when(
  #       Comp_GOF >= 0.9 ~ "Best",
  #       dplyr::between(Comp_GOF, 0.5, 0.9) ~ "Good",
  #       dplyr::between(Comp_GOF, 0.05, .5) ~ "Acceptable",
  #       .default = NA
  #     ),
  #     Fit = factor(Fit, levels = c(NA, "Acceptable", "Good", "Best"))
  #   )
  odd_idx <- seq(2, nrow(individual_paths_mat), by = 2)
  Comp_GOF <- as.numeric(individual_paths_mat[odd_idx, 2])
  GOF <- data.frame(
    Comp_GOF = Comp_GOF,
    segment = as.factor(seq_along(Comp_GOF)),
    Fit = dplyr::case_when(
      Comp_GOF >= 0.9 ~ "Best",
      Comp_GOF >= 0.5 ~ "Good",
      Comp_GOF >= 0.05 ~ "Acceptable",
      TRUE ~ NA_character_
    )
  )|>
    dplyr::mutate(Fit = factor(Fit, levels = c(NA, "Acceptable", "Good", "Best"))) #|> 
  # dplyr::arrange(Fit, Comp_GOF) 
  
  time_loc <- grep("Time", individual_paths_mat[1, ]) + 1
  # paths <- individual_paths_mat[, time_loc:ncol(individual_paths_mat)] |>
  #   as_tibble() |>
  #   combine_rows() |>
  #   # dplyr::filter(!is.na(time) & !is.na(temperature)) |>
  #   mutate(segment = assign_segment(time)) |>
  #   dplyr::right_join(GOF, by = "segment") |>
  #   dplyr::select(segment, time, temperature, Fit, Comp_GOF) |>
  #   dplyr::arrange(Fit, Comp_GOF) |>
  #   dplyr::mutate(segment = forcats::fct_reorder(segment, Comp_GOF)) |>
  #   as_tibble()
  #faster:
  path_data <- individual_paths_mat[, time_loc:ncol(individual_paths_mat)] 
  paths <- combine_rows(path_data)|>
    dplyr::mutate(
      segment = assign_segment(time)
    )|>
    dplyr::right_join(GOF, by = "segment")|>
    dplyr::select(segment, time, temperature, Fit, Comp_GOF)|>
    dplyr::arrange(Fit, Comp_GOF)|>
    dplyr::mutate(segment = forcats::fct_reorder(segment, Comp_GOF))
  
  # extract weighted mean path
  wm_path_loc <- grep("Weighted mean path", file)
  # wm <- do.call(rbind, file[wm_path_loc + c(1, 2)]) |>
  #   t() |>
  #   utils::tail(-1) |>
  #   as_tibble() |>
  #   rename(time = V1, temperature = V2) |>
  #   mutate(across(everything(), as.numeric))
  wm_data <- file[wm_path_loc + c(1, 2)]
  wm0 <- matrix(unlist(wm_data, use.names = FALSE), ncol = 2, byrow = F)
  wm_length <- nrow(wm0)
  wm <- data.frame(
    time = as.numeric(wm0[2:wm_length, 1]),
    temperature = as.numeric(wm0[2:wm_length, 2])
  )
  
  # extract constraints
  inversion_terminated_loc <- grep("Inversion", file)
  # constr1 <- do.call(rbind, file[3:inversion_terminated_loc - 1])
  # constraints <- constr1[, 1:5] |>
  #   utils::tail(-1) |>
  #   as_tibble() |>
  #   rename("constraint" = V1, "max_time" = V2, "min_time" = V3, "max_temp" = V4, "min_temp" = V5) |>
  #   mutate(across(everything(), as.numeric),
  #     constraint = forcats::as_factor(constraint)
  #   )
  # Faster:
  constraint_data <- file[3:(inversion_terminated_loc - 1)]
  constraint_mat <- matrix(unlist(constraint_data, use.names = FALSE), nrow = 5, byrow = TRUE)[, 1:5]
  constraints <- data.frame(
    constraint = as.factor(constraint_mat[, 1]),
    max_time = as.numeric(constraint_mat[, 2]),
    min_time = as.numeric(constraint_mat[, 3]),
    max_temp = as.numeric(constraint_mat[, 4]),
    min_temp = as.numeric(constraint_mat[, 5])
  )
  
  summaries_loc <- grep("Summaries", file)
  # grain_summary <- do.call(rbind, file[summaries_loc:(wm_path_loc - 1)]) |>
  #   t() |>
  #   utils::tail(-1) |>
  #   as_tibble() |>
  #   rename("grain" = V1, "mean" = V2, sd = V3, min = V4, max = V5) |>
  #   mutate(across(!matches("grain"), as.numeric))
  #   FASTER:
  summary_data <- file[summaries_loc:(wm_path_loc - 1)]
  summary_mat <- matrix(unlist(summary_data, use.names = FALSE), ncol = 5, byrow = FALSE)[-1, ]
  grain_summary <- data.frame(
    grain = summary_mat[, 1],
    mean = as.numeric(summary_mat[, 2]),
    sd = as.numeric(summary_mat[, 3]),
    min = as.numeric(summary_mat[, 4]),
    max = as.numeric(summary_mat[, 5])
  )
  
  res <- list(
    paths = paths,
    constraints = constraints,
    weighted_mean_path = wm,
    summary = grain_summary
  )
  class(res) <- append(class(res), "HeFTy")
  return(res)
}





#' Legacy version of read_hefty
#' @noRd
read_hefty_old <- function(fname) {
  Fit <- value <- Comp_GOF <- time <- segment <- temperature <- V1 <- V2 <- V3 <- V4 <- V5 <- constraint <- NULL
  file <- readLines(fname) |>
    strsplit("\t", useBytes = TRUE)

  # extract individual paths
  individual_paths_loc <- grep("Individual paths", file)
  individual_paths_mat <- do.call(rbind, file[(individual_paths_loc + 4):length(file)])

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

  time_loc <- grep("Time", individual_paths_mat[1, ]) + 1
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
  wm <- do.call(rbind, file[wm_path_loc + c(1, 2)]) |>
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
  grain_summary <- do.call(rbind, file[summaries_loc:(wm_path_loc - 1)]) |>
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
  class(res) <- append(class(res), "HeFTy")
  return(res)
}
