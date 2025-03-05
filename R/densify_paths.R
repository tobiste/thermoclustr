#' Densify Paths
#'
#' Adds extra points between the vertices of the path segments
#'
#' @param x t-T and GOF data of the modeled paths. Output of [read_hefty()].
#' @param GOF_rank numeric. Selects only the `GOF_rank`-th highest GOF ranked paths.
#' If all GOFs should be used, set to `Inf`. Default is `10`.
#' @param n integer. Adds `n` (10 by
#' default) equally-spaced extra points along each path
#' segment (between vertices).
#' @param max_distance numeric. Adds points at a maximum distance of
#' `max_distance` (in Myr) from each other. 1 by default.
#' @param samples integer or character. Number of random samples of the data.
#' This number should be less or equal then the amount of paths.
#' The default is `100`. Paths will be randomly selected **after** the data has
#' been filtered by the `GOF_rank()`.
#' Optional, set `samples` to `'all'` to consider all paths ignoring the
#' `GOF_rank()` filter (this sets `GOF` to `Inf`).
#' Set `samples` to `GOF` to consider all paths after the `GOF_rank` filter.
#' @param replace logical. Should sampling be with replacement?
#'
#' @note A large sample number `n` will require a **long(!)**
#' processing time for this function and subsequent methods such as
#' [plot_path_density()] or [cluster_paths()].
#'
#' If only paths within a specified GOF range should be densified, create a
#' subset of the data beforehand using either [subset()] or [dplyr::filter()].
#'
#' @return tibble
#'
#' @importFrom sf st_as_sf st_cast st_drop_geometry st_coordinates
#' @importFrom smoothr densify
#' @importFrom dplyr dense_rank filter distinct sample_n semi_join mutate summarise bind_cols rename select
#'
#' @export
#'
#' @examples
#' data(tT_paths1)
#' densify_paths(tT_paths1)
densify_paths <- function(x, GOF_rank = 10L, n = 10L, max_distance = 1, samples = 100L, replace = FALSE) {
  L1 <- L2 <- X <- Y <- numeric()
  segment <- time <- temperature <- NULL

  stopifnot(is.numeric(samples) | samples %in% c("all", "GOF"))
  if (samples == "all") {
    GOF_rank <- Inf
    samples <- "GOF"
  }


  # Subset Data by Highest N GOF Values (If Desired)
  x$rank <- dplyr::dense_rank(-x$Comp_GOF) # Rank Comp_GOF values with lowest rank (1) being highest GOF value
  hs.input <- dplyr::filter(x, rank <= GOF_rank) # Remove all t-T points with rank greater than N; new table shows top N GOF t-T points

  # Subset data by random N number of segments (If Desired)

  ## get unique segments
  remaining_segments <- hs.input |> 
    dplyr::distinct(segment) |> # distinct() gets unique records from the desired field, segments
    dplyr::pull(segment)


  if (samples == "GOF") samples <- length(remaining_segments)
  if (!is.integer(samples)) as.integer(samples)

  if (replace) {
    sample_size <- samples
  } else {
    sample_size <- min(samples, length(remaining_segments))
  }

  ## Randomly select N unique segments
  subset_segments <- sample(remaining_segments,
    size = sample_size,
    replace = replace
  )

  res <- hs.input |>  # filter the original data for the segments selected in unique_segments
    dplyr::filter(segment %in% subset_segments) |>  # include only rows that match the selected segments
    dplyr::mutate(x = time, y = temperature) |> 
    sf::st_as_sf(coords = c("x", "y")) |>  # make a spatial feature where x and y are the spatial coordinates
    dplyr::group_by(segment) |> 
    dplyr::summarise(do_union = FALSE) |> 
    sf::st_cast("LINESTRING") |> 
    smoothr::densify(n = n, max_distance = max_distance) # this sets the number of points that will be added to each segment. This can be changed as desired


  res_coords <- sf::st_coordinates(res)
  lookup <- data.frame(L1 = unique(res_coords[, 3]), segment = dplyr::filter(hs.input, segment %in% subset_segments) |> dplyr::pull(segment) |> unique())

  res_coords |>
    dplyr::as_tibble() |>
    dplyr::left_join(lookup, dplyr::join_by(L1)) |>
    dplyr::rename(time = X, temperature = Y)  |> 
    dplyr::select(-dplyr::any_of(c("L1", "L2")))
}

#' Densify clustered paths
#'
#' @param x clustered t-T paths. Output of [cluster_paths()].
#' @inheritParams densify_paths
#'
#' @return tibble
#' @export
#'
#' @importFrom dplyr group_by left_join distinct select join_by
#'
#' @examples
#' \dontrun{
#' data(tT_paths1)
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.5)
#' cluster_paths(tT_paths_subset, cluster = 3) |>
#'   merge(tT_paths_subset, by = "segment") |>
#'   dplyr::group_by(cluster) |>
#'   densify_cluster()
#' }
densify_cluster <- function(x, GOF_rank = Inf, n = 10L, max_distance = 1, samples = 500, replace = TRUE) {
  time <- temperature <- segment <- cluster <- NULL
  x1 <- split(x, x$cluster, drop = TRUE)  |> 
    lapply(
      FUN = densify_paths,
      GOF_rank, n, max_distance, samples, replace
    )
  
    do.call(rbind, x1) |>
    dplyr::left_join(
      dplyr::select(x, -time, -temperature) |> dplyr::distinct(),
      dplyr::join_by(segment),
      relationship = "many-to-many"
    ) |>
    dplyr::group_by(cluster)
}
