median_cluster_paths <- function(x, paths, nb) {
  cluster <- time_median <- temp_median <- NULL
  clustered <- cluster_paths(x, cluster = nb)

  clustered <- dplyr::left_join(paths, clustered, by = "segment")

  # If clustering assigned NA, drop them
  clustered <- clustered[!is.na(clustered$cluster), , drop = FALSE]

  clustered |>
    dplyr::group_by(cluster) |>
    densify_cluster() |>
    path_statistics() |>
    dplyr::mutate(
      .keep = "none",
      cluster,
      time = time_median,
      temperature = temp_median
    )
}

#' Homogenize path family names from several samples
#'
#' Compares cluster results from several samples because the [cluster_paths()]
#' algorithm assigns arbitrary cluster names to the sample that can vary for
#' each clustered sample.
#' The comparison is done by clustering the median paths of the samples' cluster.
#'
#' @param x named list containing data.frames with columns `segment`, `time`,
#' `temperature` and `cluster`.
#' @param method the cluster algorithm applied. Default is `"hdbscan"`.
#' @param minPts minimum number of samples in a cluster (default is `2`). Only
#' applied when `method = "hdbscan"` or `method = "dbscan"` is used.
#' @param ... optional arguments passed to [cluster_paths()].
#'
#' @details
#' The default cluster algorithm to assign comparable path families is `"hdbscan"`
#' trying to find clusters with at least 2 sample clusters. This algorithm may
#' also find "outliers" (`NA` values in `cluster$cluster_new`), that are clusters not
#' present in other samples.
#'
#'
#' @returns list.
#' \describe{
#' \item{'paths'}{sf object of the median cooling paths of each sample's cluster}
#' \item{'diss}{Dissimilarity matrix of the median paths}
#' \item{'hopkins'}{Hopkins statistic of the clustering}
#' \item{'dist'}{Used dissimilarity measure}
#' \item{'cluster'}{data.frame with the old and new cluster assignments of the
#' median paths (segments).}
#' \item{'mds'}{data.frame with the multidimensional scaling coordinates of the dissimilarity of the
#' median paths}
#' \item{'cluster_method'}{the clustering method used}
#' }
#' @export
#'
#' @importFrom data.table rbindlist
#' @importFrom dplyr mutate group_by
#'
#' @seealso [path_statistics()] for calculating the median paths,
#' [path_diss()] for calculating the dissimilarity between the paths, and
#' [cluster_paths()] for clustering cooling paths.
#'
#' @examples
#' # example data
#' data(tT_paths)
#' tT_paths_subset <- subset(tT_paths$paths, Comp_GOF >= 0.4)
#' # calculate dissimilarity and find 3 clusters:
#' tT_paths_subset_cl <- path_diss(tT_paths_subset) |> cluster_paths(3)
#' tT_paths_subset_all <- merge(tT_paths_subset, tT_paths_subset_cl, by = "segment")
#'
#' # create two random subsets of the same sample
#' set.seed(20250411)
#'
#' ## select 100 random path segments
#' random_segments1 <- sample(unique(tT_paths$paths$segment), size = 100)
#' tT_paths_rnd1 <- subset(tT_paths$paths, segment %in% random_segments1)
#' tT_paths_rnd1_cl <- path_diss(tT_paths_rnd1) |>
#'   cluster_paths(3) # calculate dissimilarity and find 3 clusters
#' tT_paths_rnd1 <- merge(tT_paths_rnd1, tT_paths_rnd1_cl, by = "segment")
#'
#' ## select 100 random path segments
#' random_segments2 <- sample(unique(tT_paths$paths$segment), size = 100)
#' tT_paths_rnd2 <- subset(tT_paths$paths, segment %in% random_segments2)
#' tT_paths_rnd2_cl <- path_diss(tT_paths_rnd2) |>
#'   cluster_paths(4) # calculate dissimilarity and find 4 clusters
#' tT_paths_rnd2 <- merge(tT_paths_rnd2, tT_paths_rnd2_cl, by = "segment")
#'
#' # combine all samples in a named list:
#' sample_list <- list("all" = tT_paths_subset_all, "s1" = tT_paths_rnd1, "s2" = tT_paths_rnd2)
#'
#' # check how the individial clusters compare to each other:
#' cluster_samples(sample_list)
cluster_samples <- function(x, method = "hdbscan", minPts = 2, ...) {
  temp_median <- time_median <- cluster <- segment <- geometry <- NULL
  stopifnot(
    is.list(x),
    all(sapply(x, inherits, "data.frame")),
    all(sapply(x, \(res) all(c("segment", "cluster") %in% names(res))))
  )

  median_paths <- Map(
    function(res, name) {
      # if(!('Comp_GOF' %in% names(res))) res$Comp_GOF <-  NA_real_
      res |>
        dplyr::group_by(cluster) |>
        densify_cluster() |>
        path_statistics() |>
        dplyr::mutate(.keep = "none", cluster, time = time_median, temperature = temp_median, sample = name)
    },
    x, names(x)
  )

  samples_diss <- data.table::rbindlist(median_paths) |>
    dplyr::mutate(segment = paste0(sample, "--", cluster)) |>
    path_diss()

  cluster_info <- cluster_paths(samples_diss, method = method, minPts = minPts, ...) |>
    dplyr::mutate(
      sample = sub("--.*", "", segment),
      # segment,
      cluster_new = factor(ifelse(cluster == "0", NA, as.character(cluster))),
      cluster_old = sub(".*--", "", segment),
      .keep = "none"
    )
  samples_diss$cluster <- cluster_info

  samples_diss$paths <- samples_diss$paths |>
    dplyr::mutate(
      sample = sub("--.*", "", segment),
      # segment,
      # cluster_new = factor(ifelse(cluster == "0", NA, as.character(cluster))),
      cluster_old = sub(".*--", "", segment),
      geometry,
      .keep = "none"
    )

  samples_mds <- stats::cmdscale(samples_diss$diss)
  samples_diss$mds <- data.frame(MDS1 = samples_mds[, 1], MDS2 = samples_mds[, 2], samples_diss$cluster)

  samples_diss$cluster_method <- method

  return(samples_diss)
}
