#' Plot the thermochronology paths
#'
#' Plots the time-temperature paths. If the input `x` is a `"HeFTy"` object or a
#'  data.frame containing the goodness of fit values (as. `Comp_GOF` column),
#'  the paths are color coded by the goodness of fit. If `cluster` is specified,
#'  the paths are color coded by the clusters.
#'
#' @param x either an object of class `"HeFTy"` (output of [read_hefty()]) or
#' a `data.frame` containing the `time`, `temperature`, and `segment` columns
#' of the modeled paths.
#' @param cluster (optional) Either a data.frame containing the `segment` and `cluster`
#' columns to be merged with data in `x`
#' @param do.cluster logical. Whether clustering with arguments specified in
#' `cluster.params` should be applied. This will add information about the used
#' dissimilarity measure and Hopkins statistic to the plot. Is ignored if
#' `cluster` is specified.
#' @param ... options passed to `pal`
#' @param pal color function
#' @param breaks breaks to cut the GOF values
#' @param cluster.params list. Arguments passed to [cluster_paths()] function.
#' Only effective when `do.cluster = TRUE` and `cluster = NULL`.
#'
#' @returns plot
#' @importFrom viridisLite viridis
#' @importFrom graphics plot lines legend title mtext
#' @export
#'
#' @examples
#' # example data
#' data(tT_paths)
#'
#' # Plot the paths:
#' plot_paths(tT_paths)
#'
#' # Show predefined path clusters:
#' cl <- cluster_paths(tT_paths, 2)
#' plot_paths(tT_paths, cluster = cl)
#'
#' # Calculate cluster while plotting:
#' plot_paths(tT_paths, do.cluster = TRUE, cluster.params = list(k = 3, method = "pam"))
plot_paths <- function(x, cluster = NULL, do.cluster = FALSE, cluster.params = list(), pal = viridisLite::viridis, breaks = 5, ...) {
  if (inherits(x, "HeFTy")) x <- x$paths
  stopifnot(c("segment", "time", "temperature") %in% colnames(x))
  has_info <- FALSE
  k <- method <- Comp_GOF <- time <- temperature <- NULL

  if (is.null(cluster)) {
    if (do.cluster) {
      has_info <- TRUE
      cluster.params1 <- append(cluster.params, list(x = x), 0)
      if (!is.null(cluster.params1$k)) cluster.params1$k <- NULL
      if (!is.null(cluster.params1$method)) cluster.params1$method <- NULL
      diss <- do.call(path_diss, args = cluster.params1)
      cluster <- do.call(cluster_paths, args = append(cluster.params, list(x = diss), 0))

      info_dist <- diss$dist
      info_hopkins <- round(diss$hopkins, 2)
      info_clustmeth <- ifelse(is.null(cluster.params$method), "hclust", cluster.params$method)
      do.cluster <- FALSE
    } else {
      if ("Comp_GOF" %in% colnames(x)) {
        breaks <- pretty(x$Comp_GOF, n = breaks)
        n2 <- length(breaks) - 1
        col_labels <- cut(x$Comp_GOF, breaks = breaks, include.lowest = TRUE)
        cols <- pal(n2, ...)
        col_map <- stats::setNames(cols, levels(col_labels))
        x$col <- col_map[as.character(col_labels)]
        #
      } else {
        x$col <- "grey20"
      }
    }
  }

  if (!is.null(cluster)) {
    x <- merge(x, cluster, by = "segment")
    n_cluster <- by(x, x$cluster, nrow) |> rbind()
    cl_lab <- sort(unique(x$cluster))
    cols <- stats::setNames(pal(length(cl_lab), ...), cl_lab)

    x$col <- cols[as.character(x$cluster)]
    cl_lab2 <- paste0(cl_lab, " (", n_cluster[cl_lab], ")")
  }

  graphics::plot(
    x$time,
    x$temperature,
    type = "n",
    xlim = rev(range(x$time)),
    ylim = rev(range(x$temperature)),
    xlab = "Time (Ma)", ylab = expression("Temperature (" * degree * "C)"), # main = "Multiple Lines Plot"
  )

  by(x, x$segment, function(seg) {
    graphics::lines(seg$time, seg$temperature, col = seg$col[1])
  })


  if (!is.null(cluster)) {
    graphics::legend(
      "bottomright",
      title = ifelse(has_info, info_clustmeth, "Cluster"),
      legend = cl_lab2,
      col = cols,
      pch = 19,
      bty = "o", bg = "white"
    )
    graphics::title(main = "Path families")
    if (has_info) {
      graphics::title(sub = paste0("Hopkins statistic: ", info_hopkins[1], " (p=", info_hopkins[2], ")"))
      graphics::mtext(paste0(info_dist, " distance"))
    }

    #
  } else {
    if ("Comp_GOF" %in% colnames(x)) {
      graphics::legend(
        "bottomright",
        title = "Goodness-of-Fit",
        # inset = .05, cex = .75,
        legend = levels(col_labels),
        fill = cols,
        bty = "o", bg = "white"
      )
      graphics::title(main = "Inverse model thermal history paths")
    }
  }
}
