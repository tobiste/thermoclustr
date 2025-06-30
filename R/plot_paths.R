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
#' @param cluster (optional) data.frame containing the `segment` and `cluster` 
#' columns to be merged with `x`.
#' @param ... options passed to `pal`
#' @param pal color function
#' @param breaks breaks to cut the GOF values
#'
#' @returns plot
#' @importFrom viridisLite viridis
#' @importFrom graphics plot lines legend
#' @export
#'
#' @examples
#' # example data
#' data(tT_paths1)
#' tT_paths1$paths <- subset(tT_paths1$paths, Comp_GOF >= 0.4)
#' 
#' # Plot the paths
#' plot_paths(tT_paths1)
#'
#' # Show path clusters
#' plot_paths(tT_paths1, cluster = cluster_paths(tT_paths1, k = 3))
plot_paths <- function(x, cluster = NULL, pal = viridisLite::viridis, breaks = 5, ...){
  if(inherits(x, 'HeFTy')) x <- x$paths
  stopifnot(c("segment",  "time", "temperature") %in% colnames(x))
  
  if(!is.null(cluster)){
  x <- merge(x, cluster, by = "segment" )
  cl_lab <- sort(unique(x$cluster))
  cols <- setNames(pal(length(cl_lab), ...), cl_lab)
  x$col <- cols[as.character(x$cluster)]
  #
  } else {
    if("Comp_GOF" %in% colnames(x)) {
      breaks <- pretty(x$Comp_GOF, n = breaks)
      n2 <- length(breaks) - 1
      col_labels <- cut(x$Comp_GOF, breaks = breaks, include.lowest = TRUE)
      cols <- pal(n2, ...)
      col_map <- setNames(cols, levels(col_labels))
      x$col <- col_map[as.character(col_labels)]
      
      col_legend <- data.frame(
        val = levels(col_labels),
        col = cols,
        stringsAsFactors = FALSE
      )
      #
    } else {
    x$col <- 'grey20'
    }
  }
  graphics::plot(
    x$time, 
    x$temperature, 
    type = "n", 
       xlim = rev(range(x$time)), 
       ylim = rev(range(x$temperature)), 
       # xlab = "X", ylab = "Y", main = "Multiple Lines Plot"
       )
  
  by(x, x$segment, function(seg) {
    graphics::lines(seg$time, seg$temperature, col = seg$col[1])
  })
  
  
  if(!is.null(cluster)){
  graphics::legend(
    "bottomright", 
    title = "Cluster", 
    legend = cl_lab, 
    col = cols, 
    pch = 19,
    bty = "o", bg = "white"
    )
    #
  } else {
    if("Comp_GOF" %in% colnames(x)){
      graphics::legend(
        "bottomright",
        title = "Goodness-of-Fit",
        # inset = .05, cex = .75,
        legend = col_legend$val,  
        fill = col_legend$col,
        bty = "o", bg = "white"
      )
    }
  }
}