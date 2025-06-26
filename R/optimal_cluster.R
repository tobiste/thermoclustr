#' Locate elbow in elbow curve of wss (within-cluster sum of squares) 
#'
#' @param wss_values numeric vector of wss values.
#'
#' @returns integer index of the elbow point in the wss curve.
#' @noRd
find_elbow <- function(wss_values) {
  n_points <- length(wss_values)
  all_coords <- cbind(1:n_points, wss_values)
  
  # Line from first to last point
  line_vec <- all_coords[n_points, ] - all_coords[1, ]
  line_vec_norm <- line_vec / sqrt(sum(line_vec^2))
  
  # Compute the distance from each point to the line
  vec_from_first <- sweep(all_coords, 2, all_coords[1, ])
  scalar_product <- rowSums(vec_from_first * matrix(rep(line_vec_norm, n_points), ncol = 2, byrow = TRUE))
  proj <- outer(scalar_product, line_vec_norm)
  vec_to_line <- vec_from_first - proj
  dist_to_line <- sqrt(rowSums(vec_to_line^2))
  
  which.max(dist_to_line)
}


#' Optimal number of clusters
#' 
#' Convenience function to determine the optimal number of clusters in a dissimilarity matrix using three methods
#' wss, silhouette, and gap statistic.
#'
#' @param x either an object of class `"HeFTy"` (output of [read_hefty()]), 
#' `"tTdiss` (output of [path_diss()]), or 
#' a `data.frame` containing the `time`, `temperature`, and `segment` columns of the modeled paths.
#' @param FUNcluster cluster function
#' @param k.max integer. the maximum number of clusters to consider, must be at least two.
#' @param nboot integer. integer, number of Monte Carlo ("bootstrap") samples. 
#' Used only for determining the number of clusters using gap statistic.
#' @param dim integer. Number of dimensions for [stats::cmdscale()] used for the gap statistic, default is 4.
#' @param nstart integer. number of random starts for k-means clustering.
#' @param barfill,barcolor,linecolor fill color and outline color for bars
#' @param ... optionally further arguments for `FUNcluster()`
#'
#' @returns list containing the results of the three methods, the optimal number of clusters (median of the results), and the plots for each method.
#' @export
#' 
#' @importFrom factoextra fviz_nbclust
#' @importFrom ggplot2 labs geom_vline
#'
#' @examples
#' # example data
#' data(tT_paths1)
#' tT_paths_subset <- subset(tT_paths1$paths, Comp_GOF >= 0.4)
#' res <- path_nbclust(tT_paths_subset)
#' res$optimal
path_nbclust <- function(x, FUNcluster = factoextra::hcut, k.max = 5, nboot = 50, 
                         dim = 4, 
                                   barfill = "#1D1147FF",
                                   barcolor = "#1D1147FF",
                                   linecolor = "#1D1147FF", ...) {
  
  if (inherits(x, "HeFTy")) {
    x <- path_diss(x) # calculate dissimilarities
  } else if (inherits(x, "tTdiss")) {
  } else if(inherits(x, 'data.frame')){
    stopifnot(c('time', 'temperature', 'segment') %in% colnames(x))
    x <- path_diss(x)
  }
  diss <- x$diss
  diss_mat <- as.matrix(diss)
  plot_wss <- factoextra::fviz_nbclust(diss_mat, FUNcluster = FUNcluster,
                                          method = "wss", k.max = k.max, 
                                          verbose = TRUE, print.summary = TRUE,
                                          barfill = barfill,
                                          barcolor = barcolor,
                                          linecolor = linecolor,
                                       ...
  ) + ggplot2::labs(title = NULL)
  
  nb1 <- find_elbow(plot_wss$data$y)
  
  plot_wss <- plot_wss +
    ggplot2::geom_vline(xintercept = nb1, lty = 2, color = linecolor)
  
  
  #### Silhouette method
  plot_silh <- factoextra::fviz_nbclust(diss_mat, FUNcluster = FUNcluster,
                                          method = "silhouette", k.max = k.max, 
                                          verbose = FALSE, print.summary = FALSE,
                                          barfill = barfill,
                                          barcolor = barcolor,
                                          linecolor = linecolor,
                                        ...
  ) +
    ggplot2::labs(title = NULL)
  
  nb2 <- plot_silh$layers[[3]]$data |>
    unname() |>
    as.numeric()
  
  #### Gap statistic method
  plot_gapst <- factoextra::fviz_nbclust(stats::cmdscale(diss, k=dim), FUNcluster = FUNcluster,
                                          method = "gap_stat", 
                                         nboot = nboot,
                                         k.max = k.max, 
                                         verbose = FALSE, print.summary = FALSE,
                                         barfill = barfill,
                                         barcolor = barcolor,
                                         linecolor = linecolor,
                                         ...
  ) +
    ggplot2::labs(title = NULL)
  nb3 <- plot_gapst$layers[[4]]$data |>
    unname() |>
    as.numeric()
  
  results <- c("wss" = nb1, "silhouette" = nb2, "gap_stat" = nb3)
  
  opt_nb <- median(results)
  
  list(results = results, optimal = opt_nb, wss = plot_wss, silhouette = plot_silh, gap_stat = plot_gapst)
}
