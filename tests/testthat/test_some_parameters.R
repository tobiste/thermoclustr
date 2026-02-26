data(tT_paths)


# test sampling
plot_path_density_filled(tT_paths, bins = 10, GOF_rank = 10, samples = 1, replace = FALSE, show.legend = FALSE)

plot_path_density_filled(tT_paths, bins = 10, GOF_rank = 10, samples = 1, replace = TRUE, show.legend = FALSE)

plot_path_density_filled(tT_paths, bins = 10, GOF_rank = 10, samples = 10, replace = TRUE, show.legend = FALSE)

plot_path_density_filled(tT_paths, bins = 10, GOF_rank = 10, samples = 100, replace = TRUE, show.legend = FALSE)

plot_path_density_filled(tT_paths, bins = 10, GOF_rank = 10, samples = 1000, show.legend = FALSE)

# test gof
plot_path_density_filled(tT_paths, bins = 10, GOF_rank = 1, samples = 100, show.legend = FALSE)

plot_path_density_filled(tT_paths, bins = 10, GOF_rank = 10, samples = 100, show.legend = FALSE)

plot_path_density_filled(tT_paths, bins = 10, GOF_rank = 1, samples = 100, show.legend = FALSE)

plot_path_density_filled(tT_paths, bins = 10, GOF_rank = 100, samples = 100, show.legend = FALSE)


tt_subset <- tT_paths$paths |>
  crop_paths(time = c(0, 400), temperature = c(0, 250)) |>
  dplyr::filter(Comp_GOF >= .1)

a <- path_diss(tt_subset, "Hausdorff")
b <- path_diss(tt_subset, "Frechet")

nbc <- path_nbclust(a)

cluster_paths(a, nbc$optimal, method = "hclust")
cluster_paths(a, nbc$optimal, method = "kmeans")
cluster_paths(a, nbc$optimal, method = "pam")
cluster_paths(a, nbc$optimal, method = "specc")
cluster_paths(a, nbc$optimal, method = "diana")
cluster_paths(a, nbc$optimal, method = "agnes")
cluster_paths(a, nbc$optimal, method = "clara", samples = 100)
cluster_paths(a, nbc$optimal, method = "fanny")
cluster_paths(a, method = "dbscan", minPts = 10, eps = 750)
cluster_paths(a, method = "hdbscan", minPts = 5)


path_hcut(a$diss, nbc$optimal, FUN = stats::hclust)
path_hcut(a$diss, nbc$optimal, FUN = cluster::agnes)
path_hcut(a$diss, nbc$optimal, FUN = cluster::diana)
