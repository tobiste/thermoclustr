
<!-- README.md is generated from README.Rmd. Please edit that file -->

# thermochron

<!-- badges: start -->

[![R-CMD-check](https://github.com/tobiste/thermochron/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tobiste/thermochron/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/tobiste/thermochron/branch/main/graph/badge.svg)](https://app.codecov.io/gh/tobiste/thermochron?branch=main)
<!-- badges: end -->

The goal of thermochron is to provide tools to further analyze thermal
history models in thermochronology, including estimating cooling path
densities, and path families.

## Installation

You can install the development version of thermochron from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tobiste/thermochron")
```

## Example

The following minimal example shows how to import data, visualize the
paths and the path density as well as how to filter, and cluster the tT
paths into 3 path families.

``` r
library(thermochron)
library(dplyr)
library(ggplot2)

# load example dataset of a HeFTy model output
path2myfile <- system.file("112-72_30_H1_50-inv.txt", package = "thermochron")
tT_paths <- read_hefty(path2myfile)
```

The HeFTy model contains the modeled paths, the initial model
constraints, the weighted mean path, and some summary statistics on the
mineral grains. {thermochron} imports the model as a `list` that
contains the aforementioned features as list entries.

Thus, to visualize the paths with need to extract the `path` object from
the imported model `tT_paths`:

``` r
# (Optional) Create a subset of interest of the data:
tT_paths <- crop_paths(tT_paths, time = c(0, 400), temperature = c(0, 250)) 


tT_paths$paths |>
  ggplot(aes(time, temperature, color = Comp_GOF, group = segment)) +
  geom_path() +
  scale_color_viridis_c('GOF') +
  scale_x_reverse(position = "top") +
  scale_y_reverse()
```

<img src="man/figures/README-plot-1.png" width="100%" />

The path density can be visualized using `plot_path_density_filled()`:

``` r
plot_path_density_filled(tT_paths, show.legend = FALSE) +
  scale_x_reverse(position = "top") +
  scale_y_reverse()
```

<img src="man/figures/README-density-1.png" width="100%" />

To cluster the data, the following steps are required:

``` r
# Cluster the paths
paths_cluster <- cluster_paths(tT_paths, cluster = 2)

# Join with path dataset
paths_clustered <- merge(
  tT_paths$paths,
  paths_cluster,
  by = 'segment'
)
```

Finally, the visualization of the clustered tT paths:

``` r
paths_clustered |>
  ggplot(
    aes(
      x = time,
      y = temperature,
      color = cluster,
      alpha = Comp_GOF,
      group = segment
    )
  ) +
  geom_path() +
  scale_x_reverse(position = "top") +
  scale_y_reverse() +
  scale_color_viridis_d('Cluster')
```

<img src="man/figures/README-plot2-1.png" width="100%" />
