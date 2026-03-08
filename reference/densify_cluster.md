# Densify clustered paths

Densify clustered paths

## Usage

``` r
densify_cluster(
  x,
  GOF_rank = Inf,
  n = 10L,
  max_distance = 1,
  samples = 500,
  replace = TRUE
)
```

## Arguments

- x:

  clustered t-T paths. Output of
  [`cluster_paths()`](https://tobiste.github.io/thermoclustr/reference/cluster_paths.md)
  merged (or joined) with paths.

- GOF_rank:

  numeric. Selects only the `GOF_rank`-th highest GOF ranked paths. If
  all GOFs should be used, set to `Inf`. Default is 10.

- n:

  integer. Adds `n` (10 by default) equally-spaced extra points along
  each path segment (between vertices).

- max_distance:

  numeric. Adds points at a maximum distance of `max_distance` (in Myr)
  from each other. 1 by default.

- samples:

  integer or character. Number of random samples of the data. This
  number should be less or equal then the amount of paths. The default
  is `100`. Paths will be randomly selected **after** the data has been
  filtered by the `GOF_rank()`. Optional, set `samples` to `'all'` to
  consider all paths ignoring the `GOF_rank()` filter (this sets `GOF`
  to `Inf`). Set `samples` to `GOF` to consider all paths after the
  `GOF_rank` filter.

- replace:

  logical. Should sampling be with replacement?

## Value

tibble

## Examples

``` r
if (FALSE) { # \dontrun{
data(tT_paths)
tT_paths_subset <- subset(tT_paths$paths, Comp_GOF >= 0.5)
cluster_paths(tT_paths_subset, cluster = 3) |>
  merge(tT_paths_subset, by = "segment") |>
  dplyr::group_by(cluster) |>
  densify_cluster()
} # }
```
