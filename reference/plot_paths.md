# Plot the thermochronology paths

Plots the time-temperature paths. If the input `x` is a `"HeFTy"` object
or a data.frame containing the goodness of fit values (as. `Comp_GOF`
column), the paths are color coded by the goodness of fit. If `cluster`
is specified, the paths are color coded by the clusters.

## Usage

``` r
plot_paths(
  x,
  cluster = NULL,
  do.cluster = FALSE,
  cluster.params = list(),
  pal = viridisLite::viridis,
  breaks = 5,
  ...
)
```

## Arguments

- x:

  either an object of class `"HeFTy"` (output of
  [`read_hefty()`](https://tobiste.github.io/thermoclustr/reference/read_hefty.md))
  or a `data.frame` containing the `time`, `temperature`, and `segment`
  columns of the modeled paths.

- cluster:

  (optional) Either a data.frame containing the `segment` and `cluster`
  columns to be merged with data in `x`

- do.cluster:

  logical. Whether clustering with arguments specified in
  `cluster.params` should be applied. This will add information about
  the used dissimilarity measure and Hopkins statistic to the plot. Is
  ignored if `cluster` is specified.

- cluster.params:

  list. Arguments passed to
  [`cluster_paths()`](https://tobiste.github.io/thermoclustr/reference/cluster_paths.md)
  function. Only effective when `do.cluster = TRUE` and
  `cluster = NULL`.

- pal:

  color function

- breaks:

  breaks to cut the GOF values

- ...:

  options passed to `pal`

## Value

plot

## Examples

``` r
# example data
data(tT_paths)

# Plot the paths:
plot_paths(tT_paths)


# Show predefined path clusters:
cl <- cluster_paths(tT_paths, 2)
plot_paths(tT_paths, cluster = cl)


# Calculate cluster while plotting:
plot_paths(tT_paths, do.cluster = TRUE, cluster.params = list(k = 3, method = "pam"))
```
