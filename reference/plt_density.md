# Path Density Plot

Creates a 2d kernel density estimate of the t-T paths and plots it using
ggplot

## Usage

``` r
plot_path_density_filled(
  x,
  bins = 50L,
  densify = TRUE,
  show.legend = NA,
  geom = "density_2d_filled",
  n = 100L,
  ...
)

plot_path_density(
  x,
  bins = 50L,
  densify = TRUE,
  show.legend = NA,
  n = 100L,
  ...
)

plot_path_density_filled_weighted(
  x,
  bins = 50L,
  densify = TRUE,
  show.legend = NA,
  n = 100L,
  weights = NULL,
  ...
)
```

## Arguments

- x:

  either an object of class `"HeFTy"` (output of
  [`read_hefty()`](https://tobiste.github.io/thermoclustr/reference/read_hefty.md)),
  a `data.frame` containing the `time`, `temperature` columns of the
  modeled paths, or output of
  [`densify_paths()`](https://tobiste.github.io/thermoclustr/reference/densify_paths.md).

- bins:

  integer. Amount of bins used for the kernel density estimate. 50 by
  default.

- densify:

  logical. Whether the paths in `x` should be densified first? Default
  is `TRUE`.

- show.legend:

  logical. Should this layer be included in the legends? `NA`, the
  default, includes if any aesthetics are mapped. `FALSE` never
  includes, and `TRUE` always includes. It can also be a named logical
  vector to finely select the aesthetics to display.

- geom:

  Use to override the default connection between
  [`ggplot2::geom_density_2d()`](https://ggplot2.tidyverse.org/reference/geom_density_2d.html)
  and
  [`ggplot2::stat_density_2d()`](https://ggplot2.tidyverse.org/reference/geom_density_2d.html).
  For more information at overriding these connections, see how the stat
  and geom arguments work.

- n:

  integer. Number of grid points in each direction.

- ...:

  Arguments passed on to
  [`densify_paths()`](https://tobiste.github.io/thermoclustr/reference/densify_paths.md)
  (only if `densify=TRUE`).

- weights:

  numeric vector. Weights for each path segment.

## Value

ggplot

## See also

[`gof_weighting()`](https://tobiste.github.io/thermoclustr/reference/gof_weighting.md)
for rescaling goodness-of-fit values to weights.

## Examples

``` r
data(tT_paths)
plot_path_density(tT_paths)

plot_path_density_filled(tT_paths)

plot_path_density_filled(tT_paths, geom = "raster")

plot_path_density_filled_weighted(tT_paths, weights = gof_weighting(tT_paths$paths$Comp_GOF))
```
