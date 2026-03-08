# Densify Paths

Adds extra points between the vertices of the path segments

## Usage

``` r
densify_paths(
  x,
  GOF_rank = 10L,
  n = 10L,
  max_distance = 1,
  samples = 100L,
  replace = FALSE
)
```

## Arguments

- x:

  either an object of class `"HeFTy"` (output of
  [`read_hefty()`](https://tobiste.github.io/thermoclustr/reference/read_hefty.md)),
  or a `data.frame` containing the `time`, `temperature`, and `Comp_GOF`
  columns of the modeled paths.

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

## Note

A large sample number `n` will require a **long(!)** processing time for
this function and subsequent methods such as
[`plot_path_density()`](https://tobiste.github.io/thermoclustr/reference/plt_density.md)
or
[`cluster_paths()`](https://tobiste.github.io/thermoclustr/reference/cluster_paths.md).

If only paths within a specified GOF range should be densified, create a
subset of the data beforehand using either
[`subset()`](https://rdrr.io/r/base/subset.html) or
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html).

## Examples

``` r
data(tT_paths)
densify_paths(tT_paths)
#> # A tibble: 593 × 3
#>     time temperature segment
#>    <dbl>       <dbl> <chr>  
#>  1  40          66.5 1      
#>  2  39.3        66.3 1      
#>  3  38.6        66.1 1      
#>  4  37.9        66.0 1      
#>  5  37.9        65.0 1      
#>  6  37.8        64.1 1      
#>  7  37.8        63.1 1      
#>  8  37.8        62.2 1      
#>  9  37.8        61.2 1      
#> 10  37.8        60.3 1      
#> # ℹ 583 more rows
```
