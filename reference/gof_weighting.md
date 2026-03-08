# Create weightings from goodness-of-fit values

HeFTy goodness-of-fit (GOF) values range from 0 (worst fit) to 1 (best
fit). Because, GOF values greater than 0.5 are considered equally good
("acceptable"), this function rescales the GOF values to weighting
values.

## Usage

``` r
gof_weighting(x, na.rm = FALSE)
```

## Arguments

- x:

  numeric vector of GOF values (between 0 and 1)

- na.rm:

  logical. Whether missing values should be removed (`FALSE` by
  default.)

## Value

numeric vector of weights between 0 and 1

## See also

[`plot_path_density_filled_weighted()`](https://tobiste.github.io/thermoclustr/reference/plt_density.md)
for calculating weighted density of t-T paths.

## Examples

``` r
gof_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)
gof_weighting(gof_values)
#> [1] 0.2 0.6 1.0 1.0 1.0
```
