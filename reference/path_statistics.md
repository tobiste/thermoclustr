# Path statistics

Provides binned statistics (mean, median, IQR, quantiles, etc.) of
modeled t-T paths.

## Usage

``` r
path_statistics(x, w = NULL, breaks = 50)
```

## Arguments

- x:

  either an object of class `"HeFTy"` (output of
  [`read_hefty()`](https://tobiste.github.io/thermoclustr/reference/read_hefty.md))
  or a `data.frame` containing the `time` and `temperature` columns of
  the modeled paths.

- w:

  numeric vector. Weights for each path segment.

- breaks:

  either a numeric vector of two or more unique cut points or a single
  number (greater than or equal to 2) giving the number of intervals
  into which `x` is to be cut.

## Value

tibble.

## See also

[`gof_weighting()`](https://tobiste.github.io/thermoclustr/reference/gof_weighting.md)
for rescaling goodness-of-fit values to weights.

## Examples

``` r
data(tT_paths)
path_statistics(tT_paths, w = gof_weighting(tT_paths$paths$Comp_GOF))
#> # A tibble: 50 × 12
#>    bins     time_min time_median time_max temp_mean temp_sd temp_IQR temp_median
#>    <fct>       <dbl>       <dbl>    <dbl>     <dbl>   <dbl>    <dbl>       <dbl>
#>  1 (-0.04,…    0            0       0.786      25.1    5.34     1.65        23.9
#>  2 (0.8,1.…    0.833        1.25    1.59       34.9   10.8     16.5         32.9
#>  3 (1.6,2.…    1.62         1.98    2.36       34.9    9.92    16.5         34.2
#>  4 (2.4,3.…    2.43         2.71    3.14       30.9    8.72    10.3         28.8
#>  5 (3.2,4]     3.20         3.59    3.95       33.0   10.7      8.32        28.6
#>  6 (4,4.8]     4.06         4.50    4.77       35.5   11.5     16.0         34.0
#>  7 (4.8,5.…    4.82         5.11    5.31       37.5   10.1     12.6         39.1
#>  8 (5.6,6.…    5.62         5.99    6.38       36.4   11.3     21.5         36.2
#>  9 (6.4,7.…    6.43         6.86    7.15       42.1   12.8     20.2         40.0
#> 10 (7.2,8]     7.23         7.53    7.87       43.5   10.3     23.4         46.6
#> # ℹ 40 more rows
#> # ℹ 4 more variables: temp_max <dbl>, temp_min <dbl>, temp_5 <dbl>,
#> #   temp_95 <dbl>
```
