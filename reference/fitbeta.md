# Fit of beta distribution to data

Estimates the shape parameters for a beta distribution using

## Usage

``` r
fitbeta(x, range = c(0, 1), na.rm = FALSE)
```

## Arguments

- x:

  numeric vector.

- range:

  numeric. two-element vector specifying the interval of the beta
  distribution

- na.rm:

  a logical value indicating whether `NA` values should be stripped
  before the computation proceeds.

## Value

mean, mode, and median (sorted by their values), and the variance of the
beta distribution fitted to `x`.

## References

https://en.wikipedia.org/wiki/Beta_distribution

## See also

[beta-stats](https://tobiste.github.io/thermoclustr/reference/beta-stats.md),
[stats::Beta](https://rdrr.io/r/stats/Beta.html)

## Examples

``` r
set.seed(20250411)
x <- stats::rbeta(1e5, 3.43, 5.55)
fitbeta(x)
#>   shape1   shape2 
#> 3.447885 5.586358 
```
