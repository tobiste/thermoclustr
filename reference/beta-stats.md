# Summary statistics of a beta distribution

Summary statistics of a beta distribution

## Usage

``` r
beta_mean(x, shape1 = NULL, shape2 = NULL)

beta_var(x, shape1 = NULL, shape2 = NULL)

beta_mode(x, shape1 = NULL, shape2 = NULL)

beta_summary(x, ...)
```

## Arguments

- x:

  numeric vector

- shape1, shape2:

  non-negative parameters of the Beta distribution.

- ...:

  optional arguments passed to
  [`fitbeta()`](https://tobiste.github.io/thermoclustr/reference/fitbeta.md)

## Value

numeric.

## References

https://en.wikipedia.org/wiki/Beta_distribution

## See also

[`fitbeta()`](https://tobiste.github.io/thermoclustr/reference/fitbeta.md),
[stats::Beta](https://rdrr.io/r/stats/Beta.html)

## Examples

``` r
set.seed(20250411)
x <- stats::rbeta(1e5, 3.43, 5.55)
beta_summary(x)
#>       Mode     Median       Mean   Variance 
#> 0.34799555 0.37221674 0.38164629 0.02351871 
```
