# Hausdorff Distance

Hausdorff distance (aka Hausdorff dimension)

## Usage

``` r
hausdorff(P, Q)

directed_hausdorff(P, Q)
```

## Arguments

- P, Q:

  numerical matrices, representing points in an m-dim. space.

## Value

A single scalar, the Hausdorff distance (dimension).

## Examples

``` r
P <- matrix(c(1, 1, 2, 2, 5, 4, 5, 4), 4, 2)
Q <- matrix(c(4, 4, 5, 5, 2, 1, 2, 1), 4, 2)
hausdorff(P, Q)
#> [1] 4.242641
```
