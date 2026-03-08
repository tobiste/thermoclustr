# Remove orphan paths

Sometimes HeFTy produces paths that do not start at the time = 0. These
"orphan" paths are not valid and should be removed. This function
removes paths that only exist at times greater than `min_time` (default
= 0).

## Usage

``` r
clean_HeFTy(x, min_time = 0)
```

## Arguments

- x:

  data.frame with time, temperature and segment columns

- min_time:

  numeric. Minimum time to consider a path valid (default is 0). Paths
  with segments that only exist at times greater than `min_time` will be
  removed.

## Value

data.frame with orphan paths removed.
