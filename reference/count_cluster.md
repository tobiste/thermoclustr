# Count the number of paths in each cluster

Helper function to calculate the paths in each cluster and give a
warning message if one cluster contains not enough paths to be
considered significant

## Usage

``` r
count_cluster(x)
```

## Arguments

- x:

  data.frame with a `$cluster` column

## Value

named array of integers, the number of paths in each cluster

## Examples

``` r
data(tT_paths)
tT_paths$paths <- subset(tT_paths$paths, Comp_GOF >= 0.4)

# cluster the paths
res <- cluster_paths(tT_paths, k = 3)
count_cluster(res)
#>  1  2  3 
#> 16 24 28 
```
