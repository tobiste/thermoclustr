# Hierarchical Clustering of Cooling Paths

Convenience function to compute hierarchical clustering and cut the tree
into k clusters

## Usage

``` r
path_hcut(x, k, FUN = stats::hclust, ...)
```

## Arguments

- x:

  an object of class `"HeFTy"`, `"tTdiss`" or `"dist"` (dissimilarity
  matrix).

- k:

  integer. number of clusters to be generated

- FUN:

  hierarchical clustering function to be used, i.e. one of
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) (the
  default),
  [`cluster::agnes()`](https://rdrr.io/pkg/cluster/man/agnes.html),
  [`cluster::diana()`](https://rdrr.io/pkg/cluster/man/diana.html)).

- ...:

  optional arguments past to `hc_func`

## Value

an object of class `"hcut"` containing the result of the standard
function used (read the documentation of
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html),
[`cluster::agnes()`](https://rdrr.io/pkg/cluster/man/agnes.html),
[`cluster::diana()`](https://rdrr.io/pkg/cluster/man/diana.html)).

It includes also

- cluster:

  the cluster assignment of observations after cutting the tree

- nbclust:

  the number of clusters

- size:

  the size of clusters

## Examples

``` r
data(tT_paths)
tT_paths_subset <- subset(tT_paths$paths, Comp_GOF >= 0.4)

# calculate the dissimilarities of the paths:
tT_diss <- path_diss(tT_paths_subset)
path_hcut(tT_diss$diss, 3)
#> 
#> Call:
#> FUN(d = x)
#> 
#> Cluster method   : complete 
#> Number of objects: 68 
#> 
```
