# Cluster tendency of thermochronology cooling paths

Calculate the Hopkins statistic of Hausdorff or Fréchet distance
matrices to check clusterability of thermochronologic cooling paths.
Calculated values 0-0.3 indicate regularly-spaced data. Values around
0.5 indicate random data. Values 0.7-1 indicate clustered data.

## Usage

``` r
cluster_tendency(x, m = NULL, method = c("simple", "torus"))
```

## Arguments

- x:

  either an object of class `"tTdiss"` (output of
  [`path_diss()`](https://tobiste.github.io/thermoclustr/reference/path_diss.md)),
  a distance matrix (object `"dist"`), or a `data.frame` containing
  coordinates from a dimension-reducing algorithm.

- m:

  integer. Number of rows to sample from `x`. If `NULL`, all rows are
  considered.

- method:

  character. Either `"simple"` (default) or `"torus"`.

## Value

numeric. The value of Hopkins statistic in the range of 0 and 1 and its
p-value (under the null hypothesis of spatial randomness).

## Details

Calculates the Hopkins statistic on the transformed Hausdorff or Fréchet
distance matrix using metric multidimensional scaling by default.

## Examples

``` r
data(tT_paths)
tT_paths_subset <- subset(tT_paths$paths, Comp_GOF >= 0.4)

# calculate the dissimilarities of the paths:
tT_diss <- path_diss(tT_paths_subset, densify = 1)

set.seed(20250411)
cluster_tendency(tT_diss) # H=0.77, p=2.88e-11
#>    statistic      p-value 
#> 7.655771e-01 2.883072e-11 

# Alternative dimension reducing algorithms can be used:
# Kriskal's non-metric MDS
library(MASS)
tT_nmds <- isoMDS(tT_diss$diss)$points
#> initial  value 11.693755 
#> iter   5 value 9.168128
#> iter  10 value 8.906481
#> final  value 8.850736 
#> converged
cluster_tendency(tT_nmds)
#>    statistic      p-value 
#> 0.6601201750 0.0001442679 
# Uniform Manifold Approximation and Projection (UMAP)  
library(uwot)
tT_umap <- umap2(tT_diss$diss)
cluster_tendency(tT_umap)
#> statistic   p-value 
#>  0.908872  0.000000 
```
