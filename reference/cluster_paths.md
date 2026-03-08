# Clustering thermal histories models

Groups t-T paths into "path families" based on the *Hausdorff* or
*Fréchet distance* between paths.

## Usage

``` r
cluster_paths(
  x,
  k,
  dist = c("Hausdorff", "Frechet"),
  method = c("hclust", "kmeans", "pam", "dbscan", "hdbscan", "specc", "diana", "agnes",
    "clara", "fanny"),
  naming = c("asis", "GOF", "size"),
  warn = TRUE,
  threshold = 0.01,
  ...
)
```

## Arguments

- x:

  either an object of class `"HeFTy"` (output of
  [`read_hefty()`](https://tobiste.github.io/thermoclustr/reference/read_hefty.md)),
  an object of `"tTdiss"` (output of
  [`path_diss()`](https://tobiste.github.io/thermoclustr/reference/path_diss.md)),
  or a `data.frame` containing the `time`, `temperature`, and `segment`
  columns of the modeled paths.

- k:

  an integer scalar or vector with the desired number of groups. Ignored
  when `dist` equal to `dbscan` or `hdbscan` (**see Note**).

- dist:

  character. Algorithm to calculate a dissimilarity matrix (distance)
  for lines; either `Hausdorff` (the default), or `Frechet`.

- method:

  character. Clustering method to be applied. Currently implemented are

  `"hclust"`

  :   Hierarchical Clustering using
      [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html), the
      default)

  `"kmeans"`

  :   K-Means Clustering using
      [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html))

  `"pam"`

  :   Partitioning Around Medoids using
      [`cluster::pam()`](https://rdrr.io/pkg/cluster/man/pam.html)

  `"dbscan"`

  :   Density-based Spatial Clustering of Applications with Noise
      (DBSCAN) using
      [`dbscan::dbscan()`](https://rdrr.io/pkg/dbscan/man/dbscan.html)
      (**see Note**)

  `"hdbscan"`

  :   Hierarchical DBSCAN using
      [`dbscan::hdbscan()`](https://rdrr.io/pkg/dbscan/man/hdbscan.html)
      (**see Note**)

  `"specc"`

  :   Spectral Clustering using
      [`kernlab::specc()`](https://rdrr.io/pkg/kernlab/man/specc.html)

  `"agnes"`

  :   Agglomerative hierarchical clustering using
      [`cluster::agnes()`](https://rdrr.io/pkg/cluster/man/agnes.html)

  `"diana"`

  :   Divisive hierarchical clustering using
      [`cluster::diana()`](https://rdrr.io/pkg/cluster/man/diana.html)

  `"clara"`

  :   Clustering Large Applications using
      [`cluster::clara()`](https://rdrr.io/pkg/cluster/man/clara.html)

  `"fanny"`

  :   Fuzzy Analysis Clustering using
      [`cluster::fanny()`](https://rdrr.io/pkg/cluster/man/fanny.html)

- naming:

  character. Naming scheme for clusters. One of `"asis"` (output of
  underlying cluster algorithm), `"GOF"` (ranks of the mean GOF values
  within clusters), and `"size"` (ranks of the size of clusters).
  Outliers detected by `dbscan` of `hdbscan` will always named as `0`.

- warn:

  logical. Should there be a warning message if at least one cluster
  contains less than (`threshold` \* 100)% of the total paths?

- threshold:

  numeric. The significance threshold as a fraction of the total amount
  of paths (`0.01` by default). If one cluster contains less paths per
  total paths than this value, a warning message is created. Ignored if
  `warn` is `FALSE`.

- ...:

  additional arguments passed to cluster method.

## Value

a data.frame with the path `segment` (integer) and `cluster` (factor)

## Details

If you want to use a different clustering method that is not built in
the current version of `thermoclustr`, you can use the distance matrix
produced by
[`path_diss()`](https://tobiste.github.io/thermoclustr/reference/path_diss.md)
and feed your cluster algorithm.

## Note

that `dbscan` and `hdbscan` methods require `eps` and `minPts`
arguments. Optimal `eps` values can be visually estimated from the
"knee" in a k-nearest neighbor distance plot using
[`dbscan::kNNdistplot()`](https://rdrr.io/pkg/dbscan/man/kNNdist.html)'.

## See also

[`path_diss()`](https://tobiste.github.io/thermoclustr/reference/path_diss.md)
for calculating dissimilarities and
[`path_nbclust()`](https://tobiste.github.io/thermoclustr/reference/path_nbclust.md)
for determining the optimal number of clusters.

## Examples

``` r
data(tT_paths)
tT_paths$paths <- subset(tT_paths$paths, Comp_GOF >= 0.4)

# cluster the paths
cluster_paths(tT_paths, k = 3)
#>    segment cluster
#> 1        1       1
#> 2       10       2
#> 3      104       2
#> 4       11       3
#> 5       12       2
#> 6      122       3
#> 7       13       3
#> 8       14       2
#> 9       15       1
#> 10     152       1
#> 11      16       2
#> 12     164       3
#> 13     166       3
#> 14      17       2
#> 15     172       1
#> 16     177       1
#> 17      18       2
#> 18      19       3
#> 19     195       3
#> 20       2       2
#> 21      20       2
#> 22      21       1
#> 23      22       2
#> 24      23       1
#> 25      24       2
#> 26      25       3
#> 27     259       2
#> 28      26       3
#> 29     268       3
#> 30      27       3
#> 31     274       3
#> 32     276       2
#> 33      28       1
#> 34      29       1
#> 35     290       3
#> 36     291       3
#> 37     292       3
#> 38     293       2
#> 39       3       2
#> 40      30       2
#> 41      31       2
#> 42      32       2
#> 43      33       3
#> 44      34       3
#> 45      35       1
#> 46      36       1
#> 47      37       3
#> 48      38       1
#> 49      39       2
#> 50       4       3
#> 51      40       3
#> 52      42       3
#> 53      43       3
#> 54      44       1
#> 55      45       1
#> 56      46       3
#> 57      47       3
#> 58      48       1
#> 59      49       3
#> 60       5       2
#> 61      50       3
#> 62      51       3
#> 63      55       3
#> 64       6       2
#> 65      63       2
#> 66       7       2
#> 67       8       2
#> 68       9       1
```
