# Homogenize path family names from several samples

Compares cluster results from several samples because the
[`cluster_paths()`](https://tobiste.github.io/thermoclustr/reference/cluster_paths.md)
algorithm assigns arbitrary cluster names to the sample that can vary
for each clustered sample. The comparison is done by clustering the
median paths of the samples' cluster.

## Usage

``` r
cluster_samples(x, method = "hdbscan", minPts = 2, ...)
```

## Arguments

- x:

  named list containing data.frames with columns `segment`, `time`,
  `temperature` and `cluster`.

- method:

  the cluster algorithm applied. Default is `"hdbscan"`.

- minPts:

  minimum number of samples in a cluster (default is `2`). Only applied
  when `method = "hdbscan"` or `method = "dbscan"` is used.

- ...:

  optional arguments passed to
  [`cluster_paths()`](https://tobiste.github.io/thermoclustr/reference/cluster_paths.md).

## Value

list.

- 'paths':

  sf object of the median cooling paths of each sample's cluster

- 'diss:

  Dissimilarity matrix of the median paths

- 'hopkins':

  Hopkins statistic of the clustering

- 'dist':

  Used dissimilarity measure

- 'cluster':

  data.frame with the old and new cluster assignments of the median
  paths (segments).

- 'mds':

  data.frame with the multidimensional scaling coordinates of the
  dissimilarity of the median paths

- 'cluster_method':

  the clustering method used

## Details

The default cluster algorithm to assign comparable path families is
`"hdbscan"` trying to find clusters with at least 2 sample clusters.
This algorithm may also find "outliers" (`NA` values in
`cluster$cluster_new`), that are clusters not present in other samples.

## See also

[`path_statistics()`](https://tobiste.github.io/thermoclustr/reference/path_statistics.md)
for calculating the median paths,
[`path_diss()`](https://tobiste.github.io/thermoclustr/reference/path_diss.md)
for calculating the dissimilarity between the paths, and
[`cluster_paths()`](https://tobiste.github.io/thermoclustr/reference/cluster_paths.md)
for clustering cooling paths.

## Examples

``` r
# example data
data(tT_paths)
tT_paths_subset <- subset(tT_paths$paths, Comp_GOF >= 0.4)
# calculate dissimilarity and find 3 clusters:
tT_paths_subset_cl <- path_diss(tT_paths_subset) |> cluster_paths(3)
tT_paths_subset_all <- merge(tT_paths_subset, tT_paths_subset_cl, by = "segment")

# create two random subsets of the same sample
set.seed(20250411)

## select 100 random path segments
random_segments1 <- sample(unique(tT_paths$paths$segment), size = 100)
tT_paths_rnd1 <- subset(tT_paths$paths, segment %in% random_segments1)
tT_paths_rnd1_cl <- path_diss(tT_paths_rnd1) |>
  cluster_paths(3) # calculate dissimilarity and find 3 clusters
tT_paths_rnd1 <- merge(tT_paths_rnd1, tT_paths_rnd1_cl, by = "segment")

## select 100 random path segments
random_segments2 <- sample(unique(tT_paths$paths$segment), size = 100)
tT_paths_rnd2 <- subset(tT_paths$paths, segment %in% random_segments2)
tT_paths_rnd2_cl <- path_diss(tT_paths_rnd2) |>
  cluster_paths(4) # calculate dissimilarity and find 4 clusters
tT_paths_rnd2 <- merge(tT_paths_rnd2, tT_paths_rnd2_cl, by = "segment")

# combine all samples in a named list:
sample_list <- list("all" = tT_paths_subset_all, "s1" = tT_paths_rnd1, "s2" = tT_paths_rnd2)

# check how the individial clusters compare to each other:
cluster_samples(sample_list)
#> $paths
#> Simple feature collection with 10 features and 2 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: 0 ymin: 23.24454 xmax: 40 ymax: 65.43562
#> CRS:           NA
#> # A tibble: 10 × 3
#>                                                      geometry sample cluster_old
#>  *                                               <LINESTRING> <chr>  <chr>      
#>  1 (0.2154473 24.71499, 1.1591 26.55475, 1.944927 26.93293, … all    1          
#>  2 (0.004267409 23.96469, 1.063527 24.45831, 1.951099 24.770… all    2          
#>  3 (0.2304599 27.90557, 1.162849 33.02217, 2.026578 29.54804… all    3          
#>  4 (0.2107394 24.82631, 1.139026 25.46042, 1.95558 26.16914,… s1     1          
#>  5 (0 23.24454, 0.9786171 24.58917, 1.931717 23.94287, 2.881… s1     2          
#>  6 (0.1485138 27.78524, 1.18385 31.69537, 1.981107 31.4806, … s1     3          
#>  7 (0.2105611 24.82631, 1.208798 26.05979, 2.006259 25.92004… s2     1          
#>  8 (0.1854174 29.43847, 1.185072 36.67478, 1.999527 40.91297… s2     2          
#>  9 (0 23.38913, 0.9953817 24.08435, 1.97012 24.231, 2.771527… s2     3          
#> 10 (0 23.56792, 0.9552172 23.64332, 1.903507 23.90072, 2.861… s2     4          
#> 
#> $diss
#>            1         2         3         4         5         6         7
#> 2  11.857107                                                            
#> 3  12.193871 16.691665                                                  
#> 4  15.102393  4.468874 15.424051                                        
#> 5  32.854017 21.759057 29.055547 17.751625                              
#> 6  10.710525 15.714853  4.420437 14.945145 32.557487                    
#> 7  13.220186  6.416828 12.973788  3.648220 19.633832 12.923656          
#> 8  14.057358 19.955667  5.968636 19.346207 31.613409  6.125712 16.405143
#> 9  35.584225 24.489264 31.473238 20.481832  2.730208 35.287695 22.364039
#> 10 11.473066  7.799724 19.497439  9.574394 23.593637 19.362611 10.598462
#>            8         9
#> 2                     
#> 3                     
#> 4                     
#> 5                     
#> 6                     
#> 7                     
#> 8                     
#> 9  34.343616          
#> 10 22.056694 26.323751
#> 
#> $hopkins
#> statistic   p-value 
#> 0.5545917 0.6467965 
#> 
#> $dist
#> [1] "Hausdorff"
#> 
#> $mds
#>             MDS1       MDS2 sample cluster_new cluster_old
#> 1  -11.956042312  6.0997621    all           3           1
#> 2    0.165996586  6.3121044    all           2           2
#> 3   -9.105482180 -6.7360248    all           3           3
#> 4    2.744017829  2.4918257     s1           2           1
#> 5   19.651768903 -3.0584461     s1           1           2
#> 6  -12.756873749 -3.5504257     s1           3           3
#> 7    0.475673221  0.8690642     s2           2           1
#> 8  -11.572455541 -8.7502362     s2           3           2
#> 9   22.360112228 -3.7158714     s2           1           3
#> 10  -0.006714985 10.0382478     s2           2           4
#> 
#> $cluster
#>    sample cluster_new cluster_old
#> 1     all           3           1
#> 2     all           2           2
#> 3     all           3           3
#> 4      s1           2           1
#> 5      s1           1           2
#> 6      s1           3           3
#> 7      s2           2           1
#> 8      s2           3           2
#> 9      s2           1           3
#> 10     s2           2           4
#> 
#> $cluster_method
#> [1] "hdbscan"
#> 
#> attr(,"class")
#> [1] "list"   "tTdiss"
```
