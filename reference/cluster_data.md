# Cluster Features

`cluster_data()` allows users to cluster features inside the mass data
object. This is done by creating a sparse matrix using the `distMs2()`
function and inputting that inside the clutur package. This allows us to
easily cluster features that contain an ms2 spectra.

## Usage

``` r
cluster_data(
  distance_df,
  ms2_match_data,
  cutoff = 0.3,
  cluster_method = "opticlust"
)
```

## Arguments

- distance_df:

  a distance df that was generated from the `distMs2()` function.

- ms2_match_data:

  your mass data set object generated from
  [`ms2_ms1_compare()`](https://www.mums2.org/mums2/reference/ms2_ms1_compare.md).

- cutoff:

  the cutoff value you wish to cluster to.

- cluster_method:

  the clustering algorithm you wish to use. The options are: furthest,
  nearest, weighted, average, and opticlust.

## Value

a shared `data.frame` (or a `mothur_cluster` object) displaying all the
clustered and abundance data.

## Examples

``` r
data <-
   import_all_data(peak_table =
                   mums2::mums2_example("botryllus_pt_small.csv"),
                   metadata =
                   mums2::mums2_example("boryillus_metadata.csv"),
                   format = "None")


matched_data <- ms2_ms1_compare(mums2_example("botryllus_v2.gnps.mgf"),
 data, 1, 6)
#> Reading: /home/runner/work/_temp/Library/mums2/extdata/botryllus_v2.gnps.mgf ...
#> 17/349 peaks have an MS2 spectra.

dist <- dist_ms2(data = matched_data, cutoff = 0.6, precursor_thresh = 100,
 score_params = modified_cosine_params(0.5), min_peaks = 0,
 number_of_threads = 2)

cluster_results <- cluster_data(distance_df = dist,
 ms2_match_data = matched_data, cutoff = 0.3, cluster_method = "opticlust")
```
