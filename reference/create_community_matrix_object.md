# Create Community Matrix Object.

Using the data generated from clustering or adding ms2 data to your
object, we are able to create a community matrix object. The community
matrix object stores the same data a community matrix but within a cpp
object. We use this object to conduct analysis more efficiently.

## Usage

``` r
create_community_matrix_object(data)

# S3 method for class 'mass_data'
create_community_matrix_object(data)

# S3 method for class 'mothur_cluster'
create_community_matrix_object(data)
```

## Arguments

- data:

  the result of the
  [`cluster_data()`](https://www.mums2.org/mums2/reference/cluster_data.md)
  function, or just a mass_data object created from
  [`ms2_ms1_compare()`](https://www.mums2.org/mums2/reference/ms2_ms1_compare.md).

## Value

a external pointer to an Rcpp object.

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

dist <- dist_ms2(data = matched_data, cutoff = 0.3, precursor_thresh = 2,
 score_params = modified_cosine_params(0.5), min_peaks = 0,
 number_of_threads = 2)

cluster_results <- cluster_data(distance_df = dist,
 ms2_match_data = matched_data, cutoff = 0.3, cluster_method = "opticlust")
#> Warning: [WARNING]: The mcc metric is not suitible for your data with a cutoff of 0.300000 using tptn instead.

community_with_cluster <- create_community_matrix_object(cluster_results)

community_object_mass_data <- create_community_matrix_object(matched_data)
```
