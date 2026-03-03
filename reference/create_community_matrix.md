# create community matrix

Using your community_object, we are able to convert it into a community
matrix for easier usability of the object.

## Usage

``` r
create_community_matrix(cluster_object)
```

## Arguments

- cluster_object:

  the result of the
  [`cluster_data()`](https://www.mums2.org/mums2/reference/cluster_data.md)
  function. data \<- import_all_data(peak_table =
  mums2::mums2_example("botryllus_pt_small.csv"), meta_data =
  mums2::mums2_example("meta_data_boryillus.csv"), format = "None")

  filtered_data \<- data \|\>
  filter_peak_table(filter_mispicked_ions_params()) \|\>
  filter_peak_table(filter_cv_params(cv_threshold = 0.2)) \|\>
  filter_peak_table(filter_group_params(group_threshold = 0.1,
  "Blanks")) \|\> filter_peak_table(filter_insource_ions_params())
  change_rt_to_seconds_or_minute(filtered_data, "minutes")

  matched_data \<-
  ms2_ms1_compare(mums2_example("botryllus_v2.gnps.mgf"), filtered_data,
  10, 6)

  dist \<- dist_ms2(data = matched_data, cutoff = 0.3, precursor_thresh
  = 2, score_params = modified_cosine_params(0.5), min_peaks = 0)

  cluster_results \<- cluster_data(distance_df = dist, ms2_match_data =
  matched_data, cutoff = 0.3, cluster_method = "opticlust")

  community_matrix \<- create_community_matrix_object(cluster_results)

## Value

a `data.frame` object of your community_object.
