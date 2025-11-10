# Calculate pairwise distance between MS/MS features.

`dist_ms2` calculates and stores all non-zero distance values above the
user defined cutoff (default = 0.3).

## Usage

``` r
dist_ms2(
  data,
  cutoff,
  precursor_threshold,
  score_params,
  min_peaks = 6,
  number_of_threads = detectCores()
)
```

## Arguments

- data:

  the object generated from
  [`ms2_ms1_compare()`](https://www.mums2.org/mums2/reference/ms2_ms1_compare.md).

- cutoff:

  The maximum distance value (`numeric`) to store a pairwise comparison.
  The default of .3 corresponds to a cosine score of .7, meaning pairs
  with a score of .7 or higher will be stored in the matrix.

- precursor_threshold:

  Precursor mz tolerance. MS2 scans with a difference in precursor mz
  less than or equal to this value will be scored. Disable this by
  setting this value to -1 or less.

- score_params:

  Parameters for scoring method to be applied. See
  [`modified_cosine_params()`](https://www.mums2.org/mums2/reference/modified_cosine_params.md)
  and
  [`spec_entropy_params()`](https://www.mums2.org/mums2/reference/spec_entropy_params.md)
  for more details.

- min_peaks:

  the minimum number of peaks that need to be present before you compare
  the ms2 spectra.

- number_of_threads:

  the number of threads you want to use for this calculation.

## Value

A sparse matrix of class `"data.frame"`

## Details

This function takes a `mass_data` object as input and calculates
distance between ms2 peaks. Currently, MS1 features without MS2 peaks
returns no distance value. Distance can be calculated with method
`"gnps"` or `"spectral_entropy"`. A sparse matrix is returned.

## Examples

``` r
data <-
   import_all_data(peak_table =
                   mums2::mums2_example("full_mix_peak_table_small.csv"),
                   meta_data =
                   mums2::mums2_example("full_mix_meta_data_small.csv"),
                   format = "Metaboscape")
#> If peak table has corrupted compound names they will be converted to
#>       utf-8 and if there are any commas, they will be converted to periods(.).

filtered_data <- data |>
   filter_peak_table(filter_mispicked_ions_params()) |>
   filter_peak_table(filter_cv_params(cv_threshold = 0.2)) |>
   filter_peak_table(filter_group_params(group_threshold = 0.1,
                                             "Blanks")) |>
   filter_peak_table(filter_insource_ions_params())
#> ℹ Checking 1294 peaks for mispicked peaks.
#> ℹ Argument merge_peaks is: TRUE. Merging mispicked peaks with method sum.
#> ✔ 8 ions failed the mispicked filter, 1286 ions remain.
#> ℹ Parsing 1286 peaks for replicability across technical replicates.
#> ✔ 399 ions failed the cv_filter filter, 887 ions remain.
#> ℹ Parsing 887 peaks based on the sample group: Blanks.
#> ℹ Argument remove_ions is: TRUE.Removing peaks from Blanks.
#> ✔ 494 ions failed the Blanks filter, 393 ions remain.
#> ℹ Parsing 393 peaks for insource ions.
#> ✔ 43 ions failed the insource filter, 350 ions remain.



matched_data <- ms2_ms1_compare(mums2_example("full_mix_ms2_small.mgf"),
 filtered_data, 2, 6)
#> [1] "Reading: /home/runner/work/_temp/Library/mums2/extdata/full_mix_ms2_small.mgf ..."
#> Computing                                                    | 0%  ETA: -...Computing ■                                                  | 2%  ETA: ...Computing ■■                                                 | 4%  ETA: ...Computing ■■■                                                | 6%  ETA: ...Computing ■■■■                                               | 8%  ETA: ...Computing ■■■■■                                              | 10%  ETA: ...Computing ■■■■■■                                             | 12%  ETA: ...Computing ■■■■■■■                                            | 14%  ETA: ...Computing ■■■■■■■■                                           | 16%  ETA: ...Computing ■■■■■■■■■                                          | 18%  ETA: ...Computing ■■■■■■■■■■                                         | 20%  ETA: ...Computing ■■■■■■■■■■■                                        | 22%  ETA: ...Computing ■■■■■■■■■■■■                                       | 24%  ETA: ...Computing ■■■■■■■■■■■■■                                      | 26%  ETA: ...Computing ■■■■■■■■■■■■■■                                     | 28%  ETA: ...Computing ■■■■■■■■■■■■■■■                                    | 30%  ETA: ...Computing ■■■■■■■■■■■■■■■■                                   | 32%  ETA: ...Computing ■■■■■■■■■■■■■■■■■                                  | 34%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■                                 | 36%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■                                | 38%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■                               | 40%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■                              | 42%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■                             | 44%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■                            | 46%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■                           | 48%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■                          | 50%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■                         | 52%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■                        | 54%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■                       | 56%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                      | 58%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                     | 60%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                    | 62%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                   | 64%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                  | 66%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                 | 68%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■               | 72%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■              | 74%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■             | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■            | 78%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■           | 80%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■          | 82%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■         | 84%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■        | 86%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■       | 88%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■      | 90%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     | 92%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 94%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 96%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  | 98%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
#> [1] "1/350 peaks have an MS2 spectra."

dist_gnps <- dist_ms2(data = matched_data,
 cutoff = 0.3, precursor_threshold = 2,
 score_params = modified_cosine_params(0.5), min_peaks = 0)
#> Computing                                                    | 0%  ETA: -...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...

dist_entropy <- dist_ms2(data = matched_data,
 cutoff = 0.3, precursor_threshold = 2,
 score_params = spec_entropy_params(), min_peaks = 0)
#> Computing                                                    | 0%  ETA: -...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
```
