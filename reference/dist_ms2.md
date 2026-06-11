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

  the number of threads you wish to use for this calculation. Defaults
  to the number of threads on your computer.

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
                   mums2::mums2_example("botryllus_pt_small.csv"),
                   metadata =
                   mums2::mums2_example("boryillus_metadata.csv"),
                   format = "None")

matched_data <- ms2_ms1_compare(mums2_example("botryllus_v2.gnps.mgf"),
 data, 1, 6)
#> Reading: /home/runner/work/_temp/Library/mums2/extdata/botryllus_v2.gnps.mgf ...
#> Computing                                | 0%  ETA: -...Computing ■                              | 3%  ETA: ...Computing ■■                             | 6%  ETA: 13s ...Computing ■■■                            | 10%  ETA: 8s ...Computing ■■■■                           | 13%  ETA: 6s ...Computing ■■■■■                          | 16%  ETA: 4s ...Computing ■■■■■■                         | 20%  ETA: 3s ...Computing ■■■■■■■                        | 23%  ETA: 3s ...Computing ■■■■■■■■                       | 26%  ETA: 2s ...Computing ■■■■■■■■■                      | 30%  ETA: 2s ...Computing ■■■■■■■■■■                     | 33%  ETA: 1s ...Computing ■■■■■■■■■■■                    | 36%  ETA: 1s ...Computing ■■■■■■■■■■■■                   | 40%  ETA: 1s ...Computing ■■■■■■■■■■■■■                  | 43%  ETA: 1s ...Computing ■■■■■■■■■■■■■■                 | 46%  ETA: 1s ...Computing ■■■■■■■■■■■■■■■                | 50%  ETA: ...Computing ■■■■■■■■■■■■■■■■               | 53%  ETA: ...Computing ■■■■■■■■■■■■■■■■■              | 56%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■             | 60%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■            | 63%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■           | 66%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■          | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■         | 73%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■        | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■       | 80%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■      | 83%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■     | 86%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 90%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 93%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  | 96%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
#> 17/349 peaks have an MS2 spectra.

dist_gnps <- dist_ms2(data = matched_data,
 cutoff = 0.3, precursor_threshold = 2,
 score_params = modified_cosine_params(0.5), min_peaks = 0,
 number_of_threads = 2)
#> Computing                                | 0%  ETA: -...Computing ■                              | 5%  ETA: ...Computing ■■■                            | 11%  ETA: ...Computing ■■■■■                          | 17%  ETA: ...Computing ■■■■■■■                        | 23%  ETA: ...Computing ■■■■■■■■                       | 29%  ETA: ...Computing ■■■■■■■■■■                     | 35%  ETA: ...Computing ■■■■■■■■■■■■                   | 41%  ETA: ...Computing ■■■■■■■■■■■■■■                 | 47%  ETA: ...Computing ■■■■■■■■■■■■■■■                | 52%  ETA: ...Computing ■■■■■■■■■■■■■■■■■              | 58%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■            | 64%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■          | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■         | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■       | 82%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■     | 88%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 94%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...

dist_entropy <- dist_ms2(data = matched_data,
 cutoff = 0.3, precursor_threshold = 2,
 score_params = spec_entropy_params(), min_peaks = 0,
 number_of_threads = 2)
#> Computing                                | 0%  ETA: -...Computing ■                              | 5%  ETA: ...Computing ■■■                            | 11%  ETA: ...Computing ■■■■■                          | 17%  ETA: ...Computing ■■■■■■■                        | 23%  ETA: ...Computing ■■■■■■■■                       | 29%  ETA: ...Computing ■■■■■■■■■■                     | 35%  ETA: ...Computing ■■■■■■■■■■■■                   | 41%  ETA: ...Computing ■■■■■■■■■■■■■■                 | 47%  ETA: ...Computing ■■■■■■■■■■■■■■■                | 52%  ETA: ...Computing ■■■■■■■■■■■■■■■■■              | 58%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■            | 64%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■          | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■         | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■       | 82%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■     | 88%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 94%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
```
