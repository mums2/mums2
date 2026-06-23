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
  number_of_threads = 1
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

  the number of threads you wish to use for this calculation.

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
#> 17/349 peaks have an MS2 spectra.

dist_gnps <- dist_ms2(data = matched_data,
 cutoff = 0.3, precursor_threshold = 2,
 score_params = modified_cosine_params(0.5), min_peaks = 0,
 number_of_threads = 2)

dist_entropy <- dist_ms2(data = matched_data,
 cutoff = 0.3, precursor_threshold = 2,
 score_params = spec_entropy_params(), min_peaks = 0,
 number_of_threads = 2)
```
