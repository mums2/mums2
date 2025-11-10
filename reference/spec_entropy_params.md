# Entropy similarity between two MS/MS spectra

Calculate spectral entropy similarity between two MS2 spectra

## Usage

``` r
spec_entropy_params(
  ms2_tolerance_in_da = 0.02,
  ms2_tolerance_in_ppm = -1,
  clean_spectra = TRUE,
  min_mz = 0,
  max_mz = 1000,
  noise_threshold = 0.01,
  max_peak_num = 100,
  weighted = TRUE
)
```

## Arguments

- ms2_tolerance_in_da:

  MS2 peak tolerance in Da, set to -1 to disable. Defaults to `0.02`.

- ms2_tolerance_in_ppm:

  MS2 peak tolerance in ppm, set to -1 to disable. Defaults to `-1`.

- clean_spectra:

  Either `TRUE` or `FALSE` to clean the spectra prior to calculating
  similarity. See
  [`msentropy::clean_spectrum`](https://rdrr.io/pkg/msentropy/man/clean_spectrum.html)
  for more information. Defaults to `TRUE`.

- min_mz:

  `numeric`, minimum mz to keep, set to -1 to disable. Defaults to `0`.

- max_mz:

  `numeric`, maximum mz to keep, set to -1 to disable. Defaults to
  `1000`.

- noise_threshold:

  Background intensity threshold, all peaks with intensity \<
  noise_threshold \* max_intensity are removed. Set to -1 to disable.
  Defaults to `0.01`.

- max_peak_num:

  `numeric`, maximum number of peaks to keep for score calculation. Set
  to -1 to disable. Defaults to `100`.

- weighted:

  `logical` whether weighted or unweighted entropy similarity will be
  calculated. Defaults to `TRUE`.

## Value

A parameters list for similarity scoring method "spectral_entropy"

## Details

`spec_entropy_params()` will initiate spectral entropy similarity
scoring via the `msentropy` package (Li et al. 2021). For more
information about parameters, see
[`msentropy::msentropy_similarity()`](https://rdrr.io/pkg/msentropy/man/msentropy_similarity.html).

## References

Li, Y., Kind, T., Folz, J. et al. Spectral entropy outperforms MS/MS dot
product similarity for small-molecule compound identification, Nat
Methods 18, 1524–1531 (2021). https://doi.org/10.1038/s41592-021-01331-z

## Examples

``` r
spec_entropy_params()
#> $ms2_tolerance_in_da
#> [1] 0.02
#> 
#> $ms2_tolerance_in_ppm
#> [1] -1
#> 
#> $clean_spectra
#> [1] TRUE
#> 
#> $min_mz
#> [1] 0
#> 
#> $max_mz
#> [1] 1000
#> 
#> $noise_threshold
#> [1] 0.01
#> 
#> $max_peak_num
#> [1] 100
#> 
#> $weighted
#> [1] TRUE
#> 
#> $method
#> [1] "entropy"
#> 
```
