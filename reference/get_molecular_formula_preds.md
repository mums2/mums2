# Get Molecular Formula Predictions

Returns all of the molecular formula predictions.

## Usage

``` r
get_molecular_formula_preds(mass_data)
```

## Arguments

- mass_data:

  The object generated from
  [`ms2_ms1_compare()`](https://www.mums2.org/mums2/reference/ms2_ms1_compare.md).

## Value

a `character` vector contain all of your predicted molecular formulas.

## Examples

``` r
data <-
   import_all_data(peak_table =
                   mums2::mums2_example("botryllus_pt_small.csv"),
                   metadata =
                   mums2::mums2_example("boryillus_metadata.csv"),
                   format = "None")


matched_data <- ms2_ms1_compare(mums2_example("botryllus_v2.gnps.mgf"),
 data, 0.1, 1)
#> Reading: /home/runner/work/_temp/Library/mums2/extdata/botryllus_v2.gnps.mgf ...
#> 1/349 peaks have an MS2 spectra.

matched_data <- compute_molecular_formulas(matched_data,
                                           number_of_threads = 2)
#> Calculating potential molecular formulas...
#> Calculating fragmentation trees...
#> 1/1 chemical formulas were predicted

get_molecular_formula_preds(matched_data)
#> [1] "C17H57N13O17P2S3"
```
