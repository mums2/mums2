# Convert Samples to Group Averages

To account for users measuring there data in triplicates or other forms
of measurement, we have implemented a function that can transform your
matched data object to use group averages instead of each sample
individually.

## Usage

``` r
convert_to_group_averages(matched_data, mpactr_object)
```

## Arguments

- matched_data:

  your mass data set object generated from
  [`ms2_ms1_compare()`](https://www.mums2.org/mums2/reference/ms2_ms1_compare.md).

- mpactr_object:

  The object created from
  [`import_all_data()`](https://www.mums2.org/mums2/reference/import_all_data.md).

## Value

a `mass_data` object using group averages

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

matched_data_avg <- convert_to_group_averages(matched_data,
                                              data)
```
