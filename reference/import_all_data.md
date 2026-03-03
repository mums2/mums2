# Import all data

This function is a wrapper for the mpactr import_data function. It will
import your peak table and meta data and create a mpactr_object.

## Usage

``` r
import_all_data(peak_table, meta_data, format)
```

## Arguments

- peak_table:

  The file path to your feature table file.

- meta_data:

  The file path to your meta_data file or `data.frame`.

- format:

  The expected exported type of your peak table, can be one of
  "Progenesis", "None", "None".

## Value

a `mpactr` object.

## Examples

``` r
data <-
   import_all_data(peak_table =
                   mums2::mums2_example("botryllus_pt_small.csv"),
                   meta_data =
                   mums2::mums2_example("meta_data_boryillus.csv"),
                   format = "None")
#> If peak table has corrupted compound names they will be converted to
#>       utf-8 and if there are any commas, they will be converted to periods(.).
```
