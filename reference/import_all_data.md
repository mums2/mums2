# Import all data

This function is a wrapper for the mpactr import_data function. It will
import your peak table and meta data and create a mpactr_object.

## Usage

``` r
import_all_data(peak_table, metadata, format)
```

## Arguments

- peak_table:

  The file path to your feature table file.

- metadata:

  The file path to your metadata file or `data.frame`.

- format:

  The expected exported type of your peak table, can be one of
  "Progenesis", "Metaboscape", "None".

## Value

a `mpactr` object.

## Examples

``` r
data <-
   import_all_data(peak_table =
                   mums2::mums2_example("botryllus_pt_small.csv"),
                   metadata =
                   mums2::mums2_example("boryillus_metadata.csv"),
                   format = "None")
```
