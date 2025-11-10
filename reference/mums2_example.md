# Get file paths for examples

mums2 contains a number of example files in the `inst/extdata`
directory. This function makes them accessible in documentation that
shows how file paths are used in function examples.

## Usage

``` r
mums2_example(file = NULL)
```

## Arguments

- file:

  Name of a file. If `NULL`, all examples files will be listed.

## Value

A file path to example data stored in the `inst/extdata` directory of
the package.

returns a `character` object

## Examples

``` r
mums2_example()
#>  [1] "PSU-MSMLS.msp"                 "PSUMSMLS_Adenine.csv"         
#>  [3] "demo_massdataset"              "full_mix_meta_data.csv"       
#>  [5] "full_mix_meta_data_small.csv"  "full_mix_ms2.mgf"             
#>  [7] "full_mix_ms2_small.mgf"        "full_mix_peak_table.csv"      
#>  [9] "full_mix_peak_table_small.csv" "ion_masses"                   

mums2_example("PSUMSMLS_Adenine.csv")
#> [1] "/home/runner/work/_temp/Library/mums2/extdata/PSUMSMLS_Adenine.csv"
```
