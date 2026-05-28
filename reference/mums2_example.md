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
#> [1] "MSMS-Neg-Respect.msp"               "boryillus_metadata.csv"            
#> [3] "boryillus_peaktable.csv"            "botryllus_pt_small.csv"            
#> [5] "botryllus_v2.gnps.mgf"              "massbank_example_data.msp"         
#> [7] "massbank_example_data_negative.msp"

mums2_example("massbank_example_data.msp")
#> [1] "/home/runner/work/_temp/Library/mums2/extdata/massbank_example_data.msp"
```
