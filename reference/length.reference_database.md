# Reference database length

returns the length of the database

## Usage

``` r
# S3 method for class 'reference_database'
length(x)
```

## Arguments

- x:

  reference database object.

## Value

returns the length of the regerence database

## Examples

``` r
reference <- read_msp(mums2_example("massbank_example_data.msp"))
#> Reading: /home/runner/work/_temp/Library/mums2/extdata/massbank_example_data.msp ...
length(reference)
#> [1] 8068
```
