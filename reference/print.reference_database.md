# Print reference

print reference objects.

## Usage

``` r
# S3 method for class 'reference_database'
print(x, ...)
```

## Arguments

- x:

  reference database object.

- ...:

  any extra print arguments you want to include.

## Value

prints customized message to the console

## Examples

``` r
reference <- read_msp(mums2_example("massbank_example_data.msp"))
#> Reading: /home/runner/work/_temp/Library/mums2/extdata/massbank_example_data.msp ...
print(reference)
#> You have 8068 references in this object.
```
