# Combine Databases

Add another database to your reference database.

## Usage

``` r
combined_reference_database(reference, other_reference)
```

## Arguments

- reference:

  reference database object.

- other_reference:

  your other reference database object

## Value

a `reference_database` that includes references from both reference
databases.

## Examples

``` r
reference <- read_msp(mums2_example("massbank_example_data.msp"))
#> Reading: /home/runner/work/_temp/Library/mums2/extdata/massbank_example_data.msp ...
reference2 <- read_msp(mums2_example("massbank_example_data.msp"))
#> Reading: /home/runner/work/_temp/Library/mums2/extdata/massbank_example_data.msp ...
combined_reference_database(reference, reference2)
#> You have 16136 references in this object.
```
