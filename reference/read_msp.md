# Create Reference Database

Creates a reference database by reading a download msp file. These files
can be downloaded from sites like
https://systemsomicslab.github.io/compms/msdial/main.html#MSP or
https://mona.fiehnlab.ucdavis.edu/downloads

## Usage

``` r
read_msp(msp_file)
```

## Arguments

- msp_file:

  the file path of your msp file

## Value

a `reference_database` object.

## Examples

``` r
read_msp(mums2_example("massbank_example_data.msp"))
#> Reading: /home/runner/work/_temp/Library/mums2/extdata/massbank_example_data.msp ...
#> You have 8068 references in this object.
```
