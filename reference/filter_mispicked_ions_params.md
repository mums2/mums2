# Filter Mispicked Ions Parameters

Creates a list of filter mispicked ions arguments for the
[`filter_peak_table()`](https://www.mums2.org/mums2/reference/filter_peak_table.md)
function

## Usage

``` r
filter_mispicked_ions_params(
  ringwin = 0.5,
  isowin = 0.01,
  trwin = 0.005,
  max_iso_shift = 3,
  merge_peaks = TRUE,
  merge_method = "sum",
  copy_object = FALSE
)
```

## Arguments

- ringwin:

  Ringing mass window or detector saturation mass window. Default = 0.5
  atomic mass units (AMU).

- isowin:

  Isotopic mass window. Default = 0.01 AMU.

- trwin:

  A `numeric` denoting the retention time threshold for assessing if
  ions should be merged. Default = 0.005.

- max_iso_shift:

  A `numeric`. Default = 3.

- merge_peaks:

  A `boolean` parameter to determine if peaks found to belong to the
  same ion should be merged in the feature table.

- merge_method:

  If merge_peaks is TRUE, a method for how similar peaks should be
  merged. Can be one of "sum".

- copy_object:

  A `boolean` parameter that allows users to return a copied object
  instead of modifying the object.

## Value

a `list` object of arguments needed to call the given mpactr function
when supplied to the
[`filter_peak_table()`](https://www.mums2.org/mums2/reference/filter_peak_table.md)
wrapper function.

## Examples

``` r
filter_mispicked_ions_params()
#> $ringwin
#> [1] 0.5
#> 
#> $isowin
#> [1] 0.01
#> 
#> $trwin
#> [1] 0.005
#> 
#> $max_iso_shift
#> [1] 3
#> 
#> $merge_peaks
#> [1] TRUE
#> 
#> $merge_method
#> [1] "sum"
#> 
#> $copy_object
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "filter_mispicked_ions"
```
