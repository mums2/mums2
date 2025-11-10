# Filter Group Parameters

Creates a list of filter group arguments for the
[`filter_peak_table()`](https://www.mums2.org/mums2/reference/filter_peak_table.md)
function

## Usage

``` r
filter_group_params(
  group_threshold = 0.01,
  group_to_remove,
  remove_ions = TRUE,
  copy_object = FALSE
)
```

## Arguments

- group_threshold:

  Relative abundance threshold at which to remove ions. Default = 0.01.

- group_to_remove:

  Biological group name to remove ions from.

- remove_ions:

  A `boolean` parameter. If `TRUE` failing ions will be removed from the
  peak table. Default = TRUE.

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
filter_group_params(group_to_remove = "blank")
#> $group_threshold
#> [1] 0.01
#> 
#> $group_to_remove
#> [1] "blank"
#> 
#> $remove_ions
#> [1] TRUE
#> 
#> $copy_object
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "filter_group"
```
