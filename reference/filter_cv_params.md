# Filter Cv Parameters

Creates a list of filter cv arguments for the
[`filter_peak_table()`](https://www.mums2.org/mums2/reference/filter_peak_table.md)
function

## Usage

``` r
filter_cv_params(cv_threshold = NULL, copy_object = FALSE)
```

## Arguments

- cv_threshold:

  Coefficient of variation threshold. A lower cv_threshold will result
  in more stringent filtering and higher reproducibility. Recommended
  values between 0.2 - 0.5.

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
filter_cv_params(0.2)
#> $cv_threshold
#> [1] 0.2
#> 
#> $copy_object
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "filter_cv"
filter_cv_params(0.2)
#> $cv_threshold
#> [1] 0.2
#> 
#> $copy_object
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "filter_cv"
```
