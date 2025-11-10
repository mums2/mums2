# Filter Insource Ions Parameters

Creates a list of filter insource ions arguments for the
[`filter_peak_table()`](https://www.mums2.org/mums2/reference/filter_peak_table.md)
function

## Usage

``` r
filter_insource_ions_params(cluster_threshold = 0.95, copy_object = FALSE)
```

## Arguments

- cluster_threshold:

  Cluster threshold for ion deconvolution. Default = 0.95.

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
filter_insource_ions_params()
#> $cluster_threshold
#> [1] 0.95
#> 
#> $copy_object
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "filter_insource_ions"
```
