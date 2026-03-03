# Filter Peak Table

This function is a wrapper for all of mpactr's filter functions. When
called with a list of parameters that was generated from one of the
following functions, it will call the subsequent filter:
[`filter_mispicked_ions_params()`](https://www.mums2.org/mums2/reference/filter_mispicked_ions_params.md),
[`filter_group_params()`](https://www.mums2.org/mums2/reference/filter_group_params.md),
[`filter_cv_params()`](https://www.mums2.org/mums2/reference/filter_cv_params.md),
and
[`filter_insource_ions_params()`](https://www.mums2.org/mums2/reference/filter_insource_ions_params.md).
You can also find more information on these functions in`mpactr`
documentation.

## Usage

``` r
filter_peak_table(mpactr_object, params)

# S3 method for class 'filter_mispicked_ions'
filter_peak_table(mpactr_object, params)

# S3 method for class 'filter_group'
filter_peak_table(mpactr_object, params)

# S3 method for class 'filter_cv'
filter_peak_table(mpactr_object, params)

# S3 method for class 'filter_insource_ions'
filter_peak_table(mpactr_object, params)
```

## Arguments

- mpactr_object:

  the mpactr_object is an object generated from the
  [`import_all_data()`](https://www.mums2.org/mums2/reference/import_all_data.md)
  function. This is how we begin our pipeline.

- params:

  the list of arguments generated from calling one of these functions:
  [`filter_mispicked_ions_params()`](https://www.mums2.org/mums2/reference/filter_mispicked_ions_params.md),
  [`filter_group_params()`](https://www.mums2.org/mums2/reference/filter_group_params.md),
  [`filter_cv_params()`](https://www.mums2.org/mums2/reference/filter_cv_params.md),
  and
  [`filter_insource_ions_params()`](https://www.mums2.org/mums2/reference/filter_insource_ions_params.md).

## Value

a `mpactr` object that has been filter based on the supplied parameters.

## Examples

``` r
data <-
   import_all_data(peak_table =
                   mums2::mums2_example("botryllus_pt_small.csv"),
                   meta_data =
                   mums2::mums2_example("meta_data_boryillus.csv"),
                   format = "None")
#> If peak table has corrupted compound names they will be converted to
#>       utf-8 and if there are any commas, they will be converted to periods(.).
filtered_data <- data |>
   filter_peak_table(filter_mispicked_ions_params()) |>
   filter_peak_table(filter_cv_params(cv_threshold = 0.2)) |>
   filter_peak_table(filter_group_params(group_threshold = 0.1,
                                             "Blanks")) |>
   filter_peak_table(filter_insource_ions_params())
#> ℹ Checking 1500 peaks for mispicked peaks.
#> ℹ Argument merge_peaks is: TRUE. Merging mispicked peaks with method sum.
#> ✔ 50 ions failed the mispicked filter, 1450 ions remain.
#> ℹ Parsing 1450 peaks for replicability across technical replicates.
#> ✔ 329 ions failed the cv_filter filter, 1121 ions remain.
#> ℹ Parsing 1121 peaks based on the sample group: Blanks.
#> ℹ Argument remove_ions is: TRUE.Removing peaks from Blanks.
#> ✔ 342 ions failed the Blanks filter, 779 ions remain.
#> ℹ Parsing 779 peaks for insource ions.
#> ✔ 55 ions failed the insource filter, 724 ions remain.
```
