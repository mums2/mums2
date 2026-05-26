# Get Samples

Returns a list of your samples found in the metadata file.

## Usage

``` r
get_samples(mass_data)
```

## Arguments

- mass_data:

  The object generated from
  [`ms2_ms1_compare()`](https://www.mums2.org/mums2/reference/ms2_ms1_compare.md).

## Value

a `character` vector contain all of your samples.

## Examples

``` r
data <-
   import_all_data(peak_table =
                   mums2::mums2_example("botryllus_pt_small.csv"),
                   metadata =
                   mums2::mums2_example("boryillus_metadata.csv"),
                   format = "None")


matched_data <- ms2_ms1_compare(mums2_example("botryllus_v2.gnps.mgf"),
 data, 1, 6)
#> [1] "Reading: /home/runner/work/_temp/Library/mums2/extdata/botryllus_v2.gnps.mgf ..."
#> Computing                                                    | 0%  ETA: -...Computing ■                                                  | 2%  ETA: ...Computing ■■                                                 | 4%  ETA: ...Computing ■■■                                                | 6%  ETA: ...Computing ■■■■                                               | 8%  ETA: ...Computing ■■■■■                                              | 10%  ETA: ...Computing ■■■■■■                                             | 12%  ETA: ...Computing ■■■■■■■                                            | 14%  ETA: ...Computing ■■■■■■■■                                           | 16%  ETA: ...Computing ■■■■■■■■■                                          | 18%  ETA: ...Computing ■■■■■■■■■■                                         | 20%  ETA: ...Computing ■■■■■■■■■■■                                        | 22%  ETA: ...Computing ■■■■■■■■■■■■                                       | 24%  ETA: ...Computing ■■■■■■■■■■■■■                                      | 26%  ETA: ...Computing ■■■■■■■■■■■■■■                                     | 28%  ETA: ...Computing ■■■■■■■■■■■■■■■                                    | 30%  ETA: ...Computing ■■■■■■■■■■■■■■■■                                   | 32%  ETA: ...Computing ■■■■■■■■■■■■■■■■■                                  | 34%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■                                 | 36%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■                                | 38%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■                               | 40%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■                              | 42%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■                             | 44%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■                            | 46%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■                           | 48%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■                          | 50%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■                         | 52%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■                        | 54%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■                       | 56%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                      | 58%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                     | 60%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                    | 62%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                   | 64%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                  | 66%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                 | 68%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■               | 72%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■              | 74%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■             | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■            | 78%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■           | 80%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■          | 82%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■         | 84%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■        | 86%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■       | 88%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■      | 90%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     | 92%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 94%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 96%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  | 98%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
#> [1] "17/349 peaks have an MS2 spectra."

get_samples(matched_data)
#>  [1] "221012_DGM_Blank1_1_1_390"  "221012_DGM_Blank1_1_2_391" 
#>  [3] "221012_DGM_Blank1_1_3_392"  "221012_DGM_MB1588_3_1_395" 
#>  [5] "221012_DGM_MB1588_3_2_396"  "221012_DGM_MB1588_3_3_397" 
#>  [7] "221012_DGM_MB1589_4_1_398"  "221012_DGM_MB1589_4_2_399" 
#>  [9] "221012_DGM_MB1589_4_3_400"  "221012_DGM_MB1590_5_1_401" 
#> [11] "221012_DGM_MB1590_5_2_402"  "221012_DGM_MB1590_5_3_403" 
#> [13] "221012_DGM_Blank2_1_1_404"  "221012_DGM_Blank2_1_2_405" 
#> [15] "221012_DGM_Blank2_1_3_406"  "221012_DGM_MB1591_6_1_407" 
#> [17] "221012_DGM_MB1591_6_2_408"  "221012_DGM_MB1591_6_3_409" 
#> [19] "221012_DGM_MB1592_7_1_410"  "221012_DGM_MB1592_7_2_411" 
#> [21] "221012_DGM_MB1592_7_3_412"  "221012_DGM_MB1593_8_1_413" 
#> [23] "221012_DGM_MB1593_8_2_414"  "221012_DGM_MB1593_8_3_415" 
#> [25] "221012_DGM_MB1594_9_1_416"  "221012_DGM_MB1594_9_2_417" 
#> [27] "221012_DGM_MB1594_9_3_418"  "221012_DGM_Blank3_1_1_419" 
#> [29] "221012_DGM_Blank3_1_2_420"  "221012_DGM_Blank3_1_3_421" 
#> [31] "221012_DGM_MB1595_10_1_422" "221012_DGM_MB1595_10_2_423"
#> [33] "221012_DGM_MB1595_10_3_424" "221012_DGM_MB1597_11_1_425"
#> [35] "221012_DGM_MB1597_11_2_426" "221012_DGM_MB1597_11_3_427"
#> [37] "221012_DGM_MB1598_12_1_428" "221012_DGM_MB1598_12_2_429"
#> [39] "221012_DGM_MB1598_12_3_430" "221012_DGM_MB1599_13_1_431"
#> [41] "221012_DGM_MB1599_13_2_432" "221012_DGM_MB1599_13_3_433"
#> [43] "221012_DGM_Blank4_1_1_434"  "221012_DGM_Blank4_1_2_435" 
#> [45] "221012_DGM_Blank4_1_3_436" 
```
