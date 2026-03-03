# Alpha Diversity Summary

Alpha Diversity calculates the amount of diversity in a single sample.
We can conduct this analysis using your created community object. We
support the use of Shannon and Simpson diversity index.

## Usage

``` r
alpha_summary(
  community_object,
  size,
  threshold,
  diversity_index = "shannon",
  subsample = TRUE,
  number_of_threads = detectCores(),
  iterations = 1000,
  seed = 123
)
```

## Arguments

- community_object:

  the object created from the `create_community_object()` function.

- size:

  the size you wish to rarefy your diversity matrix to.

- threshold:

  the threshold you want your species to reach before it is included in
  the rarefaction sum.

- diversity_index:

  the diversity index you wish to calculate diversity, the two options
  are shannon or simpson.

- subsample:

  if true, we will rarefy the data before we run the diversity
  calculations. Default is TRUE.

- number_of_threads:

  the amount of threads you want the calculation to use.

- iterations:

  the amount of times you wish to run your calculation.

- seed:

  the RNG (random number generator) seed you would like to use.

## Value

a `data.frame` object that shows the dissimilarity in samples.

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

 change_rt_to_seconds_or_minute(filtered_data, "minutes")
#> [1] "Changing rt values to minutes"
#> Key: <Compound, mz, kmd, RTINMINUTES>
#>                    Compound        mz     kmd RTINMINUTES
#>                      <char>     <num>   <num>       <num>
#>   1: 1000.65345 Da 418.99 s 1001.6607 0.66073   0.1163333
#>   2: 1002.23833 Da 440.75 s 1003.2456 0.24560   0.1225000
#>   3: 1002.57692 Da 414.74 s 1003.5842 0.58419   0.1151667
#>   4: 1004.19672 Da 603.96 s 1005.2040 0.20400   0.1678333
#>   5: 1004.23580 Da 548.28 s 1005.2431 0.24307   0.1523333
#>  ---                                                     
#> 720:  980.84769 Da 515.58 s  981.8550 0.85497   0.1431667
#> 721:  984.23693 Da 604.23 s  985.2442 0.24420   0.1678333
#> 722:  984.65546 Da 460.05 s  985.6627 0.66274   0.1278333
#> 723:  997.30596 Da 605.90 s  998.3132 0.31323   0.1683333
#> 724:  998.28767 Da 556.84 s  981.2843 0.28429   0.1546667
#>      221012_DGM_Blank1_1_1_390 221012_DGM_Blank1_1_2_391
#>                          <num>                     <num>
#>   1:                         0                         0
#>   2:                         0                         0
#>   3:                         0                         0
#>   4:                         0                         0
#>   5:                         0                         0
#>  ---                                                    
#> 720:                         0                         0
#> 721:                         0                         0
#> 722:                         0                         0
#> 723:                         0                         0
#> 724:                         0                         0
#>      221012_DGM_Blank1_1_3_392 221012_DGM_Blank2_1_1_404
#>                          <num>                     <num>
#>   1:                         0                         0
#>   2:                         0                         0
#>   3:                         0                         0
#>   4:                         0                         0
#>   5:                         0                         0
#>  ---                                                    
#> 720:                         0                         0
#> 721:                         0                         0
#> 722:                         0                         0
#> 723:                         0                         0
#> 724:                         0                         0
#>      221012_DGM_Blank2_1_2_405 221012_DGM_Blank2_1_3_406
#>                          <num>                     <num>
#>   1:                         0                         0
#>   2:                         0                         0
#>   3:                         0                         0
#>   4:                         0                         0
#>   5:                         0                         0
#>  ---                                                    
#> 720:                         0                         0
#> 721:                         0                         0
#> 722:                         0                         0
#> 723:                         0                         0
#> 724:                         0                         0
#>      221012_DGM_Blank3_1_1_419 221012_DGM_Blank3_1_2_420
#>                          <num>                     <num>
#>   1:                         0                         0
#>   2:                         0                         0
#>   3:                         0                         0
#>   4:                         0                         0
#>   5:                         0                         0
#>  ---                                                    
#> 720:                         0                         0
#> 721:                         0                         0
#> 722:                         0                         0
#> 723:                         0                         0
#> 724:                         0                         0
#>      221012_DGM_Blank3_1_3_421 221012_DGM_Blank4_1_1_434
#>                          <num>                     <num>
#>   1:                         0                   1538.23
#>   2:                         0                      0.00
#>   3:                         0                      0.00
#>   4:                         0                      0.00
#>   5:                         0                      0.00
#>  ---                                                    
#> 720:                         0                      0.00
#> 721:                         0                      0.00
#> 722:                         0                      0.00
#> 723:                         0                      0.00
#> 724:                         0                      0.00
#>      221012_DGM_Blank4_1_2_435 221012_DGM_Blank4_1_3_436
#>                          <num>                     <num>
#>   1:                  1201.261                  1180.144
#>   2:                     0.000                     0.000
#>   3:                     0.000                     0.000
#>   4:                     0.000                     0.000
#>   5:                     0.000                     0.000
#>  ---                                                    
#> 720:                     0.000                     0.000
#> 721:                     0.000                     0.000
#> 722:                     0.000                     0.000
#> 723:                     0.000                     0.000
#> 724:                     0.000                     0.000
#>      221012_DGM_MB1588_3_1_395 221012_DGM_MB1588_3_2_396
#>                          <num>                     <num>
#>   1:                     0.000                     0.000
#>   2:                     0.000                     0.000
#>   3:                     0.000                     0.000
#>   4:                     0.000                     0.000
#>   5:                     0.000                     0.000
#>  ---                                                    
#> 720:                     0.000                     0.000
#> 721:                     0.000                     0.000
#> 722:                     0.000                     0.000
#> 723:                  6487.566                  5098.135
#> 724:                202552.359                148361.500
#>      221012_DGM_MB1588_3_3_397 221012_DGM_MB1589_4_1_398
#>                          <num>                     <num>
#>   1:                     0.000                         0
#>   2:                     0.000                         0
#>   3:                     0.000                         0
#>   4:                     0.000                         0
#>   5:                     0.000                         0
#>  ---                                                    
#> 720:                     0.000                         0
#> 721:                     0.000                         0
#> 722:                     0.000                         0
#> 723:                  7221.981                         0
#> 724:                155704.094                         0
#>      221012_DGM_MB1589_4_2_399 221012_DGM_MB1589_4_3_400
#>                          <num>                     <num>
#>   1:                         0                         0
#>   2:                         0                         0
#>   3:                         0                         0
#>   4:                         0                         0
#>   5:                         0                         0
#>  ---                                                    
#> 720:                         0                         0
#> 721:                         0                         0
#> 722:                         0                         0
#> 723:                         0                         0
#> 724:                         0                         0
#>      221012_DGM_MB1590_5_1_401 221012_DGM_MB1590_5_2_402
#>                          <num>                     <num>
#>   1:                     0.000                     0.000
#>   2:                     0.000                     0.000
#>   3:                     0.000                     0.000
#>   4:                     0.000                     0.000
#>   5:                  1977.438                  1953.256
#>  ---                                                    
#> 720:                  5517.501                  4006.738
#> 721:                 17435.523                 17933.676
#> 722:                     0.000                     0.000
#> 723:                     0.000                     0.000
#> 724:                     0.000                     0.000
#>      221012_DGM_MB1590_5_3_403 221012_DGM_MB1591_6_1_407
#>                          <num>                     <num>
#>   1:                     0.000                         0
#>   2:                     0.000                         0
#>   3:                     0.000                         0
#>   4:                     0.000                         0
#>   5:                  2209.528                         0
#>  ---                                                    
#> 720:                  4882.215                         0
#> 721:                 18599.828                         0
#> 722:                     0.000                         0
#> 723:                     0.000                         0
#> 724:                     0.000                         0
#>      221012_DGM_MB1591_6_2_408 221012_DGM_MB1591_6_3_409
#>                          <num>                     <num>
#>   1:                         0                         0
#>   2:                         0                         0
#>   3:                         0                         0
#>   4:                         0                         0
#>   5:                         0                         0
#>  ---                                                    
#> 720:                         0                         0
#> 721:                         0                         0
#> 722:                         0                         0
#> 723:                         0                         0
#> 724:                         0                         0
#>      221012_DGM_MB1592_7_1_410 221012_DGM_MB1592_7_2_411
#>                          <num>                     <num>
#>   1:                         0                         0
#>   2:                         0                         0
#>   3:                         0                         0
#>   4:                         0                         0
#>   5:                         0                         0
#>  ---                                                    
#> 720:                         0                         0
#> 721:                         0                         0
#> 722:                         0                         0
#> 723:                         0                         0
#> 724:                         0                         0
#>      221012_DGM_MB1592_7_3_412 221012_DGM_MB1593_8_1_413
#>                          <num>                     <num>
#>   1:                         0                     0.000
#>   2:                         0                     0.000
#>   3:                         0                     0.000
#>   4:                         0                  3522.783
#>   5:                         0                     0.000
#>  ---                                                    
#> 720:                         0                     0.000
#> 721:                         0                     0.000
#> 722:                         0                 14306.895
#> 723:                         0                     0.000
#> 724:                         0                     0.000
#>      221012_DGM_MB1593_8_2_414 221012_DGM_MB1593_8_3_415
#>                          <num>                     <num>
#>   1:                      0.00                     0.000
#>   2:                      0.00                     0.000
#>   3:                      0.00                     0.000
#>   4:                   5010.69                  3815.988
#>   5:                      0.00                     0.000
#>  ---                                                    
#> 720:                      0.00                     0.000
#> 721:                      0.00                     0.000
#> 722:                  10564.97                 14566.125
#> 723:                      0.00                     0.000
#> 724:                      0.00                     0.000
#>      221012_DGM_MB1594_9_1_416 221012_DGM_MB1594_9_2_417
#>                          <num>                     <num>
#>   1:                         0                         0
#>   2:                         0                         0
#>   3:                         0                         0
#>   4:                         0                         0
#>   5:                         0                         0
#>  ---                                                    
#> 720:                         0                         0
#> 721:                         0                         0
#> 722:                         0                         0
#> 723:                         0                         0
#> 724:                         0                         0
#>      221012_DGM_MB1594_9_3_418 221012_DGM_MB1595_10_1_422
#>                          <num>                      <num>
#>   1:                         0                   9122.671
#>   2:                         0                  35174.156
#>   3:                         0                      0.000
#>   4:                         0                   1316.851
#>   5:                         0                      0.000
#>  ---                                                     
#> 720:                         0                      0.000
#> 721:                         0                      0.000
#> 722:                         0                      0.000
#> 723:                         0                      0.000
#> 724:                         0                      0.000
#>      221012_DGM_MB1595_10_2_423 221012_DGM_MB1595_10_3_424
#>                           <num>                      <num>
#>   1:                   9405.939                  10668.905
#>   2:                  32907.184                  30415.164
#>   3:                      0.000                      0.000
#>   4:                   1761.997                   1803.242
#>   5:                      0.000                      0.000
#>  ---                                                      
#> 720:                      0.000                      0.000
#> 721:                      0.000                      0.000
#> 722:                      0.000                      0.000
#> 723:                      0.000                      0.000
#> 724:                      0.000                      0.000
#>      221012_DGM_MB1597_11_1_425 221012_DGM_MB1597_11_2_426
#>                           <num>                      <num>
#>   1:                      0.000                       0.00
#>   2:                      0.000                       0.00
#>   3:                      0.000                       0.00
#>   4:                      0.000                       0.00
#>   5:                      0.000                       0.00
#>  ---                                                      
#> 720:                      0.000                       0.00
#> 721:                      0.000                       0.00
#> 722:                   7123.614                    9147.79
#> 723:                      0.000                       0.00
#> 724:                      0.000                       0.00
#>      221012_DGM_MB1597_11_3_427 221012_DGM_MB1598_12_1_428
#>                           <num>                      <num>
#>   1:                      0.000                          0
#>   2:                      0.000                          0
#>   3:                      0.000                          0
#>   4:                      0.000                          0
#>   5:                      0.000                          0
#>  ---                                                      
#> 720:                      0.000                          0
#> 721:                      0.000                          0
#> 722:                   9354.261                          0
#> 723:                      0.000                          0
#> 724:                      0.000                          0
#>      221012_DGM_MB1598_12_2_429 221012_DGM_MB1598_12_3_430
#>                           <num>                      <num>
#>   1:                          0                          0
#>   2:                          0                          0
#>   3:                          0                          0
#>   4:                          0                          0
#>   5:                          0                          0
#>  ---                                                      
#> 720:                          0                          0
#> 721:                          0                          0
#> 722:                          0                          0
#> 723:                          0                          0
#> 724:                          0                          0
#>      221012_DGM_MB1599_13_1_431 221012_DGM_MB1599_13_2_432
#>                           <num>                      <num>
#>   1:                   10114.79                   10594.81
#>   2:                       0.00                       0.00
#>   3:                   15386.47                   15706.29
#>   4:                       0.00                       0.00
#>   5:                       0.00                       0.00
#>  ---                                                      
#> 720:                       0.00                       0.00
#> 721:                       0.00                       0.00
#> 722:                       0.00                       0.00
#> 723:                       0.00                       0.00
#> 724:                       0.00                       0.00
#>      221012_DGM_MB1599_13_3_433    cor
#>                           <num> <lgcl>
#>   1:                   10425.28   TRUE
#>   2:                       0.00   TRUE
#>   3:                   16543.69   TRUE
#>   4:                       0.00   TRUE
#>   5:                       0.00   TRUE
#>  ---                                  
#> 720:                       0.00   TRUE
#> 721:                       0.00   TRUE
#> 722:                       0.00   TRUE
#> 723:                       0.00   TRUE
#> 724:                       0.00   TRUE

matched_data <- ms2_ms1_compare(mums2_example("botryllus_v2.gnps.mgf"),
 filtered_data, 10, 6)
#> [1] "Reading: /home/runner/work/_temp/Library/mums2/extdata/botryllus_v2.gnps.mgf ..."
#> Computing                                                    | 0%  ETA: -...Computing ■                                                  | 2%  ETA: ...Computing ■■                                                 | 4%  ETA: ...Computing ■■■                                                | 6%  ETA: ...Computing ■■■■                                               | 8%  ETA: ...Computing ■■■■■                                              | 10%  ETA: ...Computing ■■■■■■                                             | 12%  ETA: ...Computing ■■■■■■■                                            | 14%  ETA: ...Computing ■■■■■■■■                                           | 16%  ETA: ...Computing ■■■■■■■■■                                          | 18%  ETA: ...Computing ■■■■■■■■■■                                         | 20%  ETA: ...Computing ■■■■■■■■■■■                                        | 22%  ETA: ...Computing ■■■■■■■■■■■■                                       | 24%  ETA: ...Computing ■■■■■■■■■■■■■                                      | 26%  ETA: ...Computing ■■■■■■■■■■■■■■                                     | 28%  ETA: ...Computing ■■■■■■■■■■■■■■■                                    | 30%  ETA: ...Computing ■■■■■■■■■■■■■■■■                                   | 32%  ETA: ...Computing ■■■■■■■■■■■■■■■■■                                  | 34%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■                                 | 36%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■                                | 38%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■                               | 40%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■                              | 42%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■                             | 44%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■                            | 46%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■                           | 48%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■                          | 50%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■                         | 52%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■                        | 54%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■                       | 56%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                      | 58%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                     | 60%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                    | 62%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                   | 64%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                  | 66%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                 | 68%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■               | 72%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■              | 74%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■             | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■            | 78%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■           | 80%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■          | 82%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■         | 84%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■        | 86%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■       | 88%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■      | 90%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     | 92%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 94%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 96%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  | 98%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
#> [1] "41/724 peaks have an MS2 spectra."
dist <- dist_ms2(data = matched_data, cutoff = 0.3, precursor_thresh = 2,
 score_params = modified_cosine_params(0.5), min_peaks = 0)
#> Computing                                                    | 0%  ETA: -...Computing ■                                                  | 2%  ETA: ...Computing ■■                                                 | 4%  ETA: ...Computing ■■■                                                | 7%  ETA: ...Computing ■■■■                                               | 9%  ETA: ...Computing ■■■■■■                                             | 12%  ETA: ...Computing ■■■■■■■                                            | 14%  ETA: ...Computing ■■■■■■■■                                           | 17%  ETA: ...Computing ■■■■■■■■■                                          | 19%  ETA: ...Computing ■■■■■■■■■■                                         | 21%  ETA: ...Computing ■■■■■■■■■■■■                                       | 24%  ETA: ...Computing ■■■■■■■■■■■■■                                      | 26%  ETA: ...Computing ■■■■■■■■■■■■■■                                     | 29%  ETA: ...Computing ■■■■■■■■■■■■■■■                                    | 31%  ETA: ...Computing ■■■■■■■■■■■■■■■■■                                  | 34%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■                                 | 36%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■                                | 39%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■                               | 41%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■                              | 43%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■                            | 46%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■                           | 48%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■                          | 51%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■                         | 53%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■                       | 56%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                      | 58%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                     | 60%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                    | 63%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                   | 65%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                 | 68%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■               | 73%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■              | 75%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■            | 78%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■           | 80%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■          | 82%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■         | 85%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■        | 87%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■      | 90%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     | 92%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 95%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 97%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...

cluster_results <- cluster_data(distance_df = dist,
 ms2_match_data = matched_data,
 cutoff = 0.3, cluster_method = "opticlust")

community_object <- create_community_matrix_object(cluster_results)

alpha_summary(community_object, 4000, 100, "shannon", TRUE, iterations = 1)
#> Computing                                                    | 0%  ETA: -...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
#>      221012_DGM_Blank4_1_2_435 221012_DGM_Blank4_1_1_434
#> [1,]                 0.5215757                  0.526344
#>      221012_DGM_MB1599_13_3_433 221012_DGM_MB1599_13_1_431
#> [1,]                   1.891786                   1.867916
#>      221012_DGM_MB1598_12_3_430 221012_DGM_MB1598_12_1_428
#> [1,]                  0.1062819                  0.1181969
#>      221012_DGM_MB1597_11_2_426 221012_DGM_MB1595_10_3_424
#> [1,]                   1.350517                  0.8803178
#>      221012_DGM_MB1595_10_2_423 221012_DGM_MB1597_11_3_427
#> [1,]                   1.024649                   1.297108
#>      221012_DGM_MB1595_10_1_422 221012_DGM_Blank2_1_1_404
#> [1,]                  0.9095848                 0.3319118
#>      221012_DGM_Blank3_1_1_419 221012_DGM_MB1590_5_3_403
#> [1,]                 0.5471333                  1.006397
#>      221012_DGM_Blank1_1_2_391 221012_DGM_MB1599_13_2_432
#> [1,]                 0.8877297                   1.865068
#>      221012_DGM_MB1590_5_2_402 221012_DGM_MB1597_11_1_425
#> [1,]                 0.9873934                   1.377731
#>      221012_DGM_Blank1_1_1_390 221012_DGM_MB1590_5_1_401
#> [1,]                  1.008975                 0.9561101
#>      221012_DGM_MB1588_3_3_397 221012_DGM_MB1588_3_1_395
#> [1,]                 0.8172633                 0.8757524
#>      221012_DGM_Blank1_1_3_392 221012_DGM_MB1592_7_2_411
#> [1,]                  1.023802                 0.8823921
#>      221012_DGM_MB1591_6_1_407 221012_DGM_MB1589_4_1_398
#> [1,]                  1.139187                 0.5805109
#>      221012_DGM_Blank4_1_3_436 221012_DGM_Blank3_1_3_421
#> [1,]                 0.4856461                 0.5353043
#>      221012_DGM_MB1592_7_3_412 221012_DGM_MB1593_8_2_414
#> [1,]                 0.8636725                  1.064158
#>      221012_DGM_MB1589_4_3_400 221012_DGM_MB1598_12_2_429
#> [1,]                 0.5321427                  0.1167086
#>      221012_DGM_MB1589_4_2_399 221012_DGM_Blank2_1_2_405
#> [1,]                 0.6023095                 0.3810227
#>      221012_DGM_MB1591_6_2_408 221012_DGM_Blank2_1_3_406
#> [1,]                   1.02966                 0.4227016
#>      221012_DGM_MB1593_8_1_413 221012_DGM_MB1591_6_3_409
#> [1,]                  1.001129                  1.010464
#>      221012_DGM_MB1592_7_1_410 221012_DGM_MB1588_3_2_396
#> [1,]                 0.8478471                 0.8442407
#>      221012_DGM_MB1593_8_3_415 221012_DGM_MB1594_9_1_416
#> [1,]                  1.100267                   1.58547
#>      221012_DGM_MB1594_9_2_417 221012_DGM_MB1594_9_3_418
#> [1,]                  1.566287                  1.567183
#>      221012_DGM_Blank3_1_2_420
#> [1,]                 0.4787449
```
