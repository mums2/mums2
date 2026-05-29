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
  diversity_index = c("shannon", "simpson"),
  subsample = TRUE,
  number_of_threads = detectCores(),
  iterations = 100,
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

  the diversity index you wish to calculate diversity, the options are
  shannon, simpson, or richness. You may also compute many indexes at
  the same time using a vector (ie. c("shannon", "simpson")).

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
                   metadata =
                   mums2::mums2_example("boryillus_metadata.csv"),
                   format = "None")



matched_data <- ms2_ms1_compare(mums2_example("botryllus_v2.gnps.mgf"),
 data, 1, 6)
#> [1] "Reading: /home/runner/work/_temp/Library/mums2/extdata/botryllus_v2.gnps.mgf ..."
#> Computing                                                    | 0%  ETA: -...Computing ■                                                  | 2%  ETA: ...Computing ■■                                                 | 4%  ETA: ...Computing ■■■                                                | 6%  ETA: ...Computing ■■■■                                               | 8%  ETA: ...Computing ■■■■■                                              | 10%  ETA: ...Computing ■■■■■■                                             | 12%  ETA: ...Computing ■■■■■■■                                            | 14%  ETA: ...Computing ■■■■■■■■                                           | 16%  ETA: ...Computing ■■■■■■■■■                                          | 18%  ETA: ...Computing ■■■■■■■■■■                                         | 20%  ETA: ...Computing ■■■■■■■■■■■                                        | 22%  ETA: ...Computing ■■■■■■■■■■■■                                       | 24%  ETA: ...Computing ■■■■■■■■■■■■■                                      | 26%  ETA: ...Computing ■■■■■■■■■■■■■■                                     | 28%  ETA: ...Computing ■■■■■■■■■■■■■■■                                    | 30%  ETA: ...Computing ■■■■■■■■■■■■■■■■                                   | 32%  ETA: ...Computing ■■■■■■■■■■■■■■■■■                                  | 34%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■                                 | 36%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■                                | 38%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■                               | 40%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■                              | 42%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■                             | 44%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■                            | 46%  ETA: 1s ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■                           | 48%  ETA: 1s ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■                          | 50%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■                         | 52%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■                        | 54%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■                       | 56%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                      | 58%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                     | 60%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                    | 62%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                   | 64%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                  | 66%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                 | 68%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■               | 72%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■              | 74%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■             | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■            | 78%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■           | 80%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■          | 82%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■         | 84%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■        | 86%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■       | 88%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■      | 90%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     | 92%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 94%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 96%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  | 98%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
#> [1] "17/349 peaks have an MS2 spectra."

dist <- dist_ms2(data = matched_data, cutoff = 0.3, precursor_thresh = 2,
 score_params = modified_cosine_params(0.5), min_peaks = 0,
 number_of_threads = 2)
#> Computing                                                    | 0%  ETA: -...Computing ■■                                                 | 5%  ETA: ...Computing ■■■■■                                              | 11%  ETA: ...Computing ■■■■■■■■                                           | 17%  ETA: ...Computing ■■■■■■■■■■■                                        | 23%  ETA: ...Computing ■■■■■■■■■■■■■■                                     | 29%  ETA: ...Computing ■■■■■■■■■■■■■■■■■                                  | 35%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■                               | 41%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■                            | 47%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■                         | 52%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                      | 58%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                   | 64%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■             | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■          | 82%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■       | 88%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 94%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...

cluster_results <- cluster_data(distance_df = dist,
 ms2_match_data = matched_data,
 cutoff = 0.3, cluster_method = "opticlust")
#> Warning: [WARNING]: The mcc metric is not suitible for your data with a cutoff of 0.300000 using tptn instead.

community_object <- create_community_matrix_object(cluster_results)

alpha_summary(community_object = community_object, size = 400,
              threshold = 100,
              diversity_index = c("shannon", "simpson", "richness"),
              subsample = TRUE, iterations = 1, number_of_threads = 1)
#> Computing                                                    | 0%  ETA: -...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
#>                                               samples   shannon   simpson
#> 221012_DGM_Blank4_1_2_435   221012_DGM_Blank4_1_2_435 0.6901267 0.4997368
#> 221012_DGM_Blank4_1_1_434   221012_DGM_Blank4_1_1_434 0.0000000 0.0000000
#> 221012_DGM_MB1599_13_3_433 221012_DGM_MB1599_13_3_433 1.0308735 0.6446367
#> 221012_DGM_MB1599_13_1_431 221012_DGM_MB1599_13_1_431 1.3606209 0.7450190
#> 221012_DGM_MB1598_12_3_430 221012_DGM_MB1598_12_3_430 0.6539265 0.4812030
#> 221012_DGM_MB1598_12_1_428 221012_DGM_MB1598_12_1_428 1.0688209 0.6578568
#> 221012_DGM_MB1597_11_2_426 221012_DGM_MB1597_11_2_426 0.0000000 0.0000000
#> 221012_DGM_MB1595_10_3_424 221012_DGM_MB1595_10_3_424 0.6821825 0.4957268
#> 221012_DGM_MB1595_10_2_423 221012_DGM_MB1595_10_2_423 0.6931472 0.5012531
#> 221012_DGM_MB1597_11_3_427 221012_DGM_MB1597_11_3_427 0.0000000 0.0000000
#> 221012_DGM_MB1595_10_1_422 221012_DGM_MB1595_10_1_422 0.6901267 0.4997368
#> 221012_DGM_Blank2_1_1_404   221012_DGM_Blank2_1_1_404 0.5142872 0.4029306
#> 221012_DGM_Blank3_1_1_419   221012_DGM_Blank3_1_1_419 1.0844799 0.6633429
#> 221012_DGM_MB1590_5_3_403   221012_DGM_MB1590_5_3_403 0.6457396 0.4768800
#> 221012_DGM_Blank1_1_2_391   221012_DGM_Blank1_1_2_391 0.0000000 0.0000000
#> 221012_DGM_MB1599_13_2_432 221012_DGM_MB1599_13_2_432 1.3146591 0.7327544
#> 221012_DGM_MB1590_5_2_402   221012_DGM_MB1590_5_2_402 0.0000000 0.0000000
#> 221012_DGM_MB1597_11_1_425 221012_DGM_MB1597_11_1_425 0.0000000 0.0000000
#> 221012_DGM_Blank1_1_1_390   221012_DGM_Blank1_1_1_390 1.3293724 0.7366675
#> 221012_DGM_MB1590_5_1_401   221012_DGM_MB1590_5_1_401 0.4362592 0.3543390
#> 221012_DGM_MB1588_3_3_397   221012_DGM_MB1588_3_3_397 0.0000000 0.0000000
#> 221012_DGM_MB1588_3_1_395   221012_DGM_MB1588_3_1_395 0.4620097 0.3708023
#> 221012_DGM_Blank1_1_3_392   221012_DGM_Blank1_1_3_392 0.0000000 0.0000000
#> 221012_DGM_MB1592_7_2_411   221012_DGM_MB1592_7_2_411 0.4306658 0.3506380
#> 221012_DGM_MB1591_6_1_407   221012_DGM_MB1591_6_1_407 0.6393932 0.4735714
#> 221012_DGM_MB1589_4_1_398   221012_DGM_MB1589_4_1_398 0.4241572 0.3464192
#> 221012_DGM_Blank4_1_3_436   221012_DGM_Blank4_1_3_436 1.0325149 0.6452834
#> 221012_DGM_Blank3_1_3_421   221012_DGM_Blank3_1_3_421 0.5975411 0.4508573
#> 221012_DGM_MB1592_7_3_412   221012_DGM_MB1592_7_3_412 0.5956576 0.4499248
#> 221012_DGM_MB1593_8_2_414   221012_DGM_MB1593_8_2_414 1.3632266 0.7453416
#> 221012_DGM_MB1589_4_3_400   221012_DGM_MB1589_4_3_400 1.0166440 0.6394369
#> 221012_DGM_MB1598_12_2_429 221012_DGM_MB1598_12_2_429 0.9673828 0.6211535
#> 221012_DGM_MB1589_4_2_399   221012_DGM_MB1589_4_2_399 0.6123970 0.4590977
#> 221012_DGM_Blank2_1_2_405   221012_DGM_Blank2_1_2_405 0.0000000 0.0000000
#> 221012_DGM_MB1591_6_2_408   221012_DGM_MB1591_6_2_408 0.6925224 0.5009398
#> 221012_DGM_Blank2_1_3_406   221012_DGM_Blank2_1_3_406 0.5531454 0.4259137
#> 221012_DGM_MB1593_8_1_413   221012_DGM_MB1593_8_1_413 1.3025835 0.7295853
#> 221012_DGM_MB1591_6_3_409   221012_DGM_MB1591_6_3_409 0.6437620 0.4758772
#> 221012_DGM_MB1592_7_1_410   221012_DGM_MB1592_7_1_410 0.6629818 0.4859023
#> 221012_DGM_MB1588_3_2_396   221012_DGM_MB1588_3_2_396 0.5806894 0.4415915
#> 221012_DGM_MB1593_8_3_415   221012_DGM_MB1593_8_3_415 1.3498998 0.7420265
#> 221012_DGM_MB1594_9_1_416   221012_DGM_MB1594_9_1_416 1.3090645 0.7311946
#> 221012_DGM_MB1594_9_2_417   221012_DGM_MB1594_9_2_417 0.6227819 0.4647118
#> 221012_DGM_MB1594_9_3_418   221012_DGM_MB1594_9_3_418 1.0076959 0.6362609
#> 221012_DGM_Blank3_1_2_420   221012_DGM_Blank3_1_2_420 0.9013939 0.5951752
#>                            richness
#> 221012_DGM_Blank4_1_2_435         2
#> 221012_DGM_Blank4_1_1_434         1
#> 221012_DGM_MB1599_13_3_433        3
#> 221012_DGM_MB1599_13_1_431        4
#> 221012_DGM_MB1598_12_3_430        2
#> 221012_DGM_MB1598_12_1_428        3
#> 221012_DGM_MB1597_11_2_426        1
#> 221012_DGM_MB1595_10_3_424        2
#> 221012_DGM_MB1595_10_2_423        2
#> 221012_DGM_MB1597_11_3_427        1
#> 221012_DGM_MB1595_10_1_422        2
#> 221012_DGM_Blank2_1_1_404         2
#> 221012_DGM_Blank3_1_1_419         3
#> 221012_DGM_MB1590_5_3_403         2
#> 221012_DGM_Blank1_1_2_391         1
#> 221012_DGM_MB1599_13_2_432        4
#> 221012_DGM_MB1590_5_2_402         1
#> 221012_DGM_MB1597_11_1_425        1
#> 221012_DGM_Blank1_1_1_390         4
#> 221012_DGM_MB1590_5_1_401         2
#> 221012_DGM_MB1588_3_3_397         1
#> 221012_DGM_MB1588_3_1_395         2
#> 221012_DGM_Blank1_1_3_392         1
#> 221012_DGM_MB1592_7_2_411         2
#> 221012_DGM_MB1591_6_1_407         2
#> 221012_DGM_MB1589_4_1_398         2
#> 221012_DGM_Blank4_1_3_436         3
#> 221012_DGM_Blank3_1_3_421         2
#> 221012_DGM_MB1592_7_3_412         2
#> 221012_DGM_MB1593_8_2_414         4
#> 221012_DGM_MB1589_4_3_400         3
#> 221012_DGM_MB1598_12_2_429        3
#> 221012_DGM_MB1589_4_2_399         2
#> 221012_DGM_Blank2_1_2_405         1
#> 221012_DGM_MB1591_6_2_408         2
#> 221012_DGM_Blank2_1_3_406         2
#> 221012_DGM_MB1593_8_1_413         4
#> 221012_DGM_MB1591_6_3_409         2
#> 221012_DGM_MB1592_7_1_410         2
#> 221012_DGM_MB1588_3_2_396         2
#> 221012_DGM_MB1593_8_3_415         4
#> 221012_DGM_MB1594_9_1_416         4
#> 221012_DGM_MB1594_9_2_417         2
#> 221012_DGM_MB1594_9_3_418         3
#> 221012_DGM_Blank3_1_2_420         3
```
