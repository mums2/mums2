# Rarefy MS1 Feature Table

`rarefy_ms()` performs a single subsampling of MS1 features in sample.
Feature intensities are subsampled to the supplied `size` and accounts
for intensity thresholds due to machine limits and background noise.
Specifically, features whose abundance falls below the `threshold` after
rarefying are removed. This allows for accurate representation of
samples at different dilutions regardless of the desired submsampling
`size`.

## Usage

``` r
rarefy_ms(
  community_object,
  size,
  threshold,
  number_of_threads = detectCores(),
  seed = 123
)
```

## Arguments

- community_object:

  A `community_object`

- size:

  The desired total sample intensity to subsample to.

- threshold:

  The individual feature threshold. Each subsampled feature must be \>=
  this value to be retained.

- number_of_threads:

  the amount of threads you want the calculation to use.

- seed:

  the RNG (random number generator) seed you would like to use.

## Value

A `external_pointer` that references a community matrix of rarefied
feature intensities.

returns a `matrix` object that contains your rarefied data.

## Examples

``` r
data <-
   import_all_data(peak_table =
                   mums2::mums2_example("full_mix_peak_table_small.csv"),
                   meta_data =
                   mums2::mums2_example("full_mix_meta_data_small.csv"),
                   format = "Metaboscape")
#> If peak table has corrupted compound names they will be converted to
#>       utf-8 and if there are any commas, they will be converted to periods(.).

filtered_data <- data |>
   filter_peak_table(filter_mispicked_ions_params()) |>
   filter_peak_table(filter_cv_params(cv_threshold = 0.2)) |>
   filter_peak_table(filter_group_params(group_threshold = 0.1,
                                             "Blanks")) |>
   filter_peak_table(filter_insource_ions_params())
#> ℹ Checking 1294 peaks for mispicked peaks.
#> ℹ Argument merge_peaks is: TRUE. Merging mispicked peaks with method sum.
#> ✔ 8 ions failed the mispicked filter, 1286 ions remain.
#> ℹ Parsing 1286 peaks for replicability across technical replicates.
#> ✔ 399 ions failed the cv_filter filter, 887 ions remain.
#> ℹ Parsing 887 peaks based on the sample group: Blanks.
#> ℹ Argument remove_ions is: TRUE.Removing peaks from Blanks.
#> ✔ 494 ions failed the Blanks filter, 393 ions remain.
#> ℹ Parsing 393 peaks for insource ions.
#> ✔ 43 ions failed the insource filter, 350 ions remain.

change_rt_to_seconds_or_minute(filtered_data, "minutes")
#> [1] "Changing rt values to minutes"
#> Key: <Compound, mz, kmd, RTINMINUTES>
#>      Compound        mz      kmd RTINMINUTES 50MeOH_10570 50MeOH_10571
#>        <char>     <num>    <num>       <num>        <num>        <num>
#>   1:     1000 749.51361 0.513606    3.622667            0            0
#>   2:     1002 734.46721 0.467206    3.725667            0            0
#>   3:     1003 749.51369 0.513686    3.730833            0            0
#>   4:     1006  42.03365 0.033646    3.802167            0            0
#>   5:     1007 137.07102 0.071016    3.807167            0            0
#>  ---                                                                  
#> 346:      989  83.04884 0.048836    3.307833            0            0
#> 347:      992 749.51424 0.514236    3.343167            0            0
#> 348:      996 716.45731 0.457306    3.447000            0            0
#> 349:      997 749.51417 0.514166    3.524000            0            0
#> 350:      999 558.36355 0.363546    3.556833            0            0
#>      50MeOH_10572 FullMix_10574 FullMix_10577 FullMix_10580    cor
#>             <num>         <num>         <num>         <num> <lgcl>
#>   1:            0     226277.39     187860.42     166603.72   TRUE
#>   2:            0    1285652.63    1244527.38    1222279.13   TRUE
#>   3:            0     325578.88     357518.06     316568.66   TRUE
#>   4:            0     270776.94     258514.84     251178.77   TRUE
#>   5:            0     116058.97     118762.07     102857.34   TRUE
#>  ---                                                              
#> 346:            0     198559.53     206690.77     217850.67   TRUE
#> 347:            0      86088.67      83416.42      65714.98   TRUE
#> 348:            0    1578886.38    1743622.13    1736709.88   TRUE
#> 349:            0     129910.10     137885.48     128256.25   TRUE
#> 350:            0     192658.81     218406.34     205511.38   TRUE

matched_data <- ms2_ms1_compare(mums2_example("full_mix_ms2_small.mgf"),
 filtered_data, 2, 6)
#> [1] "Reading: /home/runner/work/_temp/Library/mums2/extdata/full_mix_ms2_small.mgf ..."
#> Computing                                                    | 0%  ETA: -...Computing ■                                                  | 2%  ETA: ...Computing ■■                                                 | 4%  ETA: ...Computing ■■■                                                | 6%  ETA: ...Computing ■■■■                                               | 8%  ETA: ...Computing ■■■■■                                              | 10%  ETA: ...Computing ■■■■■■                                             | 12%  ETA: ...Computing ■■■■■■■                                            | 14%  ETA: ...Computing ■■■■■■■■                                           | 16%  ETA: ...Computing ■■■■■■■■■                                          | 18%  ETA: ...Computing ■■■■■■■■■■                                         | 20%  ETA: ...Computing ■■■■■■■■■■■                                        | 22%  ETA: ...Computing ■■■■■■■■■■■■                                       | 24%  ETA: ...Computing ■■■■■■■■■■■■■                                      | 26%  ETA: ...Computing ■■■■■■■■■■■■■■                                     | 28%  ETA: ...Computing ■■■■■■■■■■■■■■■                                    | 30%  ETA: ...Computing ■■■■■■■■■■■■■■■■                                   | 32%  ETA: ...Computing ■■■■■■■■■■■■■■■■■                                  | 34%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■                                 | 36%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■                                | 38%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■                               | 40%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■                              | 42%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■                             | 44%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■                            | 46%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■                           | 48%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■                          | 50%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■                         | 52%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■                        | 54%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■                       | 56%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                      | 58%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                     | 60%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                    | 62%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                   | 64%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                  | 66%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                 | 68%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■               | 72%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■              | 74%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■             | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■            | 78%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■           | 80%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■          | 82%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■         | 84%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■        | 86%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■       | 88%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■      | 90%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     | 92%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 94%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 96%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  | 98%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
#> [1] "214/350 peaks have an MS2 spectra."

dist <- dist_ms2(data = matched_data, cutoff = 0.3, precursor_thresh = 2,
 score_params = modified_cosine_params(0.5), min_peaks = 0)
#> Computing                                                    | 0%  ETA: -...Computing ■                                                  | 2%  ETA: ...Computing ■■                                                 | 4%  ETA: ...Computing ■■■                                                | 6%  ETA: ...Computing ■■■■                                               | 8%  ETA: ...Computing ■■■■■                                              | 10%  ETA: ...Computing ■■■■■■                                             | 12%  ETA: ...Computing ■■■■■■■                                            | 14%  ETA: ...Computing ■■■■■■■■                                           | 16%  ETA: ...Computing ■■■■■■■■■                                          | 18%  ETA: ...Computing ■■■■■■■■■■                                         | 20%  ETA: ...Computing ■■■■■■■■■■■                                        | 22%  ETA: ...Computing ■■■■■■■■■■■■                                       | 24%  ETA: ...Computing ■■■■■■■■■■■■■                                      | 26%  ETA: ...Computing ■■■■■■■■■■■■■■                                     | 28%  ETA: ...Computing ■■■■■■■■■■■■■■■                                    | 30%  ETA: ...Computing ■■■■■■■■■■■■■■■■                                   | 32%  ETA: ...Computing ■■■■■■■■■■■■■■■■■                                  | 34%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■                                 | 36%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■                                | 38%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■                               | 40%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■                              | 42%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■                             | 44%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■                            | 46%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■                           | 48%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■                          | 50%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■                         | 52%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■                        | 54%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■                       | 56%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                      | 58%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                     | 60%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                    | 62%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                   | 64%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                  | 66%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                 | 68%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■               | 72%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■              | 74%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■             | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■            | 78%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■           | 80%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■          | 82%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■         | 84%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■        | 86%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■       | 88%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■      | 90%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     | 92%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 94%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 96%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  | 98%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...

cluster_results <- cluster_data(distance_df = dist,
 ms2_match_data = matched_data, cutoff = 0.3, cluster_method = "opticlust")

community_object <- create_community_matrix_object(cluster_results)
rarefy_ms(community_object, 4000, 100)
#>               omu1 omu2 omu3 omu4 omu5 omu6 omu7 omu8 omu9 omu10 omu11 omu12
#> FullMix_10577  141    0    0    0    0    0  160    0    0     0     0     0
#> FullMix_10580  136    0    0    0    0    0  125    0    0     0     0     0
#> FullMix_10574  134  101    0    0    0    0  140    0    0     0     0     0
#> 50MeOH_10572     0    0    0    0    0    0    0    0 2604     0     0     0
#> 50MeOH_10571     0    0    0    0    0    0    0    0 2560     0     0     0
#> 50MeOH_10570     0    0    0    0    0    0    0    0 2419     0     0     0
#>               omu13 omu14 omu15 omu16 omu17 omu18 omu19 omu20 omu21 omu22 omu23
#> FullMix_10577     0     0     0     0     0     0     0     0     0     0     0
#> FullMix_10580     0     0     0     0     0     0     0     0     0     0     0
#> FullMix_10574     0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10572      0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10571      0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10570      0     0     0     0     0     0     0     0     0     0     0
#>               omu24 omu25 omu26 omu27 omu28 omu29 omu30 omu31 omu32 omu33 omu34
#> FullMix_10577     0     0     0     0     0     0     0     0     0     0     0
#> FullMix_10580     0     0     0     0     0     0     0     0     0     0     0
#> FullMix_10574     0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10572      0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10571      0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10570      0     0     0     0     0     0     0     0     0     0     0
#>               omu35 omu36 omu37 omu38 omu39 omu40 omu41 omu42 omu43 omu44 omu45
#> FullMix_10577     0     0     0     0     0     0     0     0     0     0     0
#> FullMix_10580     0     0     0     0     0     0     0     0     0     0     0
#> FullMix_10574     0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10572      0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10571      0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10570      0     0     0     0     0     0     0     0     0     0     0
#>               omu46 omu47 omu48 omu49 omu50 omu51 omu52 omu53 omu54 omu55 omu56
#> FullMix_10577     0     0     0     0     0     0     0     0     0     0     0
#> FullMix_10580     0     0     0     0     0     0     0     0     0     0     0
#> FullMix_10574     0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10572      0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10571      0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10570      0     0     0     0     0     0     0     0     0     0     0
#>               omu57 omu58 omu59 omu60 omu61 omu62 omu63 omu64 omu65 omu66 omu67
#> FullMix_10577     0     0     0     0     0     0     0     0     0     0     0
#> FullMix_10580     0     0     0     0     0     0     0     0     0     0     0
#> FullMix_10574     0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10572      0     0     0     0     0     0     0     0     0   349     0
#> 50MeOH_10571      0     0     0     0     0     0     0   111     0   382     0
#> 50MeOH_10570      0     0     0     0     0     0     0   100     0   376     0
#>               omu68 omu69 omu70 omu71 omu72 omu73 omu74 omu75 omu76 omu77 omu78
#> FullMix_10577     0     0     0     0     0     0     0     0     0     0     0
#> FullMix_10580     0     0     0     0     0     0     0     0     0     0     0
#> FullMix_10574     0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10572      0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10571      0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10570      0     0     0     0     0     0     0     0     0     0     0
#>               omu79 omu80 omu81 omu82 omu83 omu84 omu85 omu86 omu87 omu88 omu89
#> FullMix_10577     0     0     0     0     0     0     0     0     0     0     0
#> FullMix_10580     0     0     0     0     0     0     0     0     0     0     0
#> FullMix_10574     0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10572      0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10571      0     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10570      0     0     0     0     0     0     0     0     0     0     0
#>               omu90 omu91 omu92 omu93 omu94 omu95 omu96 omu97 omu98 omu99
#> FullMix_10577     0     0     0     0     0     0     0     0     0     0
#> FullMix_10580     0     0     0     0     0     0     0     0     0     0
#> FullMix_10574     0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10572      0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10571      0     0     0     0     0     0     0     0     0     0
#> 50MeOH_10570      0     0     0     0     0     0     0     0     0     0
#>               omu100 omu101 omu102 omu103 omu104 omu105 omu106 omu107 omu108
#> FullMix_10577      0      0      0      0      0      0      0      0      0
#> FullMix_10580      0      0      0      0      0      0      0      0      0
#> FullMix_10574      0      0      0      0      0      0      0      0      0
#> 50MeOH_10572       0      0      0      0      0      0      0      0      0
#> 50MeOH_10571       0      0      0      0      0      0      0      0      0
#> 50MeOH_10570       0      0      0      0      0      0      0      0      0
#>               omu109 omu110 omu111 omu112 omu113 omu114 omu115 omu116 omu117
#> FullMix_10577      0      0      0      0      0      0      0      0      0
#> FullMix_10580      0      0      0      0      0      0      0      0      0
#> FullMix_10574      0      0      0      0      0      0      0      0      0
#> 50MeOH_10572       0    424      0      0      0      0      0      0      0
#> 50MeOH_10571       0    393      0      0      0      0      0      0      0
#> 50MeOH_10570       0    518      0      0      0      0      0      0      0
#>               omu118 omu119 omu120 omu121 omu122 omu123 omu124 omu125 omu126
#> FullMix_10577      0    104      0      0      0      0      0      0   1339
#> FullMix_10580      0      0      0      0      0      0      0      0   1327
#> FullMix_10574      0      0      0      0      0      0      0      0   1264
#> 50MeOH_10572     166      0      0      0      0      0      0      0      0
#> 50MeOH_10571     160      0      0      0      0      0      0      0      0
#> 50MeOH_10570     111      0      0      0      0      0      0      0      0
#>               omu127 omu128 omu129 omu130 omu131 omu132 omu133 omu134 omu135
#> FullMix_10577    458    421    436      0    190      0      0      0      0
#> FullMix_10580    482    452    522      0    206      0      0      0      0
#> FullMix_10574    490    388    524      0    214      0      0      0      0
#> 50MeOH_10572       0      0      0      0      0    457      0      0      0
#> 50MeOH_10571       0      0      0      0      0    394      0      0      0
#> 50MeOH_10570       0      0      0      0      0    573      0      0      0
#>               omu136 omu137 omu138 omu139 omu140 omu141 omu142 omu143 omu144
#> FullMix_10577      0      0      0      0      0      0    173      0    301
#> FullMix_10580      0      0      0      0      0      0    156      0    277
#> FullMix_10574      0      0      0      0      0      0    160      0    327
#> 50MeOH_10572       0      0      0      0      0      0      0      0      0
#> 50MeOH_10571       0      0      0      0      0      0      0      0      0
#> 50MeOH_10570       0      0      0      0      0      0      0      0      0
#>               omu145 omu146 omu147 omu148 omu149 omu150 omu151 omu152 omu153
#> FullMix_10577      0      0      0      0      0      0      0    165    116
#> FullMix_10580      0      0      0      0      0      0      0    171    146
#> FullMix_10574      0      0      0      0      0      0      0    170    155
#> 50MeOH_10572       0      0      0      0      0      0      0      0      0
#> 50MeOH_10571       0      0      0      0      0      0      0      0      0
#> 50MeOH_10570       0      0      0      0      0      0      0      0      0
```
