# Compute Molecular formula

de novo algorithm for computing molecular formulas. Using fragmentation
trees we are able to generate a resultant molecular formula. To ensure
efficient we are using a greedy heurstic to generate the resultant
formula. Although this may not always result in the correct prediction,
it allows us to efficiently calculate a multitudeof chemical formulas.

## Usage

``` r
compute_molecular_formulas(
  mass_data,
  parent_ppm = 3,
  num_threads = detectCores()
)
```

## Arguments

- mass_data:

  your mass_data object generated from
  [`ms2_ms1_compare()`](https://www.mums2.org/mums2/reference/ms2_ms1_compare.md)

- parent_ppm:

  the ppm you wish to generate the candidate molecular formulas.

- num_threads:

  the amount of threads we the algorithm will use.

## Value

your mass_data object with an additional `character` vector of all the
predicted formulas.

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


matched_data <- ms2_ms1_compare(mums2_example("full_mix_ms2_small.mgf"),
 filtered_data, 2, 6)
#> [1] "Reading: /home/runner/work/_temp/Library/mums2/extdata/full_mix_ms2_small.mgf ..."
#> Computing                                                    | 0%  ETA: -...Computing ■                                                  | 2%  ETA: ...Computing ■■                                                 | 4%  ETA: ...Computing ■■■                                                | 6%  ETA: ...Computing ■■■■                                               | 8%  ETA: ...Computing ■■■■■                                              | 10%  ETA: ...Computing ■■■■■■                                             | 12%  ETA: ...Computing ■■■■■■■                                            | 14%  ETA: ...Computing ■■■■■■■■                                           | 16%  ETA: ...Computing ■■■■■■■■■                                          | 18%  ETA: ...Computing ■■■■■■■■■■                                         | 20%  ETA: ...Computing ■■■■■■■■■■■                                        | 22%  ETA: ...Computing ■■■■■■■■■■■■                                       | 24%  ETA: ...Computing ■■■■■■■■■■■■■                                      | 26%  ETA: ...Computing ■■■■■■■■■■■■■■                                     | 28%  ETA: ...Computing ■■■■■■■■■■■■■■■                                    | 30%  ETA: ...Computing ■■■■■■■■■■■■■■■■                                   | 32%  ETA: ...Computing ■■■■■■■■■■■■■■■■■                                  | 34%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■                                 | 36%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■                                | 38%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■                               | 40%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■                              | 42%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■                             | 44%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■                            | 46%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■                           | 48%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■                          | 50%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■                         | 52%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■                        | 54%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■                       | 56%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                      | 58%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                     | 60%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                    | 62%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                   | 64%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                  | 66%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                 | 68%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■               | 72%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■              | 74%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■             | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■            | 78%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■           | 80%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■          | 82%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■         | 84%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■        | 86%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■       | 88%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■      | 90%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     | 92%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 94%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 96%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  | 98%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
#> [1] "1/350 peaks have an MS2 spectra."
compute_molecular_formulas(matched_data)
#> Computing                                                    | 100%  ETA: -...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
#> 1/1 chemical formulas were predicted
#> $ms2_matches
#>         mz   rt ms1_compound_id spectra_index   ms2_spectrum_id
#> 1 239.1618 1.87             718             1 mz239.16175rt1.87
#> 
#> $peak_data
#> $peak_data[[1]]
#> $peak_data[[1]]$mz
#>  [1]  38.01506  38.96281  39.02327  41.03875  42.03413  42.03972  43.01874
#>  [8]  43.05478  44.04912  44.99170  45.03389  47.99072  52.94750  53.93871
#> [15]  55.01857  55.05485  55.93472  56.94236  56.96490  57.07004  57.93536
#> [22]  58.99330  59.04886  59.93003  61.01041  69.03312  69.06974  72.93790
#> [29]  73.03217  74.00629  97.00795 105.03361 105.95877 120.01301 122.02762
#> [36] 129.01374 132.02772 133.90649 136.03776 137.96378 138.05578 138.96609
#> [43] 140.97971 145.93059 151.01753 153.06692 156.95778 160.91077 164.02204
#> [50] 168.92493 170.99648 174.02347 174.96760 174.98656 179.97517 181.95017
#> [57] 181.95569 184.96368 186.93666 192.98391 193.93319 195.96563 196.10127
#> [64] 196.95447 196.97234 197.04933 197.98410 203.01392 204.01280 204.94744
#> [71] 206.93863 212.96473 214.91363 215.99458 217.01871 219.00773 221.02622
#> [78] 222.08246 223.17514 239.05184 239.08599 239.09347 261.89221
#> 
#> $peak_data[[1]]$intensity
#>  [1]  7.805710  4.273187 32.497013 49.697014  6.400919  6.736593 38.550970
#>  [8] 49.660650 12.503839 19.869379 36.359043 12.356153  5.141713  5.457959
#> [15] 32.076270 22.621857 42.517624  7.303315 19.283361 55.480408 16.385298
#> [22]  4.888481 20.684738 10.297341 10.151362  6.662028 10.159239  7.992649
#> [29]  6.681326  6.308631 22.593502  8.005645  9.072005  6.373744  6.433869
#> [36]  5.599212  5.778799  7.512570 11.863341  6.909616  6.481654  5.300295
#> [43]  6.205842  5.172432  6.773351  6.719396 12.599935  7.117033  3.997243
#> [50]  6.117493  9.481720  8.508172 10.193501  6.657302  4.136527  3.825796
#> [57]  8.147555  9.578996  4.409321 15.075025  5.853627  8.694979  8.776896
#> [64]  5.096160  6.562389  6.529307 12.022185 17.018707  9.602232  5.086577
#> [71] 12.687627  4.044897  4.732786 12.064457 10.283426  4.658615 13.005842
#> [78]  3.757663  4.169741  6.364818  6.746964  4.704956  5.234788
#> 
#> 
#> 
#> $ms1_data
#> Key: <Compound, mz, kmd, rt>
#>      Compound        mz      kmd     rt 50MeOH_10570 50MeOH_10571 50MeOH_10572
#>        <char>     <num>    <num>  <num>        <num>        <num>        <num>
#>   1:     1000 749.51361 0.513606 217.36            0            0            0
#>   2:     1002 734.46721 0.467206 223.54            0            0            0
#>   3:     1003 749.51369 0.513686 223.85            0            0            0
#>   4:     1006  42.03365 0.033646 228.13            0            0            0
#>   5:     1007 137.07102 0.071016 228.43            0            0            0
#>  ---                                                                          
#> 346:      989  83.04884 0.048836 198.47            0            0            0
#> 347:      992 749.51424 0.514236 200.59            0            0            0
#> 348:      996 716.45731 0.457306 206.82            0            0            0
#> 349:      997 749.51417 0.514166 211.44            0            0            0
#> 350:      999 558.36355 0.363546 213.41            0            0            0
#>      FullMix_10574 FullMix_10577 FullMix_10580    cor
#>              <num>         <num>         <num> <lgcl>
#>   1:     226277.39     187860.42     166603.72   TRUE
#>   2:    1285652.63    1244527.38    1222279.13   TRUE
#>   3:     325578.88     357518.06     316568.66   TRUE
#>   4:     270776.94     258514.84     251178.77   TRUE
#>   5:     116058.97     118762.07     102857.34   TRUE
#>  ---                                                 
#> 346:     198559.53     206690.77     217850.67   TRUE
#> 347:      86088.67      83416.42      65714.98   TRUE
#> 348:    1578886.38    1743622.13    1736709.88   TRUE
#> 349:     129910.10     137885.48     128256.25   TRUE
#> 350:     192658.81     218406.34     205511.38   TRUE
#> 
#> $samples
#> [1] "50MeOH_10570"  "50MeOH_10571"  "50MeOH_10572"  "FullMix_10574"
#> [5] "FullMix_10577" "FullMix_10580"
#> 
#> $predicted_molecular_formulas
#> [1] "C4H33NO3S3"
#> 
#> attr(,"class")
#> [1] "mass_data"
```
