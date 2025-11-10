# Change RT time to minutes or seconds

This function will change your ms1 peak table rt time to rt time in
seconds or minutes. This modification happens in place (or by
reference), so the `mpactr_object` will be updated.

## Usage

``` r
change_rt_to_seconds_or_minute(mpactr_object, rt_type = "seconds")
```

## Arguments

- mpactr_object:

  The object created from
  [`import_all_data()`](https://www.mums2.org/mums2/reference/import_all_data.md)
  file.

- rt_type:

  how you want to convert your retention time, your options are minutes,
  or seconds. defaults to seconds.

## Value

a modified `mpactr` object.

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
change_rt_to_seconds_or_minute(data, "minutes")
#> [1] "Changing rt values to minutes"
#>       Compound        mz RTINMINUTES 50MeOH_10570 50MeOH_10571 50MeOH_10572
#>         <char>     <num>       <num>        <num>        <num>        <num>
#>    1:        1 111.11648  0.04683333      4166.40     10104.69      7403.40
#>    2:        2 140.06824  0.04816667      3674.63      5520.20      3303.08
#>    3:        3 173.07892  0.05783333      1059.17       150.01       647.51
#>    4:        4 158.00288  0.05983333     11022.72      7634.81     10692.11
#>    5:        5  57.06977  0.06116667     10547.08     19923.64     15259.67
#>   ---                                                                      
#> 1290:     1209 287.88849  6.76383333     57581.37     56411.48     51881.45
#> 1291:     1210 246.86210  6.75750000     26356.65     26702.41     21675.92
#> 1292:     1259 241.88354  7.02666667     46642.71     44577.70     23018.57
#> 1293:     1275 236.93797  7.14633333     70893.87    100240.23     33815.00
#> 1294:     1287 163.94000  7.30300000     10838.34      5769.84     25185.58
#>       FullMix_10574 FullMix_10577 FullMix_10580      kmd
#>               <num>         <num>         <num>    <num>
#>    1:       8752.06      12820.89       1558.26 0.116476
#>    2:      15006.73       1182.62       8596.09 0.068236
#>    3:       8780.84       2584.23      35997.61 0.078916
#>    4:       4627.74      12896.38       2925.34 0.002876
#>    5:       2398.49      11353.59       1914.96 0.069766
#>   ---                                                   
#> 1290:      46938.96      41305.46      31569.37 0.888490
#> 1291:      22908.83      16813.92      17843.02 0.862100
#> 1292:      40485.21      41501.13      40877.48 0.883540
#> 1293:      31953.40     101686.62       9737.37 0.937970
#> 1294:     118191.05       9190.82     145368.44 0.940000
```
