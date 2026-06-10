# Match your ms1 spectra to a ms2

We are matching your ms1 to your supplied ms2 by looking at the
difference between the mz and rt.

## Usage

``` r
ms2_ms1_compare(ms2_files, mpactr_object, mz_tolerance, rt_tolerance)
```

## Arguments

- ms2_files:

  a list of all your mgf, mzml, or mzxml files.

- mpactr_object:

  your mpactr object created from
  [`import_all_data()`](https://www.mums2.org/mums2/reference/import_all_data.md)

- mz_tolerance:

  your mass-charge ratio tolerance in ppm (parts per million).

- rt_tolerance:

  your retention time tolerance.

## Value

returns a `mass_data` object of all of the ms2 and ms1 matches.

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
#> Reading: /home/runner/work/_temp/Library/mums2/extdata/botryllus_v2.gnps.mgf ...
#> Computing                                | 0%  ETA: -...Computing ■                              | 3%  ETA: ...Computing ■■                             | 6%  ETA: ...Computing ■■■                            | 10%  ETA: ...Computing ■■■■                           | 13%  ETA: ...Computing ■■■■■                          | 16%  ETA: ...Computing ■■■■■■                         | 20%  ETA: ...Computing ■■■■■■■                        | 23%  ETA: ...Computing ■■■■■■■■                       | 26%  ETA: ...Computing ■■■■■■■■■                      | 30%  ETA: ...Computing ■■■■■■■■■■                     | 33%  ETA: ...Computing ■■■■■■■■■■■                    | 36%  ETA: ...Computing ■■■■■■■■■■■■                   | 40%  ETA: ...Computing ■■■■■■■■■■■■■                  | 43%  ETA: ...Computing ■■■■■■■■■■■■■■                 | 46%  ETA: ...Computing ■■■■■■■■■■■■■■■                | 50%  ETA: ...Computing ■■■■■■■■■■■■■■■■               | 53%  ETA: ...Computing ■■■■■■■■■■■■■■■■■              | 56%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■             | 60%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■            | 63%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■           | 66%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■          | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■         | 73%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■        | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■       | 80%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■      | 83%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■     | 86%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 90%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 93%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  | 96%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
#> 17/349 peaks have an MS2 spectra.
```
