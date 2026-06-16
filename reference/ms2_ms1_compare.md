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
#> 17/349 peaks have an MS2 spectra.
```
