# Read HMDB database

This function allows you to create an hmdb database. However you are
required to supply an xml hmdb file and a folder path that contains all
of the ms2 spectras from the hmdb download page
https://www.hmdb.ca/downloads.

## Usage

``` r
read_hmdb(hmdb_file, ms2_folder)
```

## Arguments

- hmdb_file:

  the xml hmdb file.

- ms2_folder:

  the folder path of your ms2 spectra files.

## Value

a `reference_database` object.

## References

Wishart DS, Guo A, Oler E, Wang F, Anjum A, Peters H, Dizon R, Sayeeda
Z, Tian S, Lee BL, Berjanskii M, Mah R, Yamamoto M, Jovel J,
Torres-Calzada C, Hiebert-Giesbrecht M, Lui VW, Varshavi D, Varshavi D,
Allen D, Arndt D, Khetarpal N, Sivakumaran A, Harford K, Sanford S, Yee
K, Cao X, Budinski Z, Liigand J, Zhang L, Zheng J, Mandal R, Karu N,
Dambrova M, Schiöth HB, Greiner R, Gautam V. HMDB 5.0: the Human
Metabolome Database for 2022. Nucleic Acids Res. 2022 Jan
7;50(D1):D622-D631. doi: 10.1093/nar/gkab1062. PMID: 34986597; PMCID:
PMC8728138.

## Examples

``` r
read_msp(mums2_example("massbank_example_data.msp" ))
#> Reading: /home/runner/work/_temp/Library/mums2/extdata/massbank_example_data.msp ...
#> Computing                                | 0%  ETA: -...Computing ■                              | 3%  ETA: ...Computing ■■                             | 6%  ETA: ...Computing ■■■                            | 10%  ETA: ...Computing ■■■■                           | 13%  ETA: ...Computing ■■■■■                          | 16%  ETA: ...Computing ■■■■■■                         | 20%  ETA: 3s ...Computing ■■■■■■■                        | 23%  ETA: 3s ...Computing ■■■■■■■■                       | 26%  ETA: 2s ...Computing ■■■■■■■■■                      | 30%  ETA: 2s ...Computing ■■■■■■■■■■                     | 33%  ETA: 1s ...Computing ■■■■■■■■■■■                    | 36%  ETA: 1s ...Computing ■■■■■■■■■■■■                   | 40%  ETA: 1s ...Computing ■■■■■■■■■■■■■                  | 43%  ETA: 1s ...Computing ■■■■■■■■■■■■■■                 | 46%  ETA: 1s ...Computing ■■■■■■■■■■■■■■■                | 50%  ETA: ...Computing ■■■■■■■■■■■■■■■■               | 53%  ETA: ...Computing ■■■■■■■■■■■■■■■■■              | 56%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■             | 60%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■            | 63%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■           | 66%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■          | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■         | 73%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■        | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■       | 80%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■      | 83%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■     | 86%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 90%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 93%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  | 96%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
#> You have 8068 references in this object.
```
