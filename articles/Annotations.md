# Annotations

``` r
library(mums2)
```

To allow users to access data inside of reference databases without
causing lag, we opted to created reference databases as a Rcpp object.
This means that users are not able to access them unless you use
specific getter functions
([`get_reference_data()`](https://www.mums2.org/mums2/reference/get_reference_data.md),
and
[`combined_reference_database()`](https://www.mums2.org/mums2/reference/combined_reference_database.md)).

To create a reference database object, you can either use the
[`read_msp()`](https://www.mums2.org/mums2/reference/read_msp.md)
function or the
[`read_hmdb()`](https://www.mums2.org/mums2/reference/read_hmdb.md)
function. Using the [`length()`](https://rdrr.io/r/base/length.html)
function will allow you to see the amount of references.

``` r
reference_db <- read_msp(msp_file = mums2_example("massbank_example_data.msp"))
#> [1] "Reading: /home/runner/work/_temp/Library/mums2/extdata/massbank_example_data.msp ..."
#> Computing                                                    | 0%  ETA: -...Computing ■                                                  | 2%  ETA: ...Computing ■■                                                 | 4%  ETA: ...Computing ■■■                                                | 6%  ETA: ...Computing ■■■■                                               | 8%  ETA: ...Computing ■■■■■                                              | 10%  ETA: ...Computing ■■■■■■                                             | 12%  ETA: ...Computing ■■■■■■■                                            | 14%  ETA: ...Computing ■■■■■■■■                                           | 16%  ETA: ...Computing ■■■■■■■■■                                          | 18%  ETA: ...Computing ■■■■■■■■■■                                         | 20%  ETA: ...Computing ■■■■■■■■■■■                                        | 22%  ETA: ...Computing ■■■■■■■■■■■■                                       | 24%  ETA: ...Computing ■■■■■■■■■■■■■                                      | 26%  ETA: ...Computing ■■■■■■■■■■■■■■                                     | 28%  ETA: ...Computing ■■■■■■■■■■■■■■■                                    | 30%  ETA: ...Computing ■■■■■■■■■■■■■■■■                                   | 32%  ETA: ...Computing ■■■■■■■■■■■■■■■■■                                  | 34%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■                                 | 36%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■                                | 38%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■                               | 40%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■                              | 42%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■                             | 44%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■                            | 46%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■                           | 48%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■                          | 50%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■                         | 52%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■                        | 54%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■                       | 56%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                      | 58%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                     | 60%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                    | 62%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                   | 64%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                  | 66%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                 | 68%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■               | 72%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■              | 74%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■             | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■            | 78%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■           | 80%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■          | 82%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■         | 84%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■        | 86%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■       | 88%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■      | 90%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     | 92%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 94%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 96%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
print(reference_db)
#> [1] "You have 8068 references in this object."
length(reference_db)
#> [1] 8068
```

The
[`get_reference_data()`](https://www.mums2.org/mums2/reference/get_reference_data.md)
function will return the reference at the specified index you input.
This allows you to see the data that is inside of your database.

``` r
get_reference_data(reference_db, 1)
#> $info
#> $info$keys
#>  [1] "name"            "precursormz"     "precursortype"   "formula"        
#>  [5] "ontology"        "inchikey"        "smiles"          "retentiontime"  
#>  [9] "ccs"             "ionmode"         "instrumenttype"  "instrument"     
#> [13] "collisionenergy" "comment"         "num peaks"      
#> 
#> $info$values
#>  [1] "Veratramine; LC-ESI-TOF; MS2; CE\r"                                                           
#>  [2] "410.305356\r"                                                                                 
#>  [3] "[M+H]+\r"                                                                                     
#>  [4] "C27H39NO2\r"                                                                                  
#>  [5] "Fluorenes\r"                                                                                  
#>  [6] "MALFODICFSIXPO-KFKQDBFTSA-N\r"                                                                
#>  [7] "C[C@H]([C@@H]1NC[C@@H](C)C[C@H]1O)C1=C(C)C2=C(C=C1)[C@@H]1CC=C3C[C@@H](O)CC[C@]3(C)[C@H]1C2\r"
#>  [8] "\r"                                                                                           
#>  [9] "207.3842046\r"                                                                                
#> [10] "Positive\r"                                                                                   
#> [11] "LC-ESI-QTOF\r"                                                                                
#> [12] "\r"                                                                                           
#> [13] "\r"                                                                                           
#> [14] "registered in MassBank\r"                                                                     
#> [15] "83\r"                                                                                         
#> 
#> 
#> $spec
#> $spec$mz
#>  [1]  84.1 105.1 107.1 114.1 115.1 119.1 121.1 124.1 125.1 129.1 131.1 132.1
#> [13] 133.1 134.1 144.1 145.1 146.1 147.1 151.1 155.1 156.1 157.1 158.1 159.1
#> [25] 160.1 161.1 167.1 168.1 169.1 170.1 171.1 172.1 173.1 175.1 181.1 182.1
#> [37] 183.1 184.1 185.1 193.1 194.1 195.1 196.1 197.1 198.1 206.1 207.1 208.1
#> [49] 209.1 211.1 212.2 219.1 220.1 221.1 222.1 223.1 233.1 234.1 235.2 236.2
#> [61] 237.2 247.2 248.2 249.2 251.2 262.2 263.2 277.2 278.2 280.2 281.2 295.2
#> [73] 296.2 297.2 309.2 319.2 320.2 333.2 375.3 392.3 393.3 396.3 410.3
#> 
#> $spec$intensity
#>  [1]   32    6    6  105    7    8    6   83    8    7  129   13  212   23    5
#> [16]  108   12   14    8   26   17  372   45  216   29    6    6    7  121   17
#> [31]  350   47   10    5   27   12  108   18   18   16    6   37   12   41    7
#> [46]    8   28    8   23   51    9    6   11   35   13   11   16    8   34   13
#> [61]   12   16   18   17    6   42   11   77   17   17   22 1000  225   26    7
#> [76]   10    6    9    6   62   18    8   16
get_reference_data(reference_db, 2)
#> $info
#> $info$keys
#>  [1] "name"            "precursormz"     "precursortype"   "formula"        
#>  [5] "ontology"        "inchikey"        "smiles"          "retentiontime"  
#>  [9] "ccs"             "ionmode"         "instrumenttype"  "instrument"     
#> [13] "collisionenergy" "comment"         "num peaks"      
#> 
#> $info$values
#>  [1] "Cyclopamine; LC-ESI-TOF; MS2; CE\r"                                                             
#>  [2] "412.321006\r"                                                                                   
#>  [3] "[M+H]+\r"                                                                                       
#>  [4] "C27H41NO2\r"                                                                                    
#>  [5] "Jerveratrum-type alkaloids\r"                                                                   
#>  [6] "QASFUMOKHFSJGL-LAFRSMQTSA-N\r"                                                                  
#>  [7] "C[C@@H]1[C@@H]2NC[C@@H](C)C[C@H]2O[C@]11CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC2=C1C\r"
#>  [8] "\r"                                                                                             
#>  [9] "207.2302643\r"                                                                                  
#> [10] "Positive\r"                                                                                     
#> [11] "LC-ESI-QTOF\r"                                                                                  
#> [12] "\r"                                                                                             
#> [13] "\r"                                                                                             
#> [14] "registered in MassBank\r"                                                                       
#> [15] "99\r"                                                                                           
#> 
#> 
#> $spec
#> $spec$mz
#>  [1]  67.1  81.1  84.1  85.1  93.1  96.1 102.1 105.1 107.1 109.1 110.1 112.1
#> [13] 114.1 115.1 119.1 121.1 124.1 125.1 126.1 127.1 129.1 131.1 133.1 134.1
#> [25] 135.1 140.1 143.1 144.1 145.1 146.1 147.1 148.1 149.1 150.1 151.1 155.1
#> [37] 156.1 157.1 158.1 159.1 160.1 161.1 162.1 163.1 167.1 169.1 170.1 171.1
#> [49] 172.1 173.1 174.1 175.1 177.1 181.1 182.1 183.1 184.1 185.1 186.1 187.1
#> [61] 195.1 196.1 197.1 198.1 199.1 200.2 201.2 207.1 209.1 210.1 211.1 212.2
#> [73] 213.2 214.2 221.1 222.1 223.1 224.2 225.2 227.2 235.1 237.2 239.2 249.2
#> [85] 251.2 253.2 264.2 269.2 279.2 280.2 281.2 287.2 295.2 297.2 321.2 322.2
#> [97] 377.3 394.3 412.3
#> 
#> $spec$intensity
#>  [1]   29   47   87   70   30   30   68   52   82  724   62   53  693   41  110
#> [16]   94   50   27  295   25   32  151  607   60  149  185  277   45  889   99
#> [31]  271   31   48   46   76  107   55 1000  132  860  104  334   37   80   24
#> [46]  162   37  566   82  330   43   57   26   58   38  218   60  244   37   76
#> [61]  140   56  254   47  138   28   26   29  113   35  174   28  352   53   31
#> [76]   28   79   25   61   27   43   71   62   30   58   79   27   28  100   28
#> [91]   35   42   58   43  130   38   44   55  113
```

The
[`combined_reference_database()`](https://www.mums2.org/mums2/reference/combined_reference_database.md)
function allows you to combined two different databases ontop of
eachother to make one reference database.

``` r
reference_db_2 <- read_msp(msp_file = mums2_example("massbank_example_data_negative.msp"))
#> [1] "Reading: /home/runner/work/_temp/Library/mums2/extdata/massbank_example_data_negative.msp ..."
#> Computing                                                    | 0%  ETA: -...Computing ■                                                  | 2%  ETA: ...Computing ■■                                                 | 4%  ETA: ...Computing ■■■                                                | 6%  ETA: ...Computing ■■■■                                               | 8%  ETA: ...Computing ■■■■■                                              | 10%  ETA: ...Computing ■■■■■■                                             | 12%  ETA: ...Computing ■■■■■■■                                            | 14%  ETA: ...Computing ■■■■■■■■                                           | 16%  ETA: ...Computing ■■■■■■■■■                                          | 18%  ETA: ...Computing ■■■■■■■■■■                                         | 20%  ETA: ...Computing ■■■■■■■■■■■                                        | 22%  ETA: ...Computing ■■■■■■■■■■■■                                       | 24%  ETA: ...Computing ■■■■■■■■■■■■■                                      | 26%  ETA: ...Computing ■■■■■■■■■■■■■■                                     | 28%  ETA: ...Computing ■■■■■■■■■■■■■■■                                    | 30%  ETA: ...Computing ■■■■■■■■■■■■■■■■                                   | 32%  ETA: ...Computing ■■■■■■■■■■■■■■■■■                                  | 34%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■                                 | 36%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■                                | 38%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■                               | 40%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■                              | 42%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■                             | 44%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■                            | 46%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■                           | 48%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■                          | 50%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■                         | 52%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■                        | 54%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■                       | 56%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                      | 58%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                     | 60%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                    | 62%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                   | 64%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                  | 66%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                 | 68%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■                | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■               | 72%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■              | 74%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■             | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■            | 78%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■           | 80%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■          | 82%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■         | 84%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■        | 86%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■       | 88%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■      | 90%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     | 92%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 94%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 96%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
length(reference_db_2)
#> [1] 4782

combined_db <- combined_reference_database(reference_db, reference_db_2)
length(combined_db)
#> [1] 12850
```
