# Print reference

print reference objects.

## Usage

``` r
# S3 method for class 'reference_database'
print(x, ...)
```

## Arguments

- x:

  reference database object.

- ...:

  any extra print arguments you want to include.

## Value

prints customized message to the console

## Examples

``` r
reference <- read_msp(mums2_example("massbank_example_data.msp"))
#> Reading: /home/runner/work/_temp/Library/mums2/extdata/massbank_example_data.msp ...
#> Computing                                | 0%  ETA: -...Computing ■                              | 3%  ETA: ...Computing ■■                             | 6%  ETA: ...Computing ■■■                            | 10%  ETA: ...Computing ■■■■                           | 13%  ETA: ...Computing ■■■■■                          | 16%  ETA: ...Computing ■■■■■■                         | 20%  ETA: ...Computing ■■■■■■■                        | 23%  ETA: ...Computing ■■■■■■■■                       | 26%  ETA: ...Computing ■■■■■■■■■                      | 30%  ETA: ...Computing ■■■■■■■■■■                     | 33%  ETA: ...Computing ■■■■■■■■■■■                    | 36%  ETA: 1s ...Computing ■■■■■■■■■■■■                   | 40%  ETA: 1s ...Computing ■■■■■■■■■■■■■                  | 43%  ETA: 1s ...Computing ■■■■■■■■■■■■■■                 | 46%  ETA: 1s ...Computing ■■■■■■■■■■■■■■■                | 50%  ETA: ...Computing ■■■■■■■■■■■■■■■■               | 53%  ETA: ...Computing ■■■■■■■■■■■■■■■■■              | 56%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■             | 60%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■            | 63%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■           | 66%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■          | 70%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■         | 73%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■        | 76%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■       | 80%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■      | 83%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■     | 86%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■    | 90%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■   | 93%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  | 96%  ETA: ...Computing ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■ | 100%  ETA: ...
print(reference)
#> You have 8068 references in this object.
```
