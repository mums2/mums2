# GNPS-like similarity between two MS/MS spectra

`modified_cosine_params()` generates a parameter list to perform
GNPS-like cosine similarity score calculation between two MS2 spectra.

## Usage

``` r
modified_cosine_params(frag_tolerance)
```

## Arguments

- frag_tolerance:

  The mz fragment tolerance threshold for aligning fragment peaks from
  two ms2 spectra. GNPS default = 0.5.

## Value

A parameters list for similarity scoring method "gnps"

## Details

`modified_cosine_params()` will initiate cosine scoring based on the
Python code by Wang et al. (2016), which is currently used for cosine
scoring in GNPS, to calculate similarity between two MS2 spectra. This
scoring method will compare peaks data, apply a square root
normalization to peak intensities, align peaks both with and without
correction for mass shifts, and calculate similarity.

## References

Mingxun Wang, Jeremy J. Carver, Vanessa V. Phelan, Laura M. Sanchez,
Neha Garg, Yao Peng, Don Duy Nguyen et al. "Sharing and community
curation of mass spectrometry data with Global Natural Products Social
Molecular Networking." Nature biotechnology 34, no. 8 (2016): 828. PMID:
27504778

## Examples

``` r
modified_cosine_params(0.5)
#> $tolerance
#> [1] 0.5
#> 
#> $method
#> [1] "gnps"
#> 
```
