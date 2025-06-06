
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mums2

<!-- badges: start -->

<!-- badges: end -->

mums2 is a package that was created to supply a collection of
metabolomic analysis tools. This package will supply tools to analysis
Mass spectrometry (MS) data in a easy and efficent manner. As of current
we support importation and filteration of mass spectrometry (MS1)
feature tables, annotations, de-novo molecular formula predictions,
scoring of ms2 spectra, clustering, alpha and beta diversity
calculations, and average subsampled distance calculations using
rarefaction. We are using the mpactR package our team created to import
and filter ms1 data. This allows us to import a number of differentially
formatted peak tables and run numerous filters on the data to ensure for
high-quality data. In addition to mpactR, our team created clustur, a
package that has integrated the mothur clustering algorithms inside of
an R package. With this collection of methods, we belive we can create a
powerful package that can rival other non open-source alternatives.

## Installation

You can install the CRAN version with:

``` r
install.packages("mums2")
```

You can install the development version of mpactr from
[GitHub](https://github.com/mums2/mums2) with:

``` r
devtools::install_github("mums2/mums2")
```

## Get started

See the [Getting
Started](https://www.mums2.org/mpactr/articles/mpactr.html) page to get
started.

## Getting help

If you encounter an issue, please file an issue on
[GitHub](https://github.com/mums2/mums2/issues). Please include a
minimal reproducible example with your issue.

## Contributing

Is there a feature you’d like to see included, please let us know! Pull
requests are welcome on [GitHub](https://github.com/mums2/mums2/pulls).
