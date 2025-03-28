
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mums2

<!-- badges: start -->

<!-- badges: end -->

mums2 is a package that was created to supply a collection of
metabolomic analysis tools. This package will supply tools to quantify
Mass spectrometry (MS) feature abundance, analysis untargeted
metabolites, and find similarities in metabolomic data all while
ensuring it is open source. As of current we support importation and
filteration of mass spectrometry (MS1) feature tables, annotations,
scoring of ms2 spectra, clustering, alpha and beta diversity
calculations, and average subsampled distance calculations using
rarefaction. We are using the mpactR package our team created to import
and filter ms1 data. This allows us to import a number of differentially
formatted peak tables and run numerous filters on the data to ensure for
high-quality data. In addition to mpactR, our team created clustur, a
package that has integrated the mothur clustering algorithms inside of
an R package. Therefore, we can easily incorapte these algorithms inside
of our package to cluster data MS data. In order to anaylsis ms2 data,
we are using the `tidymass` package. This allows us to match ms1
features to ms2 features. We have a created a simple converter for users
to easily transform there data into a tidymass (massdataset) object
aswell (`convert_mpactr_object_to_mass_data_set()`).

Users are able to import and filter their data using the function
wrappers we have created for mpactr. To import data we have created a
`import_all_data()` function. It takes a peak table, metadata table and
a peaktable type (just like in mpactr) to create a mpactr object.
Following this approach we have created a `filter_peak_table()`
function. This function will filter your data based on the paramter
function you supply. And are four paramters functions:
`filter_mispicked_ions_parameters()`, `filter_group_parameters()`,
`filter_cv_parameters()`, and `filter_insource_ions_parameters()`. Each
one of these paramter functions have identical parameters to that of the
mpactr function. Once we have preprocessed all of our data (filtered
added a ms2 spectra, and have converted our mpactr (or mpactr filtered
object) to a massdataset using
`convert_mpactr_object_to_mass_data_set()`), we can annotate ms2 spectra
given a reference database. Using supporting functions for `tidymass` we
were able to created a simplistic function to allow users to annotate
their data. This will give user an idea about the compounds or
metabolites inside of their spectra. To support scoring/distance
calculations or spectra similarity, we have implemeted a `dist_ms2()`
function. This function uses two methods to score spectra, GNPS (Global
Natural Products Social Molecular Networking) or Spectral Entropy. GNPS
uses a cosine scoring algorith to calculate the similarity between two
spectra [Cosine Score
Paper](https://doi.org/10.1016/1044-0305(94)87009-8). In order to ensure
a user could run through our pipeline in a resonable amount of time, we
modified the algorithm to make it much more efficent than that of gnps.
We have confirmed our results to be identical as it still uses the same
algorithm. The other algorithm that we have implemented is Spectral
Entropy. Spectral entropy uses a different algorithm to compute
similarity between spectras (more to be found
[here](https://doi.org/10.1038/s41592-021-01331-z)), but the results are
comparable to that of GNPS’ method.

After the data is scored, we can take advantage of the mothur aglorithms
that are implemented inside of the clustur package. There are five
clustering algorithms in total, furthest, nearest, average, weighted,
and opticlust ([for more on
opticlust](https://doi.org/10.1128/mspheredirect.00073-17)). This allows
us to cluster our data into Operational Metageonomic Units (OMUs) which
are needed to conduct diversity calculations on the data. In order to
cluster the data, we have implemented a wrapper over clusturs
`cluster()` function, and we call it `cluster_data()`. Using this
function you are able to determine how you want to cluster your data and
it returns a shared dataframe (using this shared data frame you have the
option to save it as a csv and use other software like `mothur` for
continued analysis as well). Once this dataframe is created, you have to
convert it into a “community matrix”, before you can actually run
diversity calculations. Therefore, we have created a
`create_commmunity_object()` function. Depending on the size of your
data, the created object may require a substantial amount of memory.

With your community object we can began to calculate diversity metrics.
We have implemented a function called `diversity()`. This function will
conduct a diversity index calculation based on the index you pass to it.
You can conduct alpha diversity (simpson or shannon) or beta diversity
(bray curtis). We also have implemented a function called
`rarefaction()` for rarefaction calculations of the data. This will
allow users to conduct more efficent calculations. Finally, we have
created the `averaged_subsampled_dissimilarity()` function. This
function will run a diversity index calculation and rarefy the data “n”
number of times, then it will calculate the average of all the runs.
This allows users to ensure all of their data was accounted for when
attempting to analyze diversity.

Thus we believe that mums2 will benefit those attempting to analyze
untargeted metabolites inside of their data and other calculations
related to metabolomic analysis. Whether you want to run this pipeline
up to step three or as a whole we belive this package will allow
scientist to easily analyze their data. This package is also open source
and free to use, which means scientist are able to create reproudicible
analysis that others can copy. They also have the ability to create pull
requests or ask for the creation of new features and methods that can be
incorpated inside of mums2. This not only opens the door for a rich
community of scientist in this area of research to join together, but
will allow advances in the community of researchers as a whole.

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
