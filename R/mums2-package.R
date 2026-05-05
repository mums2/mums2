#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom clustur cluster
#' @importFrom clustur create_sparse_matrix
#' @importFrom clustur get_abundance
#' @importFrom clustur read_count
#' @importFrom clustur read_dist
#' @importFrom data.table .SD
#' @importFrom data.table as.data.table
#' @importFrom data.table data.table
#' @importFrom mpactr filter_cv
#' @importFrom mpactr filter_group
#' @importFrom mpactr filter_insource_ions
#' @importFrom mpactr filter_mispicked_ions
#' @importFrom mpactr get_meta_data
#' @importFrom mpactr get_peak_table
#' @importFrom mpactr import_data
#' @importFrom msentropy msentropy_similarity
#' @importFrom mzR close
#' @importFrom mzR header
#' @importFrom mzR openMSfile
#' @importFrom mzR peaks
#' @importFrom parallel clusterExport
#' @importFrom parallel detectCores
#' @importFrom parallel stopCluster
#' @importFrom pbapply pblapply
#' @importFrom pbapply pboptions
#' @importFrom progress progress_bar
#' @importFrom Rcpp sourceCpp
#' @importFrom Rdisop decomposeIsotopes
#' @importFrom Rdisop decomposeMass
#' @importFrom stats as.dist
#' @importFrom utils head
#' @importFrom utils tail
#' @importFrom xml2 read_xml
#' @importFrom xml2 xml_children
#' @importFrom xml2 xml_find_all
#' @importFrom xml2 xml_name
#' @importFrom xml2 xml_text
#' @useDynLib mums2, .registration = TRUE
## usethis namespace: end
NULL
