#' @title Print reference
#' @export
#' @description print reference objects.
#' @param x reference database object.
#' @examples
#' reference <- read_msp(mums2_example("PSU-MSMLS.msp"))[[1]]
#' print(reference)
#'
#' @return prints customized message to the console
print.reference_database <- function(x, ...) {
  node_count <- GetNodeCount(x)
  print(paste0("You have ", node_count, " references in this object."))
}


#' @title Get Reference Data
#' @export
#' @description Will return the data inside the reference object based on the index given.
#' @param x reference database object.
#' @param index the index of the data.
#' @examples
#' reference <- read_msp(mums2_example("PSU-MSMLS.msp"))
#' get_reference_data(reference, 1)
#'
#' @return prints customized message to the console
get_reference_data <- function(reference, index) {
  return(GetNode(reference, index))
}


#' @title Reference database length
#' @export
#' @description returns the length of the database
#' @param x reference database object.
#' @examples
#' reference <- read_msp(mums2_example("PSU-MSMLS.msp"))
#' length(references)
#'
#' @return returns the length of the regerence database
length.reference_database <- function(reference) {
  return(GetNodeCount(reference))
}

