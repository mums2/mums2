
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
  if(!"reference_database" %in% class(reference)) {
    stop("Ensure reference is the object generated from `read_msp()` or `read_hmdb()`")
  }
  if(class(index) != "numeric") {
    stop("index has to be a numeric")
  }
  return(GetNode(reference, index))
}


#' @title Combine Databases
#' @export
#' @description Add another database to your reference database.
#' @param reference reference database object.
#' @param other_reference your other reference database object
#' @examples
#' reference <- read_msp(mums2_example("PSU-MSMLS.msp"))
#' reference2 <- read_msp(mums2_example("PSU-MSMLS.msp"))
#' add_references(reference, reference2)
#'
#' @return prints customized message to the console
combined_reference_database <- function(reference, other_reference) {
  if(!"reference_database" %in% class(reference) || 
     !"reference_database" %in% class(other_reference)) {
      stop("Ensure reference is the object generated from `read_msp()` or `read_hmdb()`")
  }
  new_reference_db <- CombineReferenceDatabases(reference, other_reference)
  class(new_reference_db) <- "reference_database"
  return(new_reference_db)
}


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

