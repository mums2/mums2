
#' @title Get Reference Data
#' @export
#' @description Will return the data inside the reference object
#' based on the index given.
#' @param reference reference database object.
#' @param index the index of the data. The index starts at 1.
#' @examples
#' reference <- read_msp(mums2_example("massbank_example_data.msp"))
#' get_reference_data(reference, 1)
#'
#' @return returns a `list` object with all of the reference data at
#' the specified index.
get_reference_data <- function(reference, index) {
  if (!inherits(reference, "reference_database")) {
    stop(paste0("Ensure reference is the object generated from",
                " `read_msp()` or `read_hmdb()`"))
  }
  if (!inherits(index, "numeric")) {
    stop("Index has to be a numeric")
  }
  if (index <= 0 || index > length(reference)) {
    stop("Index must be less than the size of the database")
  }
  index <- index - 1
  GetNode(reference, index)
}


#' @title Combine Databases
#' @export
#' @description Add another database to your reference database.
#' @param reference reference database object.
#' @param other_reference your other reference database object
#' @examples
#' reference <- read_msp(mums2_example("massbank_example_data.msp"))
#' reference2 <- read_msp(mums2_example("massbank_example_data.msp"))
#' combined_reference_database(reference, reference2)
#'
#' @return a `reference_database` that includes references from both
#' reference databases.
combined_reference_database <- function(reference, other_reference) {
  if (!inherits(reference, "reference_database") ||
        !inherits(other_reference, "reference_database")) {
    stop(paste0("Ensure reference is the object generated from",
                "`read_msp()` or `read_hmdb()`"))
  }
  new_reference_db <- CombineReferenceDatabases(reference, other_reference)
  class(new_reference_db) <- "reference_database"
  new_reference_db
}


#' @title Print reference
#' @export
#' @description print reference objects.
#' @param x reference database object.
#' @param ... any extra print arguments you want to include.
#' @examples
#' reference <- read_msp(mums2_example("massbank_example_data.msp"))
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
#' reference <- read_msp(mums2_example("massbank_example_data.msp"))
#' length(reference)
#'
#' @return returns the length of the regerence database
length.reference_database <- function(x) {
  GetNodeCount(x)
}
