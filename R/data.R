#' Get file paths for examples
#'
#' mums2 contains a number of example files in the `inst/extdata` directory.
#' This function makes them accessible in documentation that shows how file
#' paths are used in function examples.
#'
#' @param file Name of a file. If `NULL`, all examples files will be listed.
#'
#' @export
#' @return A file path to example data stored in the `inst/extdata` directory
#' of the package.
#' @examples
#' mums2_example()
#'
#' mums2_example("massbank_example_data.msp")
#' @return returns a `character` object
mums2_example <- function(file = NULL) {
  path <- ""
  if (is.null(file)) {
    path <- dir(system.file("extdata", package = "mums2"))
  } else {
    path <- system.file("extdata", file, package = "mums2", mustWork = TRUE)
  }
  return(path)

}
