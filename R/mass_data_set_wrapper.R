#' @export
#' @title Add MS2 Spectras
#' @description
#' A wrapper for the `mutate_ms2()` function from mass_data_set. It will match your ms1 features inside of your massdataset object and add ms2 features to it.
#' @param mass_data_set your massdataset object.
#' @param path the path to your ms2 files.
#' @param polarity the polarity of your mgf file, it can be "negative" or "positive".
#' @param column the column type of your data, can choose between "rp" or "hilic".
#' @param ms1.ms2.match.mz.tol your mz (mass to charge ratio) tolerance for matching ms1 and ms2.
#' @param ms1.ms2.match.rt.tol your retention time tolerance for ms1 and s2 matching.
add_ms2 <- function(mass_data_set, path, polarity, column, ms1.ms2.match.mz.tol = 15, ms1.ms2.match.rt.tol = 30) {
  browser()
  return(mutate_ms2(
    object = mass_data_set, 
    path = path,
    polarity = polarity,
    column = column,
    ms1.ms2.match.mz.tol = ms1.ms2.match.mz.tol,
    ms1.ms2.match.rt.tol = ms1.ms2.match.rt.tol
  ))
}