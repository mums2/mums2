#' Annotate LC-MS/MS features
#' @description
#' Annotate query LC-MS/MS features in a `mass_dataset` object given a
#' reference list.
#'
#' `annotate_ms2()` allows for annotation of mass spectrometry features.
#' Similaritity between query and refrence level 2 spectra are determined via
#' spectral scoring methods. Currently scoring methods `"gnps"` and
#' `"spectral_entropy"` are supported. The scoring method is specified by the
#' `score_params` argument. `score_params` is a list of parameters for the
#' chosen scoring method. Parameters for "gnps" and "spectral_entropy" can be
#' created with notur functions [gnps_params()] and [spec_entropy_params()],
#' respectively.
#'
#'
#' @param query A `mass_dataset` object containing ms2 data.
#' @param reference A list of reference data downloaded from
#'  \href{https://massdatabase.tidymass.org}{massdatabase}.
#' for more information about this class.
#' @param score_params Parameters for scoring method to be applied.
#' @param precursor_tolerance Precursor mz tolerance. MS2 scans with a
#'  difference in precursor mz less than or equal to this value will be scored.
#' @param min_score Similarity score threshold to determine a match for
#'  annotation. Comparisons with scores below this value will not be reported.
#'
#' @return A `data.frame` with all comparisons with scores above the threshold.
#'  Information for the query scan include `query_ms1_id` (the variable_id
#'  for features in expression_data of the `mass_data` object)
#'  `"query_ms2_id"` (the `ms2_spectrum_id` in the `query`
#'  object), `"query_mz"` (the precursor mz for the scan), and `"query_rt"`
#'  (the retention time for the scan). `query_mz` and `query_rt` are derived
#'  from the ms2 matches data. A column (`"ref_idx`) is included to report the
#'  location for the matching reference molecule in `"reference"`. Scores
#'  are reported in the `"score"` column. Annotation information is returned
#'  given the information provided in the reference used as input.
#'
#' @importFrom stats setNames
#'
#' @usage annotate_ms2(query, reference, score_params,
#' precursor_tolerance, min_score)
#'
#' @export
#' @rdname annotate_ms2
annotate_ms2 <- function(query, reference, score_params,
                         precursor_tolerance, min_score) {
  UseMethod("annotate_ms2", query)
}

#' @method annotate_ms2 mass_data
#' @export
annotate_ms2.mass_data <- function(query, reference, score_params,
  precursor_tolerance, min_score) {
  ms2 <- query$ms2_matches
  matches <- AnnotateMs2Features(ms2$ms1_compound_id, ms2$ms2_spectrum_id,
                                 ms2$mz, ms2$rt, 
                                 query$peak_data,
                                 reference, score_params, precursor_tolerance,
                                 min_score)

annotations <- add_annotations(matches, reference)

return(annotations)
}

#' @importFrom stats setNames
add_annotations <- function(matches, reference) {
  ref_info <- data.frame()

  for (ref_idx in matches$ref_idx) {
    new_ann <- stats::setNames(data.frame(t(reference[[ref_idx]]$info[, -1])),
                               reference[[ref_idx]]$info[, 1])

    ref_info <- rbind(ref_info, new_ann)
  }

  matches <- cbind(matches, ref_info)

  return(matches)
}
