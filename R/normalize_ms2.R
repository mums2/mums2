#' Normalize ms2 peak intensities
#'
#' @description
#' `normalize_ms2` normalizes peak intensities within a single ms2
#' scan of tandem mass spectromentry.
#'
#' @details
#' `normalize_ms2` takes a vector of intensity values from a single
#' ms2 scan and returns normalized values.
#'
#' normalization methods are controlled within the `method` argument,
#' which currently supports square root (`method = "sqrt"`) and scale
#' (`method = "scale"`) normalization.
#'
#' Method `"sqrt"` normalizes square root peak intensities relative to
#' the square root of summed intensities in the scan using the following
#' equation:
#'
#' \deqn{sqrt(intensity_i) / swrt(sum(intensity))}.
#'
#' Square root normalization is used to calculate the gnps-like cosine score.
#'
#' Method `"scale"` normalzies peak intensities relative to the highest peak
#' in the scan using the following equation:
#'
#' \deqn{intensity_i / max(intensity)}
#'
#'
#' @param data A `numeric` vector of intensity values.
#' @param method A method for normalization. May be one of `c("sqrt", "scale")`.
#'
#' @return A data structure of the same class as input with normalized
#' intensities in the `intensity` column.
#'
#' @examples
#' intensity <- c(20, 2, 103, 6)
#'
#' dat_sqrt <- normalize_ms2(intensity, method = "sqrt")
#' dat_sqrt
#'
#' dat_scale <- normalize_ms2(intensity, method = "scale")
#' dat_scale
#'
#' @noRd
normalize_ms2 <- function(data, method) {
  if (method == "sqrt") {
    return(squareRootNormalize(vec = data))
  }

  if (method == "scale") {
    return(scaleNormalize(vec = data))
  }
}
