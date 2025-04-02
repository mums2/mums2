#include "NormalizeMs2/SquareRootNormalize.h"
#include "NormalizeMs2/ScaleNormalize.h"
#include <Rcpp.h>

// [[Rcpp::export]]
std::vector<double> squareRootNormalize(std::vector<double>& vec) {
    SquareRootNormalize normalize;
    vec = normalize.Normalize(vec);

    return(vec);
}

// [[Rcpp::export]]
std::vector<double> scaleNormalize(std::vector<double>& vec) {
    ScaleNormalize normalize;
    vec = normalize.Normalize(vec);

    return(vec);
}

