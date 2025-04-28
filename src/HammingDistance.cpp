//
// Created by Gregory Johnson on 4/24/25.
//

#include "DiversityMetrics/BetaDiversityCalculators/HammingDistance.h"

double HammingDistance::Calculate(const Rcpp::List& abundanceVectors) const {
    const Rcpp::NumericVector& sampleOne = abundanceVectors[0];
    const Rcpp::NumericVector& sampleTwo = abundanceVectors[1];
    return Rcpp::sum(Rcpp::abs(sampleOne - sampleTwo));
}
