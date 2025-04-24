//
// Created by Gregory Johnson on 4/24/25.
//

#include "DiversityMetrics/BetaDiversityCalculators/EuclideanDistance.h"

double EuclideanDistance::Calculate(const Rcpp::List& abundanceVectors) const {
    const Rcpp::NumericVector& sampleOne = abundanceVectors[0];
    const Rcpp::NumericVector& sampleTwo = abundanceVectors[1];
    double distance = 0;
    for(int i = 0; i < sampleOne.size(); i++) {
        distance += pow(sampleOne[i] - sampleTwo[i], 2);
    }
    return std::sqrt(distance);
}
