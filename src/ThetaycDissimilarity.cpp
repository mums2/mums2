//
// Created by Gregory Johnson on 4/25/25.
//

#include "DiversityMetrics/BetaDiversityCalculators/ThetaycDissimilarity.h"

double ThetaycDissimilarity::Calculate(const Rcpp::List& abundanceVectors) const {
    const Rcpp::NumericVector& sampleOne = abundanceVectors[0]; // always the same size
    const Rcpp::NumericVector& sampleTwo = abundanceVectors[1];
    double totalSampleOne = Rcpp::sum(sampleOne);
    double totalSampleTwo = Rcpp::sum(sampleTwo);
    double relativeSummationOfSamples = 0;
    double summationOfDifferenceOfSamples = 0;
    for(int i = 0; i < sampleOne.size(); i++) {
        double relativeSampleOne = sampleOne[i]/totalSampleOne;
        double relativeSampleTwo = sampleTwo[i]/totalSampleTwo;
        relativeSummationOfSamples += relativeSampleOne * relativeSampleTwo;
        summationOfDifferenceOfSamples += std::pow(relativeSampleOne - relativeSampleTwo, 2);
    }
    return 1 - (relativeSummationOfSamples / (summationOfDifferenceOfSamples + relativeSummationOfSamples));
}
