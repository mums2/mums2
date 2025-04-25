//
// Created by Gregory Johnson on 4/25/25.
//

#include "DiversityMetrics/BetaDiversityCalculators/MorisitahornIndex.h"

double MorisitahornIndex::Calculate(const Rcpp::List& abundanceVectors) const {
    const Rcpp::NumericVector& sampleOne = abundanceVectors[0];
    const Rcpp::NumericVector& sampleTwo = abundanceVectors[1];
    double totalSampleOne = Rcpp::sum(sampleOne);
    double totalSampleTwo = Rcpp::sum(sampleTwo);
    double summationOfSampleOneAndTwo = 0;
    double summationSampleOne = 0;
    double summationSampleTwo = 0;
    for(int i = 0; i < sampleOne.size(); i++) {
        double currentSampleOne = sampleOne[i];
        double currentSampleTwo = sampleTwo[i];
        summationOfSampleOneAndTwo += (currentSampleOne / totalSampleOne)
        * (currentSampleTwo / totalSampleTwo);

        summationSampleOne += std::pow(currentSampleOne/totalSampleOne, 2);
        summationSampleTwo += std::pow(currentSampleTwo/totalSampleTwo, 2);
    }
    return 1 - 2 * (summationOfSampleOneAndTwo / (summationSampleOne + summationSampleTwo));
}
