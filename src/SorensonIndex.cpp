//
// Created by Gregory Johnson on 4/25/25.
//

#include "DiversityMetrics/BetaDiversityCalculators/SorensonIndex.h"

double SorensonIndex::Calculate(const Rcpp::List& abundanceVectors) const {
    const Rcpp::NumericVector& sampleOne = abundanceVectors[0];
    const Rcpp::NumericVector& sampleTwo = abundanceVectors[1];
    double sharedOmus = 0;
    double sampleOneOmus = 0;
    double sampleTwoOmus = 0;
    for(int i = 0; i < sampleOne.size(); i++) {
        double firstCurrentOmu = sampleOne[i];
        double secondCurrentOmu = sampleTwo[i];
        if(firstCurrentOmu > 0 && secondCurrentOmu > 0) {
            sharedOmus++;
            sampleOneOmus++;
            sampleTwoOmus++;
            continue;
        }
        if(firstCurrentOmu > 0) sampleOneOmus++;
        if(secondCurrentOmu > 0) sampleTwoOmus++;
    }
    return 1 - 2 * sharedOmus / (sampleOneOmus + sampleTwoOmus);
}
