//
// Created by Gregory Johnson on 4/24/25.
//

#include "DiversityMetrics/BetaDiversityCalculators/JaccardDistance.h"

double JaccardDistance::Calculate(const Rcpp::List& abundanceVectors) const {
    const Rcpp::NumericVector& sampleOne = abundanceVectors[0];
    const Rcpp::NumericVector& sampleTwo = abundanceVectors[1];
    double sumOfAllOmus = Rcpp::sum(sampleOne) + Rcpp::sum(sampleTwo);
    double amountOfSpeciesInBothSamples = 0;
    for(int i = 0; i < abundanceVectors.size(); i++) {
        if(sampleOne[i] > 0 && sampleTwo[i] > 0) {
            amountOfSpeciesInBothSamples += sampleOne[i] + sampleTwo[i];
        }
    }
    return amountOfSpeciesInBothSamples/sumOfAllOmus;
}
