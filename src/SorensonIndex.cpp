//
// Created by Gregory Johnson on 4/25/25.
//

#include "DiversityMetrics/BetaDiversityCalculators/SorensonIndex.h"

double SorensonIndex::Calculate(const std::vector<std::vector<double>>& abundanceVectors) const {
    const std::vector<double>& sampleOne = abundanceVectors[0];
    const std::vector<double>& sampleTwo = abundanceVectors[1];
    double sharedOmus = 0;
    double sampleOneOmus = 0;
    double sampleTwoOmus = 0;
    for(size_t i = 0; i < sampleOne.size(); i++) {
        const double firstCurrentOmu = sampleOne[i];
        const double secondCurrentOmu = sampleTwo[i];
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
