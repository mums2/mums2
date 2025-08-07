//
// Created by Gregory Johnson on 12/12/24.
//

#include "DiversityMetrics/BetaDiversityCalculators/BrayCurtisDissimilarity.h"

#include <numeric>

// Assume there are only two vectors for this calculation
double BrayCurtisDissimilarity::Calculate(const std::vector<std::vector<double>>& abundanceVectors) const {
    const std::vector<double>& sampleOne = abundanceVectors[0];
    const std::vector<double>& sampleTwo = abundanceVectors[1];
    const size_t size = sampleOne.size();
    const double abundanceSum = std::accumulate(sampleOne.begin(), sampleOne.end(), 0.0) +
        std::accumulate(sampleTwo.begin(), sampleTwo.end(), 0.0);
    double lesserAbundances = 0;
    for(size_t i = 0; i < size; i++) {
        if(sampleOne[i] < sampleTwo[i]) {
            lesserAbundances += sampleOne[i];
            continue;
        }
        lesserAbundances += sampleTwo[i];
    }
    return 1 - 2 * lesserAbundances / abundanceSum;
}
