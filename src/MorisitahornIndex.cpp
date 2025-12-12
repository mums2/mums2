//
// Created by Gregory Johnson on 4/25/25.
//
#include <numeric>
#include <cmath>
#include "DiversityMetrics/BetaDiversityCalculators/MorisitahornIndex.h"

double MorisitahornIndex::Calculate(const std::vector<std::vector<double>>& abundanceVectors) const {
    const std::vector<double>& sampleOne = abundanceVectors[0];
    const std::vector<double>& sampleTwo = abundanceVectors[1];
    const double totalSampleOne = std::accumulate(sampleOne.begin(), sampleOne.end(), 0.0);
    const double totalSampleTwo = std::accumulate(sampleTwo.begin(), sampleTwo.end(), 0.0);
    double summationOfSampleOneAndTwo = 0;
    double summationSampleOne = 0;
    double summationSampleTwo = 0;
    for(std::size_t i = 0; i < sampleOne.size(); i++) {
        const double currentSampleOne = sampleOne[i];
        const double currentSampleTwo = sampleTwo[i];
        summationOfSampleOneAndTwo += (currentSampleOne / totalSampleOne)
        * (currentSampleTwo / totalSampleTwo);

        summationSampleOne += std::pow(currentSampleOne/totalSampleOne, 2);
        summationSampleTwo += std::pow(currentSampleTwo/totalSampleTwo, 2);
    }
    return 1 - 2 * (summationOfSampleOneAndTwo / (summationSampleOne + summationSampleTwo));
}
