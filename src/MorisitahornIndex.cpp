//
// Created by Gregory Johnson on 4/25/25.
//

#include "DiversityMetrics/BetaDiversityCalculators/MorisitahornIndex.h"

double MorisitahornIndex::Calculate(const std::vector<std::vector<double>>& abundanceVectors) const {
    const std::vector<double>& sampleOne = abundanceVectors[0];
    const std::vector<double>& sampleTwo = abundanceVectors[1];
    double totalSampleOne = std::accumulate(sampleOne.begin(), sampleOne.end(), 0.0);
    double totalSampleTwo = std::accumulate(sampleTwo.begin(), sampleTwo.end(), 0.0);
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
