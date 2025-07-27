//
// Created by Gregory Johnson on 4/25/25.
//

#include "DiversityMetrics/BetaDiversityCalculators/ThetaycDissimilarity.h"

double ThetaycDissimilarity::Calculate(const std::vector<std::vector<double>>& abundanceVectors) const {
    const std::vector<double>& sampleOne = abundanceVectors[0]; // always the same size
    const std::vector<double>& sampleTwo = abundanceVectors[1];
    double totalSampleOne = std::accumulate(sampleOne.begin(), sampleOne.end(), 0.0);
    double totalSampleTwo = std::accumulate(sampleTwo.begin(), sampleTwo.end(), 0.0);
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
