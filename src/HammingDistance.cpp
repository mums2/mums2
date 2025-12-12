//
// Created by Gregory Johnson on 4/24/25.
//
#include <numeric>
#include <cmath>
#include "DiversityMetrics/BetaDiversityCalculators/HammingDistance.h"

double HammingDistance::Calculate(const std::vector<std::vector<double>>& abundanceVectors) const {
    const std::vector<double>& sampleOne = abundanceVectors[0];
    const std::vector<double>& sampleTwo = abundanceVectors[1];
    const size_t size = abundanceVectors.size();
    std::vector<double> subtractedVector(size,0);
    for (size_t i = 0; i < size; i++) {
        subtractedVector[i] = std::abs(sampleOne[i] - sampleTwo[i]);
    }
    return std::accumulate(subtractedVector.begin(), subtractedVector.end(), 0.0);
}
