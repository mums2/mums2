//
// Created by Gregory Johnson on 12/9/24.
//
#include <numeric>
#include "DiversityMetrics/AlphaDiversityCalculators/ShannonDiversityIndex.h"


double ShannonDiversityIndex::Calculate(const std::vector<std::vector<double>>& abundanceVectors) const {
    const std::vector<double>& abundanceList = abundanceVectors[0];
    const double sum = std::accumulate(abundanceList.begin(), abundanceList.end(), 0.0);
    double speciesProportion = 0;
    for (const auto &abundance : abundanceList) {
        const double result = abundance / sum;
        speciesProportion += result * result;
    }
    return -std::log(speciesProportion);
}
