//
// Created by Gregory Johnson on 12/9/24.
//

#include "DiversityMetrics/AlphaDiversityCalculators/ShannonDiversityIndex.h"

#include <numeric>

double ShannonDiversityIndex::Calculate(const std::vector<double> &abundances) const {
    const double sum = std::accumulate(abundances.begin(), abundances.end(), 0.0);
    double speciesProportion = 0;
    for (const auto &abundance : abundances) {
        const double result = abundance / sum;
        speciesProportion += result * result;
    }
    return -std::log(speciesProportion);
}
