//
// Created by gregj on 5/11/2026.
//

#include "DiversityMetrics/AlphaDiversityCalculators/SpeciesRichnessDiversity.h"

#include <algorithm>

double SpeciesRichnessDiversity::Calculate(const std::vector<std::vector<double>>& abundanceVectors) const {
    const std::vector<double>& abundanceList = abundanceVectors[0];
    return static_cast<double>(std::count_if(abundanceList.begin(), abundanceList.end(),
        [](const double value) {
        return value != 0;
    }));
}
