//
// Created by Gregory Johnson on 12/20/24.
//

#include "DiversityMetricFactory.h"
#include "AlphaDiversityCalculators/ShannonDiversityIndex.h"
#include "AlphaDiversityCalculators/SimpsonsDiversityIndex.h"
#include "BetaDiversityCalculators/BrayCurtisDissimilarity.h"

DiversityCalculator * DiversityMetricFactory::ChooseDiversityMetricBasedOnName(const std::string &metricIndex) {
    std::transform(metricIndex.begin(), metricIndex.end(), metricIndex.begin(), tolower);
    if(metricIndex == "shannon") {
        return new ShannonDiversityIndex();
    }
    else if(metricIndex == "simpson") {
        return new SimpsonsDiversityIndex();
    }
    else if(metricIndex == "bray") {
        return new BrayCurtisDissimilarity();
    }
    return nullptr;
}
