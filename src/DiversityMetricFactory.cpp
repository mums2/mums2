//
// Created by Gregory Johnson on 12/20/24.
//

#include "DiversityMetrics/DiversityMetricFactory.h"
#include "DiversityMetrics/AlphaDiversityCalculators/ShannonDiversityIndex.h"
#include "DiversityMetrics/AlphaDiversityCalculators/SimpsonsDiversityIndex.h"
#include "DiversityMetrics/BetaDiversityCalculators/BrayCurtisDissimilarity.h"

DiversityCalculator* DiversityMetricFactory::ChooseDiversityMetricBasedOnName(const std::string &metricIndex) {

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
