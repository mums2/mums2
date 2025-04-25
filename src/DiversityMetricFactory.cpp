//
// Created by Gregory Johnson on 12/20/24.
//

#include "DiversityMetrics/DiversityMetricFactory.h"

#include "DiversityMetrics/AlphaDiversityCalculators/AlphaDiversity.h"
#include "DiversityMetrics/AlphaDiversityCalculators/ShannonDiversityIndex.h"
#include "DiversityMetrics/AlphaDiversityCalculators/SimpsonsDiversityIndex.h"
#include "DiversityMetrics/BetaDiversityCalculators/BetaDiversity.h"
#include "DiversityMetrics/BetaDiversityCalculators/BrayCurtisDissimilarity.h"
#include "DiversityMetrics/BetaDiversityCalculators/HammingDistance.h"
#include "DiversityMetrics/BetaDiversityCalculators/JaccardDistance.h"
#include "DiversityMetrics/BetaDiversityCalculators/MorisitahornIndex.h"
#include "DiversityMetrics/BetaDiversityCalculators/SorensonIndex.h"
#include "DiversityMetrics/BetaDiversityCalculators/ThetaycDissimilarity.h"

DiversityCalculator* DiversityMetricFactory::ChooseDiversityMetricBasedOnName(const std::string &metricIndex) {

    if(metricIndex == "shannon" ) {
        return new ShannonDiversityIndex();
    }
    if(metricIndex == "simpson") {
        return new SimpsonsDiversityIndex();
    }
    if(metricIndex == "bray") {
        return new BrayCurtisDissimilarity();
    }
    if(metricIndex == "jaccard") {
        return new JaccardDistance();
    }
    if(metricIndex == "hamming") {
        return new HammingDistance();
    }
    if(metricIndex == "soren") {
        return new SorensonIndex();
    }
    if(metricIndex == "morista") {
        return new MorisitahornIndex();
    }
    if(metricIndex == "thetayc") {
        return new ThetaycDissimilarity();
    }
    return nullptr;
}

Diversity* DiversityMetricFactory::ChooseDiversityBasedOnIndex(const std::string& index) {
    std::unordered_map<std::string, std::string> indexMap;
    indexMap["shannon"] = "alpha";
    indexMap["simpson"] = "alpha";
    indexMap["bray"] = "beta";
    indexMap["jaccard"] = "beta";
    indexMap["hamming"] = "beta";
    indexMap["soren"] = "beta";
    indexMap["morista"] = "beta";
    indexMap["thetayc"] = "beta";

    if(indexMap[index] == "alpha") {
        return new AlphaDiversity();
    }
    if(indexMap[index] == "beta") {
        return new BetaDiversity();
    }
    return nullptr;
}
