//
// Created by Gregory Johnson on 12/20/24.
//

#include "DiversityMetrics/DiversityMetricFactory.h"

#include "DiversityMetrics/AlphaDiversityCalculators/AlphaDiversity.h"
#include "DiversityMetrics/AlphaDiversityCalculators/ShannonDiversityIndex.h"
#include "DiversityMetrics/AlphaDiversityCalculators/SimpsonsDiversityIndex.h"
#include "DiversityMetrics/BetaDiversityCalculators/BetaDiversity.h"
#include "DiversityMetrics/BetaDiversityCalculators/BrayCurtisDissimilarity.h"
#include "DiversityMetrics/BetaDiversityCalculators/EuclideanDistance.h"
#include "DiversityMetrics/BetaDiversityCalculators/HammingDistance.h"
#include "DiversityMetrics/BetaDiversityCalculators/JaccardDistance.h"

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
    if(metricIndex == "euclidean") {
        return new EuclideanDistance();
    }
    if(metricIndex == "hamming") {
        return new HammingDistance();
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
    indexMap["euclidean"] = "beta";

    if(indexMap[index] == "alpha") {
        return new AlphaDiversity();
    }
    if(indexMap[index] == "beta") {
        return new BetaDiversity();
    }
    return nullptr;
}
