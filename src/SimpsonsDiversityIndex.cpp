//
// Created by gregj on 12/3/2024.
//

#include "DiversityMetrics/AlphaDiversityCalculators/SimpsonsDiversityIndex.h"

double SimpsonsDiversityIndex::Calculate(const std::vector<double>& abundances) const {

    double sumOfParticularSpecies = 0;
    double sumOfSpeciesInTotalPopulation = 0;
    for (const auto& abundance : abundances) {
        sumOfParticularSpecies += abundance * (abundance - 1);
        sumOfSpeciesInTotalPopulation += abundance;
    }
    return 1 - sumOfParticularSpecies / (sumOfSpeciesInTotalPopulation * (sumOfSpeciesInTotalPopulation - 1));
}

