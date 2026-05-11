//
// Created by gregj on 5/11/2026.
//

#ifndef SPECIESRICHNESSDIVERSITY_H
#define SPECIESRICHNESSDIVERSITY_H
#include "../DiversityCalculator.h"


class SpeciesRichnessDiversity : public DiversityCalculator {
public:
    SpeciesRichnessDiversity() = default;
    double Calculate(const std::vector<std::vector<double>> &) const override;
    ~SpeciesRichnessDiversity() override = default;
};



#endif //SPECIESRICHNESSDIVERSITY_H
