//
// Created by Gregory Johnson on 12/12/24.
//

#ifndef BRAYCURTISDISSIMILARITY_H
#define BRAYCURTISDISSIMILARITY_H
#include "../DiversityCalculator.h"


class BrayCurtisDissimilarity final : public DiversityCalculator{
public:
    BrayCurtisDissimilarity() = default;
    ~BrayCurtisDissimilarity() override = default;

    double Calculate(const std::vector<std::vector<double>> &) const override;
};



#endif //BRAYCURTISDISSIMILARITY_H
