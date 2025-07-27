//
// Created by Gregory Johnson on 12/9/24.
//

#ifndef SHANNONDIVERSITYINDEX_H
#define SHANNONDIVERSITYINDEX_H
#include "../DiversityCalculator.h"
#include <numeric>
#include <cmath>

class ShannonDiversityIndex final : public DiversityCalculator{
public:
    ShannonDiversityIndex() = default;
    ~ShannonDiversityIndex() override = default;
    double Calculate(const std::vector<std::vector<double>>&) const override;
};



#endif //SHANNONDIVERSITYINDEX_H
