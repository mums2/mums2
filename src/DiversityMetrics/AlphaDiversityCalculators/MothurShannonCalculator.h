//
// Created by gregj on 1/9/2025.
//

#ifndef MOTHURSHANNONCALCULATOR_H
#define MOTHURSHANNONCALCULATOR_H
#include "../DiversityCalculator.h"


class MothurShannonCalculator final : DiversityCalculator {
public:
    ~MothurShannonCalculator() override;
    double Calculate(const Rcpp::List &) const override;
};



#endif //MOTHURSHANNONCALCULATOR_H
