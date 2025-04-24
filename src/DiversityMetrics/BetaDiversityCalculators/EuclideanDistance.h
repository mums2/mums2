//
// Created by Gregory Johnson on 4/24/25.
//

#ifndef EUCLIDEANDISTANCE_H
#define EUCLIDEANDISTANCE_H
#include "../DiversityCalculator.h"


class EuclideanDistance final : DiversityCalculator {
public:
    EuclideanDistance() = default;
    ~EuclideanDistance() override = default;
    double Calculate(const Rcpp::List &) const override;
};



#endif //EUCLIDEANDISTANCE_H
