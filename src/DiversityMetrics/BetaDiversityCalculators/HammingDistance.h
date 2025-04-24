//
// Created by Gregory Johnson on 4/24/25.
//

#ifndef HAMMINGDISTANCE_H
#define HAMMINGDISTANCE_H
#include "../DiversityCalculator.h"


class HammingDistance final : public DiversityCalculator {
public:
    HammingDistance() = default;
    ~HammingDistance() override = default;
    double Calculate(const Rcpp::List&) const override;

};



#endif //HAMMINGDISTANCE_H
