//
// Created by Gregory Johnson on 4/25/25.
//

#ifndef SORENSONINDEX_H
#define SORENSONINDEX_H
#include "../DiversityCalculator.h"


class SorensonIndex final : public DiversityCalculator{
public:
    SorensonIndex() = default;
    ~SorensonIndex() override = default;
    double Calculate(const Rcpp::List &) const override;
};



#endif //SORENSONINDEX_H
