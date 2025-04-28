//
// Created by Gregory Johnson on 4/25/25.
//

#ifndef THETAYCDISSIMILARITY_H
#define THETAYCDISSIMILARITY_H
#include "../DiversityCalculator.h"


class ThetaycDissimilarity final : public DiversityCalculator{
public:
    ThetaycDissimilarity() = default;
    ~ThetaycDissimilarity() override = default;
    double Calculate(const Rcpp::List &) const override;
};



#endif //THETAYCDISSIMILARITY_H
