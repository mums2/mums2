//
// Created by Gregory Johnson on 1/16/25.
//

#ifndef BETADIVERSITY_H
#define BETADIVERSITY_H

#include <Rcpp.h>

#include "../Diversity.h"


class BetaDiversity final : public Diversity {
public:
    Rcpp::NumericMatrix CalculateDiversity(const Rcpp::NumericMatrix& communityMatrix,
        const std::string& index) override;
    ~BetaDiversity() override = default;
};



#endif //BETADIVERSITY_H
