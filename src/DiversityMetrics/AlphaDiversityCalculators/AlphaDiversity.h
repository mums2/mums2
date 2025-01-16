//
// Created by Gregory Johnson on 1/16/25.
//

#ifndef ALPHADIVERSITY_H
#define ALPHADIVERSITY_H
#include <vector>
#include <Rcpp.h>

#include "../Diversity.h"


class AlphaDiversity : public Diversity {
public:
    Rcpp::NumericMatrix CalculateDiversity(const Rcpp::NumericMatrix& communityMatrix,
        const std::string& index) override;
};



#endif //ALPHADIVERSITY_H
