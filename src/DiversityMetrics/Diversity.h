//
// Created by Gregory Johnson on 1/16/25.
//

#ifndef DIVERSITY_H
#define DIVERSITY_H
#include <Rcpp.h>

class Diversity {
public:
    virtual ~Diversity();
    virtual Rcpp::NumericMatrix CalculateDiversity(const Rcpp::NumericMatrix& communityMatrix,
        const std::string& index) = 0;
};



#endif //DIVERSITY_H
