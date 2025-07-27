//
// Created by Gregory Johnson on 1/16/25.
//

#ifndef DIVERSITY_H
#define DIVERSITY_H
#include <Rcpp.h>
#include "../DataStructures/CppMatrix.h"
class Diversity {
public:
    virtual ~Diversity();
    virtual CppMatrix CalculateDiversity(const CppMatrix& communityMatrix,
        const std::string& index) = 0;
};



#endif //DIVERSITY_H
