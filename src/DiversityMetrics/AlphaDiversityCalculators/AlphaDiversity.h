//
// Created by Gregory Johnson on 1/16/25.
//

#ifndef ALPHADIVERSITY_H
#define ALPHADIVERSITY_H
#include <vector>
#include <Rcpp.h>
#include "../Diversity.h"
#include "../../DataStructures/CppMatrix.h"

class AlphaDiversity final : public Diversity {
public:
    CppMatrix CalculateDiversity(const CppMatrix& communityMatrix,
        const std::string& index) override;
};



#endif //ALPHADIVERSITY_H
