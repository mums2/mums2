//
// Created by Gregory Johnson on 1/16/25.
//

#ifndef BETADIVERSITY_H
#define BETADIVERSITY_H

#include <Rcpp.h>

#include "../Diversity.h"
#include "../../DataStructures/CppMatrix.h"

class BetaDiversity final : public Diversity {
public:
    CppMatrix CalculateDiversity(const CppMatrix& communityMatrix,
        const std::string& index) override;
    ~BetaDiversity() override = default;
};



#endif //BETADIVERSITY_H
