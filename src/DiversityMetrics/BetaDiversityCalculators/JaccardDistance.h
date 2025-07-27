//
// Created by Gregory Johnson on 4/24/25.
//

#ifndef JACCARDDISTANCE_H
#define JACCARDDISTANCE_H
#include "../DiversityCalculator.h"


class JaccardDistance final : public DiversityCalculator{
public:
    JaccardDistance() = default;
    ~JaccardDistance() override = default;
    double Calculate(const std::vector<std::vector<double>>&) const override;
};



#endif //JACCARDDISTANCE_H
// https://www.researchgate.net/profile/Raimundo-Real/publication/239604848_The_Probabilistic_Basis_of_Jaccard%27s_Index_of_Similarity/links/0c9605268d8ff04ab1000000/The-Probabilistic-Basis-of-Jaccards-Index-of-Similarity.pdf