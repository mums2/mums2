//
// Created by Gregory Johnson on 4/25/25.
//

#ifndef MORISITAHORNINDEX_H
#define MORISITAHORNINDEX_H
#include "../DiversityCalculator.h"


class MorisitahornIndex final : public DiversityCalculator{
public:
    MorisitahornIndex() = default;
    ~MorisitahornIndex() override = default;
    double Calculate(const std::vector<std::vector<double>>&) const override;
};



#endif //MORISITAHORNINDEX_H
