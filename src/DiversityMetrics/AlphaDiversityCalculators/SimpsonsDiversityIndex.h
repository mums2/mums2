//
// Created by gregj on 12/3/2024.
//

#ifndef SIMPSONSDIVERSITYINDEX_H
#define SIMPSONSDIVERSITYINDEX_H
#include "../DiversityCalculator.h"


class SimpsonsDiversityIndex final : public DiversityCalculator {
public:
    double Calculate(const std::vector<double> &) override;
    ~SimpsonsDiversityIndex() override = default;
};



#endif //SIMPSONSDIVERSITYINDEX_H
