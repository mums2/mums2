//
// Created by gregj on 12/3/2024.
//

#ifndef DIVERSITYCALCULATOR_H
#define DIVERSITYCALCULATOR_H
#include <vector>


class DiversityCalculator {
public:

    virtual ~DiversityCalculator() = default;
    virtual double Calculate(const std::vector<std::vector<double>>&) const = 0;
};



#endif //DIVERSITYCALCULATOR_H
