//
// Created by Gregory Johnson on 2/14/24.
//

#ifndef SQUAREROOTNORMALIZE_H
#define SQUAREROOTNORMALIZE_H
#include <complex>
#include <vector>


class SquareRootNormalize {
public:
    inline std::vector<double> Normalize(const std::vector<double> value)
    {
        const size_t count = value.size();
        std::vector<double> new_values(value.size());
        double summation = 0;
        for(size_t i = 0; i < count; i++)
        {
            summation += value[i];
        }
        summation = std::sqrt(summation);
        for(size_t i = 0; i < count; i++)
        {
            new_values[i] = (std::sqrt(value[i])/ summation);
        }
        return new_values;
    }
};



#endif // SQUAREROOTNORMALIZE_H
