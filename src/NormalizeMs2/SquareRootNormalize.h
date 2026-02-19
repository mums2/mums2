//
// Created by Gregory Johnson on 2/14/24.
//

#ifndef SQUAREROOTNORMALIZE_H
#define SQUAREROOTNORMALIZE_H
#include <complex>
#include <numeric>
#include <vector>



class SquareRootNormalize {
public:
    static void Normalize(std::vector<double>& value)
    {
        const size_t count = value.size();
        double summation = std::accumulate(value.begin(), value.end(), 0.0);
        summation = std::sqrt(summation);
        for(size_t i = 0; i < count; i++)
        {
            value[i] = std::sqrt(value[i])/ summation;
        }
    }
};



#endif // SQUAREROOTNORMALIZE_H
