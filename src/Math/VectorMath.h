//
// Created by Gregory Johnson on 5/27/25.
//

#ifndef VECTORMATH_H
#define VECTORMATH_H
#include <Rcpp.h>
class VectorMath {
    public:
    static double CosineScore(const std::vector<double>& x, const std::vector<double>& y) {
        double dotValue = 0;
        double magnitudeOne = 0;
        double magnitudeTwo = 0;
        for (size_t i = 0; i < x.size(); i++) {
            dotValue += x[i] * y[i];
            magnitudeOne += std::pow(x[i], 2);
            magnitudeTwo += std::pow(y[i], 2);
        }
        return dotValue / std::sqrt(magnitudeOne * magnitudeTwo);
    }
};
#endif //VECTORMATH_H
