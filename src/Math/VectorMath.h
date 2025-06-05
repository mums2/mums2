//
// Created by Gregory Johnson on 5/27/25.
//

#ifndef VECTORMATH_H
#define VECTORMATH_H
#include <Rcpp.h>
class VectorMath {
    public:
    static double CosineScore(Rcpp::NumericVector x, Rcpp::NumericVector y) {
        double dotValue = Rcpp::sum(x * y);
        double magnitudeOne = std::sqrt(Rcpp::sum(Rcpp::pow(x, 2)));
        double magnitudeTwo = std::sqrt(Rcpp::sum(Rcpp::pow(y, 2)));
        return dotValue / (magnitudeOne * magnitudeTwo);
    }
};
#endif //VECTORMATH_H
