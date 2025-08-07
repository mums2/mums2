//
// Created by Gregory Johnson on 7/25/25.
//

#ifndef PARALLELRANDOMNUMBERSITMO_H
#define PARALLELRANDOMNUMBERSITMO_H
#include <vector>
#include <sitmo.h>
class ParallelRandomNumberSitmo {
public:
    // [[Rcpp::depends(sitmo)]]
    explicit ParallelRandomNumberSitmo(int seed): rngEngine(sitmo::prng(seed)) {}
    ParallelRandomNumberSitmo() = default;
    // [[Rcpp::depends(sitmo)]]
    size_t NextRandomValue(const double min = 0.0, const double max = 1.0) {
        double dis = max - min;
        return min + (static_cast<double>(rngEngine()) / (sitmo::prng::max())) * (dis);
    }

private:
    sitmo::prng rngEngine{};
};
#endif //PARALLELRANDOMNUMBERSITMO_H
