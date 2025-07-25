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
    static sitmo::prng CreateRandomNumberGenerateSitmo(uint32_t seed) {
        return sitmo::prng(seed);
    }
    //
    // static size_t runif_sitmo(unsigned int n, double min = 0.0, double max = 1.0, uint32_t seed = 1) {
    //     std::vector<size_t> o(n);
    //
    //     // Create a prng engine
    //     sitmo::prng eng(seed);
    //     // Obtain the range between max and min
    //     double dis = max - min;
    //
    //     for(int i = 0; i < n; ++i) {
    //         // Sample from the RNG and divide it by the maximum value possible (can also use SITMO_RAND_MAX, which is 4294967295)
    //         // Apply appropriate scale (MAX-MIN)
    //         o[i] = min + (static_cast<double>(eng()) / (sitmo::prng::max())) * (dis);
    //     }
    //
    //     return o;
    // }
    // [[Rcpp::depends(sitmo)]]
    static size_t GetRandomValue(sitmo::prng& engine, double min = 0.0, double max = 1.0) {
        double dis = max - min;
        return min + (static_cast<double>(engine()) / (sitmo::prng::max())) * (dis);
    }
};
#endif //PARALLELRANDOMNUMBERSITMO_H
