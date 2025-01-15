//
// Created by gregj on 1/14/2025.
//
#include <Rcpp.h>
#include "Random/RandomizationMethods.h"

size_t RandomizationMethods::GetRandomNumberIndex(const std::vector<int64_t> &weightedToPull, const int64_t sum) {
    auto randomNumber = static_cast<int64_t>(R::runif(0, static_cast<double>(sum)));
    for (size_t i = 0; i < weightedToPull.size(); i++) {
        if (randomNumber < weightedToPull[i]) {
            return i;
        }
        randomNumber -= weightedToPull[i];
    }
    return 0;
}
