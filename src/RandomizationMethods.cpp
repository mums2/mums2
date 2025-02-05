//
// Created by gregj on 1/14/2025.
//
#include <Rcpp.h>
#include "Random/RandomizationMethods.h"
#include "Random/ReservoirPairs.h"
#include <complex>
#include <algorithm>
size_t RandomizationMethods::GetRandomNumberIndex(const std::vector<int64_t> &weightedToPull, const size_t vectorSize,
    const int64_t sum) {
    // O(N) worse case O(log(N)) best case
    auto randomNumber = static_cast<int64_t>(R::runif(0, static_cast<double>(sum)));
    for (size_t i = 0; i < vectorSize; i++) {
        const int64_t currentWeight = weightedToPull[i];
        if (randomNumber < currentWeight) {
            return i;
        }
        randomNumber -= currentWeight;
    }
    return 0;
}

std::vector<size_t> RandomizationMethods::GetRandomVectorWithoutReplacement(const std::vector<int64_t> &weights,
    const int64_t sizeToPull, const int64_t sum) {
    std::vector<size_t> indexes(sizeToPull);
    std::vector<int64_t> ranges(weights.size());
    ranges[0] = weights[0];
    for(int i = 1; i < ranges.size(); i++) {
        ranges[i] = ranges[i - 1] + weights[i];
    }
    for(int i = 0; i < sizeToPull; i++) {
        const auto randomNum = static_cast<int64_t>(R::runif(0, static_cast<double>(sum - i)));
        const size_t index = std::lower_bound(ranges.begin(),
            ranges.end(), randomNum) - ranges.begin();
        indexes[i] = index;
        ranges[index]--;
    }
    return indexes;
}

// If batched -> Worse: O(N^2)  Best: O(Nlog(N))
// 66, 167