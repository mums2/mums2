//
// Created by gregj on 1/14/2025.
//
#include <Rcpp.h>
#include "Random/RandomizationMethods.h"
#include <algorithm>
// Have to deal with ranges overlapping...How to properly do this
std::vector<size_t> RandomizationMethods::GetRandomVectorWithoutReplacement(std::vector<int64_t> &weightRanges,
    const int64_t sizeToPull, const int64_t sum) {
    std::vector<size_t> indexes(sizeToPull);
    for(int i = 0; i < sizeToPull; i++) {
        const auto randomNum = static_cast<int64_t>(R::runif(0, static_cast<double>(sum - i)));
        const size_t index = (std::upper_bound(weightRanges.begin(),
            weightRanges.end(), randomNum) - weightRanges.begin()) - 1;
        indexes[i] = index;
        weightRanges[index]--;
    }
    return indexes;
}
