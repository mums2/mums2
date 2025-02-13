//
// Created by gregj on 1/14/2025.
//
#include <Rcpp.h>
#include "Random/RandomizationMethods.h"
#include <algorithm>
// Have to deal with ranges overlapping...How to properly do this
std::vector<size_t> RandomizationMethods::GetRandomVectorWithoutReplacement(const std::vector<int64_t> &weightRanges,
    const int64_t sizeToPull, const int64_t sum) {
    std::set<CountIndexPair> vals;
    for(int i = 0; i < weightRanges.size(); i++) {
        vals.insert(CountIndexPair{i - 1,weightRanges[i]});
    }
    std::vector<size_t> indexes(sizeToPull);
    for(int i = 0; i < sizeToPull; i++) {
        const CountIndexPair randomNum{-1,
            static_cast<int64_t>(R::runif(1, static_cast<double>(sum - i)))};
        const auto val = vals.upper_bound(randomNum);
        const CountIndexPair updatedPair{val->index , val->abundance - 1};
        indexes[i] = val->index;
        vals.erase(val);
        vals.insert(updatedPair);
    }
    return indexes;
}

