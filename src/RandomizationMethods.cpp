//
// Created by gregj on 1/14/2025.
//
#include <Rcpp.h>
#include "Random/RandomizationMethods.h"
#include <algorithm>
#include <set>
// Have to deal with ranges overlapping...How to properly do this
std::vector<size_t> RandomizationMethods::GetRandomVectorWithoutReplacement(std::set<CountIndexPair> &weightRanges,
    const int64_t sizeToPull, const int64_t sum) {
    // std::set<CountIndexPair> vals;
    // for(int i = 0; i < weightRanges.size(); i++) {
    //     vals.insert(CountIndexPair{i - 1,weightRanges[i]});
    // }
    std::vector<size_t> indexes(sizeToPull);
    for(int i = 0; i < sizeToPull; i++) {
        const CountIndexPair randomNum{-1,
            static_cast<int64_t>(R::runif(1, static_cast<double>(sum - i)))};
        const auto val = weightRanges.upper_bound(randomNum);
        const CountIndexPair updatedPair{val->index , val->abundance - 1};
        indexes[i] = val->index;
        weightRanges.erase(val);
        weightRanges.insert(updatedPair);
    }
    return indexes;
}

std::vector<size_t> RandomizationMethods::GetRandomIndexVector(std::vector<int64_t> &weights, const int64_t sizeToPull,
    const int64_t sum, const size_t vectorSize) {
    std::vector<size_t> indexes(sizeToPull);
    for(size_t i = 0; i < sizeToPull; i++) {
        auto randomNumber = static_cast<int64_t>(R::runif(1, static_cast<double>(sum)));
        for (size_t j = 0; j < vectorSize; j++) {
            const int64_t currentWeight = weights[j];
            if(currentWeight <= 0) continue;
            if (randomNumber < currentWeight) {
                indexes[i] = j;
                weights[j]--;
                break;
            }
            randomNumber -= currentWeight;
        }
    }
    return indexes;
}

