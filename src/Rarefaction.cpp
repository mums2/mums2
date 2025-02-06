//
// Created by gregj on 1/14/2025.
//

#include "Rarefy/Rarefaction.h"

#include <numeric>
#include <unordered_map>

#include "Random/RandomizationMethods.h"

std::vector<int64_t> Rarefaction::Rarefy(const std::vector<int> &feature, std::vector<int64_t> &abund,
                                    const std::vector<int64_t>& eligibleIndexes,
                                    const std::vector<int64_t>& eligibleAbundances,
                                    std::vector<int64_t>& abundanceRanges,
                                    const int64_t size, const int64_t threshold) {
    if(eligibleIndexes.empty()) return abund;
    int64_t sum = std::accumulate(abund.begin(), abund.end(), 0LL);
    if (sum <= size) return abund;
    int64_t grandTotal = 0;
    int64_t incrementer = size;
    const size_t abundSize = abund.size();
    // const size_t eligibleIndexSize = eligibleIndexes.size();
    std::vector<int64_t> counter(abundSize, 0);
    std::vector<size_t> indexes = RandomizationMethods::GetRandomVectorWithoutReplacement(abundanceRanges,
        size, sum);
    while(grandTotal <= size) {
        const int64_t currentSize = incrementer;
        for(const auto& index : indexes) {
            abund[eligibleIndexes[index]]--;
            counter[eligibleIndexes[index]]++;
        }
        sum -= currentSize;
        for(size_t i = 0; i < abundSize; i++) {
            const auto value = counter[i];
            if(value < threshold)
                continue;
            grandTotal += value;
        }
        if(grandTotal >= size) {
            break;
        }
        incrementer = size - grandTotal;
        grandTotal = 0;
        indexes = RandomizationMethods::GetRandomVectorWithoutReplacement(abundanceRanges,
                                                                          incrementer, sum);
    }
    //Set counter size to 0 if they do not pass the threshold
    for(auto& num : counter) {
        if(num < threshold) {
            num = 0;
        }
    }
    return counter;
}

