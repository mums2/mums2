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
                                    const std::vector<int64_t>& abundanceRanges,
                                    const int64_t size, const int64_t threshold) {
    if(eligibleIndexes.empty()) return abund;
    // O(N)
    int64_t sum = std::accumulate(abund.begin(), abund.end(), 0LL);
    if (sum <= size) return abund;

    std::set<RandomizationMethods::CountIndexPair> vals;
    // O(N)
    for(int i = 0; i < abundanceRanges.size(); i++) {
        vals.insert({i - 1,abundanceRanges[i]});
    }
    int64_t grandTotal = 0;
    int64_t incrementer = size;
    const size_t abundSize = abund.size();
    // const size_t eligibleIndexSize = eligibleIndexes.size();
    std::vector<int64_t> counter(abundSize, 0);
    // (O(NlogN))
    std::vector<size_t> indexes = RandomizationMethods::GetRandomVectorWithoutReplacement(vals,
        size, sum);
    // O(N)
    while(grandTotal <= size) {
        const int64_t currentSize = incrementer;
        for(const auto& index : indexes) {
            const auto eligibleIndex = eligibleIndexes[index];
            abund[eligibleIndex]--;
            counter[eligibleIndex]++;
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
        indexes = RandomizationMethods::GetRandomVectorWithoutReplacement(vals,
                                                                          incrementer, sum);
    }
    //Set counter size to 0 if they do not pass the threshold
    // O(N)
    for(auto& num : counter) {
        if(num < threshold) {
            num = 0;
        }
    }
    return counter;
}

std::vector<int64_t> Rarefaction::Rarefy2(const std::vector<int64_t>& abundance,
    const std::vector<int64_t>& eligibleIndex,
    std::vector<int64_t>& eligibleAbundances,
    const int64_t size, const int64_t threshold) {

    if(eligibleIndex.empty()) return abundance;
    const int64_t sum = std::accumulate(abundance.begin(), abundance.end(), 0LL);
    if(sum <= size) return abundance;
    const size_t vectorSize = eligibleAbundances.size();
    int64_t grandTotal = 0;
    int64_t incrementer = size;
    std::vector<int64_t> counter(abundance.size(), 0);
    std::vector<size_t> indexes = RandomizationMethods::GetRandomIndexVector(eligibleAbundances, incrementer,
                                                                             sum, vectorSize);
    while(grandTotal <= size) {
        for(const auto& index : indexes) {
            counter[eligibleIndex[index]]++;
        }
        for(const auto& abund : counter) {
            const auto value = abund;
            if(value >= threshold) {
                grandTotal += value;
            }
        }
        if(grandTotal >= size) {
            break;
        }

        incrementer = size - grandTotal;
        indexes = RandomizationMethods::GetRandomIndexVector(eligibleAbundances, incrementer,
                                                             sum, vectorSize);
        grandTotal = 0;
    }
    for(auto& value : counter) {
        if(value < threshold)
            value = 0;
    }
    return counter;
}

