//
// Created by gregj on 1/14/2025.
//

#include "Rarefy/Rarefaction.h"

#include <numeric>
#include <unordered_map>

#include "Random/RandomizationMethods.h"

Rcpp::DataFrame Rarefaction::Rarefy(const std::vector<int> &feature, std::vector<int64_t> &abund,
                                    const int64_t size, const int threshold) {
    const int64_t sum = std::accumulate(abund.begin(), abund.end(), 0LL);
    int64_t grandTotal = 0;
    int64_t incrementer = size;
    const size_t abundSize = abund.size();
    std::unordered_map<int, int64_t> filtered;
    std::vector<int> counter(abundSize, 0);
    while(grandTotal <= size) {
        std::unordered_map<int, int64_t> finalMap;
        const int64_t currentSize = incrementer;
        for(int64_t i = 0; i < currentSize; i++) {
            const size_t index = RandomizationMethods::GetRandomNumberIndex(abund, sum - i);
            abund[index]--;
            counter[index]++;
        }
        int count = 0;
        for(const auto& feat: feature) {
            const int value = counter[count++];
            if(value < threshold)
                continue;
            grandTotal += value;
            finalMap[feat] = value;
        }
        if(grandTotal >= size) {
            filtered = finalMap;
            break;
        }
        incrementer = size - grandTotal;
        grandTotal = 0;
    }
    std::vector<int> rare_mz(filtered.size(), 0);
    std::vector<int64_t> rare_abund(filtered.size(), 0);
    int pos = 0;
    for (const auto & it : filtered) {
        rare_mz[pos] = it.first;
        rare_abund[pos] = it.second;
        pos ++;
    }

    // delete calculator;
    // This is where we compute alpha and beta diversity
    return Rcpp::DataFrame::create(Rcpp::Named("mz") = rare_mz,
                             Rcpp::Named("abund") = rare_abund);
}

