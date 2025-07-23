//
// Created by gregj on 1/14/2025.
//

#include "Rarefy/Rarefaction.h"

#include <numeric>
#include <unordered_map>

std::vector<uint32_t> Rarefaction::Rarefy(const std::vector<uint32_t>& abundance,
    const std::vector<uint32_t>& eligibleIndex,
    std::vector<uint32_t>& abundancesRanges,
    const uint32_t size, const uint32_t sum,
    const uint32_t threshold) {

    if(eligibleIndex.empty()) return abundance;
    uint32_t aboveThresholdSum = 0;
    for(const auto& abund : abundance) {
        if(abund >= threshold)
            aboveThresholdSum += abund;
    }
    if(sum <= size || aboveThresholdSum <= size)
        return abundance;

    const size_t vectorSize = abundance.size();
    uint32_t grandTotal = 0;
    uint32_t incrementer = size;

    std::vector<uint32_t> counter(vectorSize, 0);
    std::unordered_map<size_t, size_t> indexSwap;
    size_t currentIndex = 0;
    size_t currentIndexSwap = 0;
    while(grandTotal <= size) {

        const auto maxValue = incrementer + currentIndex;
        for(size_t i = currentIndex; i < maxValue; i++) {
            auto randomValue = static_cast<size_t>(R::runif(i, sum));
            if (indexSwap.find(randomValue) != indexSwap.end()) {
                randomValue = indexSwap[randomValue];
                // Take the number that was swapped with the original number
                // And give it its own value;
                size_t tempNumber = indexSwap[randomValue];
                indexSwap[tempNumber] = currentIndexSwap++;
            }
            else
                indexSwap[randomValue] = currentIndexSwap++;

            const auto index = std::lower_bound(abundancesRanges.begin(),
                abundancesRanges.end(), randomValue) - abundancesRanges.begin();
            counter[eligibleIndex[index]]++;
        }
        if(currentIndex <= 0) currentIndex += size;
        else currentIndex += incrementer;

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
        grandTotal = 0;
    }
    for(auto& value : counter) {
        if(value < threshold)
            value = 0;
    }
    return counter;
}


