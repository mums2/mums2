//
// Created by gregj on 1/14/2025.
//

#include "Rarefy/Rarefaction.h"

#include <numeric>
#include <unordered_map>

#include "Math/ParallelRandomNumberSitmo.h"

std::vector<uint64_t> Rarefaction::Rarefy(const std::vector<uint64_t>& abundance,
                                          const std::vector<uint64_t>& eligibleIndex,
                                          const std::vector<uint64_t>& abundancesRanges,
                                          ParallelRandomNumberSitmo& rngEngine,
                                          const uint64_t size, const uint64_t sum,
                                          const uint64_t threshold) {

    if(eligibleIndex.empty()) return abundance;
    uint64_t aboveThresholdSum = 0;
    for(const auto& abund : abundance) {
        if(abund >= threshold)
            aboveThresholdSum += abund;
    }
    if(sum <= size || aboveThresholdSum <= size)
        return abundance;

    const size_t vectorSize = abundance.size();
    uint64_t grandTotal = 0;
    uint64_t incrementer = size;

    std::vector<uint64_t> counter(vectorSize, 0);
    std::unordered_map<size_t, size_t> indexSwap;
    size_t currentIndex = 0;
    while(grandTotal <= size) {

        const auto maxValue = incrementer + currentIndex;
        std::vector<size_t> randomNumbers(maxValue - currentIndex);
        for(size_t i = currentIndex, index = 0; i < maxValue; ++i) {
            randomNumbers[index++] = rngEngine.NextRandomValue(i, sum);

        }
        for(size_t i = 0; i < randomNumbers.size(); i++) {
            auto randomValue = randomNumbers[i];
            if (indexSwap.find(randomValue) != indexSwap.end()) {
                // Set the random number to the next index
                size_t currentRandomValue = randomValue;
                randomValue = indexSwap[randomValue];
                if (indexSwap.find(randomValue) != indexSwap.end())
                    randomValue = indexSwap[randomValue];
                indexSwap[currentRandomValue] = i;
            }
            else
                indexSwap[randomValue] = i;

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


