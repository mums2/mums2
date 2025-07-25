//
// Created by gregj on 1/14/2025.
//

#include "Rarefy/Rarefaction.h"

#include <numeric>
#include <unordered_map>

#include "Math/ParallelRandomNumberSitmo.h"

std::vector<uint32_t> Rarefaction::Rarefy(const std::vector<uint32_t>& abundance,
                                          const std::vector<uint32_t>& eligibleIndex,
                                          const std::vector<uint32_t>& abundancesRanges,
                                          ParallelRandomNumberSitmo& rngEngine,
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
    std::mutex mutex;
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


