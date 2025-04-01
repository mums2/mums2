//
// Created by gregj on 1/14/2025.
//

#include "Rarefy/Rarefaction.h"

#include <numeric>
#include <unordered_map>

std::vector<uint32_t> Rarefaction::Rarefy(const std::vector<uint32_t>& abundance,
    const std::vector<uint32_t>& eligibleIndex,
    std::vector<uint32_t>& availableIndexValues,
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
    std::deque<std::pair<size_t, size_t>> indexSwap;

    size_t currentIndex = 0;
    while(grandTotal <= size) {

        const auto maxValue = incrementer + currentIndex;
        for(size_t i = currentIndex; i < maxValue; i++) {
            const auto randomIndex = static_cast<size_t>(R::runif(i, sum));
            const auto index = availableIndexValues[randomIndex];
            std::swap(availableIndexValues[randomIndex], availableIndexValues[i]);
            indexSwap.emplace_front(i, randomIndex);
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
    while(!indexSwap.empty()) {
        const auto pair = indexSwap.front();
        indexSwap.pop_front();
        std::swap(availableIndexValues[pair.first],availableIndexValues[pair.second]);
    }
    return counter;
}

