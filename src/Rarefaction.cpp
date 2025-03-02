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
    std::vector<int64_t>& shuffledVector,
    const std::vector<int64_t>& eligibleRanges,
    const int64_t size, const int64_t sum,
    const int64_t threshold) {

    if(eligibleIndex.empty()) return abundance;
    if(sum <= size) return abundance;
    const size_t vectorSize = abundance.size();
    int64_t grandTotal = 0;
    int64_t incrementer = size;

    std::vector<int64_t> counter(vectorSize, 0);
    std::deque<std::pair<size_t, size_t>> indexSwap;
    // std::vector<size_t> indexes (shuffledVector.begin(), shuffledVector.begin()  + size);
    // std::vector<bool> hasChosen(shuffledVector.size(), false);
    size_t currentIndex = 0;
    while(grandTotal <= size) {
        // const auto startTime2 = std::chrono::steady_clock::now();
        const auto maxValue = incrementer + currentIndex;
        for(size_t i = currentIndex; i < maxValue; i++) {
            const auto randomIndex = static_cast<size_t>(R::runif(i, sum));
            const auto index = shuffledVector[randomIndex];
            std::swap(shuffledVector[randomIndex], shuffledVector[i]);
            indexSwap.emplace_front(i, randomIndex);
            // if(number >= sum) continue;

            // if(hasChosen[number]) continue;
            // hasChosen[number] = true;
            // const size_t index = std::upper_bound(eligibleRanges.begin(),
            //     eligibleRanges.end(), number) - eligibleRanges.begin();
            counter[eligibleIndex[index]]++;
        }
        // const auto endTime2 = std::chrono::steady_clock::now();
        // const auto time = std::chrono::duration_cast<std::chrono::microseconds>(endTime2-startTime2).count();
        // Rcpp::Rcout << "It took: " << time << " microseconds to finish an iteration of the loop" << std::endl;
        if(currentIndex <= 0) currentIndex += size;
        else currentIndex += incrementer;
        // for(const auto& number : indexes) {
        //     if(number >= sum) continue;
        //     if(hasChosen[number]) continue;
        //     hasChosen[number] = true;
        //     const size_t index = std::upper_bound(eligibleRanges.begin(),
        //         eligibleRanges.end(), number) - eligibleRanges.begin();
        //     counter[eligibleIndex[index - 1]]++;
        // }
        // may need to make an early exit
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
        // indexes = std::vector<size_t>(shuffledVector.begin() + currentIndex,
        //     shuffledVector.begin() + currentIndex + incrementer);
        // currentIndex += incrementer;
        grandTotal = 0;
    }
    for(auto& value : counter) {
        if(value < threshold)
            value = 0;
    }
    while(!indexSwap.empty()) {
        const auto pair = indexSwap.front();
        indexSwap.pop_front();
        std::swap(shuffledVector[pair.first],shuffledVector[pair.second]);
    }
    // for(const auto pair : indexSwap) {
    //     std::swap(shuffledVector[pair.first],shuffledVector[pair.second]);
    // }
    return counter;
}

