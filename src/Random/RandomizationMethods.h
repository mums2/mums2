//
// Created by gregj on 1/14/2025.
//

#ifndef RANDOMIZATIONMETHODS_H
#define RANDOMIZATIONMETHODS_H
#include <cstdint>
#include <vector>


class RandomizationMethods final {
public:
    RandomizationMethods() = default;
    ~RandomizationMethods() = default;
    static size_t GetRandomNumberIndex(const std::vector<int64_t> &weightedToPull, size_t vectorSize, int64_t sum);
    static std::vector<size_t> GetRandomVectorWithoutReplacement(
        std::vector<int64_t> &weightRanges, int64_t sizeToPull, int64_t sum);
};



#endif //RANDOMIZATIONMETHODS_H
