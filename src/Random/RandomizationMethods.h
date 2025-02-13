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
    static std::vector<size_t> GetRandomVectorWithoutReplacement(
        const std::vector<int64_t> &weightRanges, int64_t sizeToPull, int64_t sum);
private:
    struct CountIndexPair {
        int64_t index;
        int64_t abundance;
        bool operator<(const CountIndexPair& other) const {
            return abundance < other.abundance;
        }
    };
};



#endif //RANDOMIZATIONMETHODS_H
