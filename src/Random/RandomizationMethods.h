//
// Created by gregj on 1/14/2025.
//

#ifndef RANDOMIZATIONMETHODS_H
#define RANDOMIZATIONMETHODS_H
#include <cstdint>
#include <vector>


class RandomizationMethods final {
public:
    struct CountIndexPair {
        int64_t index;
        int64_t abundance;
        bool operator<(const CountIndexPair& other) const {
            return abundance < other.abundance;
        }
    };
    RandomizationMethods() = default;
    ~RandomizationMethods() = default;
    static std::vector<size_t> GetRandomVectorWithoutReplacement(
        std::set<CountIndexPair> &weightRanges, int64_t sizeToPull, int64_t sum);
    static std::vector<size_t> GetRandomIndexVector(std::vector<int64_t> &weights,
                                                    int64_t sizeToPull, int64_t sum, size_t vectorSize);

};



#endif //RANDOMIZATIONMETHODS_H
