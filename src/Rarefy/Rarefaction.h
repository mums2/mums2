//
// Created by gregj on 1/14/2025.
//

#ifndef RAREFACTION_H
#define RAREFACTION_H
#include <cstdint>
#include <vector>
#include <Rcpp.h>

class Rarefaction {
public:
    std::vector<int64_t> Rarefy(const std::vector<int> &feature, std::vector<int64_t> &abund, const std::vector<int64_t> &eligibleIndexes,
                                const std::vector<int64_t> &eligibleAbundances, const std::vector<int64_t> &abundanceRanges, int64_t size, int64_t threshold);
    std::vector<int64_t> Rarefy2(const std::vector<int64_t> &, const std::vector<int64_t> &,
                                 std::vector<int64_t> &shuffledVector, const std::vector<int64_t> &eligibleRanges, int64_t size, int64_t sum, int64_t
                                 threshold);
};



#endif //RAREFACTION_H
