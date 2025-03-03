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
    std::vector<int64_t> Rarefy(const std::vector<int64_t> &abundance, const std::vector<int64_t> &eligibleIndex,
                                 std::vector<int64_t> &availableIndexValues, int64_t size, int64_t sum, int64_t
                                 threshold);
};



#endif //RAREFACTION_H
