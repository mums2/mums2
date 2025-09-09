//
// Created by gregj on 1/14/2025.
//

#ifndef RAREFACTION_H
#define RAREFACTION_H
#include <cstdint>
#include <vector>
#include <Rcpp.h>
#include "../Math/ParallelRandomNumberSitmo.h"

class Rarefaction {
public:
    static std::vector<uint64_t> Rarefy(const std::vector<uint64_t> &abundance, const std::vector<uint64_t> &eligibleIndex,
                                        const std::vector<uint64_t> &abundancesRanges,
                                        ParallelRandomNumberSitmo& rngEngine,
                                        uint64_t size, uint64_t sum,
                                        uint64_t threshold);
};



#endif //RAREFACTION_H
