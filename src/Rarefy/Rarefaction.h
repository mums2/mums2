//
// Created by gregj on 1/14/2025.
//

#ifndef RAREFACTION_H
#define RAREFACTION_H
#include <cstdint>
#include <vector>
#include <Rcpp.h>
#include <sitmo.h>
class Rarefaction {
public:
    // [[Rcpp::depends(sitmo)]]
    static std::vector<uint32_t> Rarefy(const std::vector<uint32_t> &abundance, const std::vector<uint32_t> &eligibleIndex,
                                        const std::vector<uint32_t> &abundancesRanges, sitmo::prng& rngEngine,
                                        uint32_t size, uint32_t sum,
                                        uint32_t threshold);
};



#endif //RAREFACTION_H
