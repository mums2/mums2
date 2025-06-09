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
    static std::vector<uint32_t> Rarefy(const std::vector<uint32_t> &abundance, const std::vector<uint32_t> &eligibleIndex,
                                        std::vector<uint32_t> &availableIndexValues, uint32_t size, uint32_t sum, uint32_t
                                        threshold);
};



#endif //RAREFACTION_H
