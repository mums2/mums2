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
    std::vector<uint> Rarefy(const std::vector<uint> &abundance, const std::vector<uint> &eligibleIndex,
                                 std::vector<uint> &availableIndexValues, uint size, uint sum, uint
                                 threshold);
};



#endif //RAREFACTION_H
