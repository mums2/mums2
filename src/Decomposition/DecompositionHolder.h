//
// Created by gregj on 4/30/2026.
//

#ifndef DECOMPOSITIONHOLDER_H
#define DECOMPOSITIONHOLDER_H
#include "../RdisopHeaderFiles/composedelement.h"

struct DecompositionHolder {
    ims::ComposedElement element;
    double score;
};

struct CompareDecompositions {
    bool operator()(DecompositionHolder const& s1, DecompositionHolder const & s2) const {
        return s1.score > s2.score;
    }
};
#endif //DECOMPOSITIONHOLDER_H
