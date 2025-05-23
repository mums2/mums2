//
// Created by gregj on 5/14/2025.
//

#ifndef FRAGMENTATIONNODE_H
#define FRAGMENTATIONNODE_H
#include "../Chemicals/MolecularFormula/MolecularFormula.h"

struct FragmentationNode {
    int color{};
    size_t index{};
    double score{};
    double subTreeScore{};
    MolecularFormula formula;
    std::list<int> parentIndexes;
};

struct CompareFragmentationNodes {
    bool operator()(FragmentationNode const& s1, FragmentationNode const & s2) const {
        return s1.score > s2.score;
    }
};
#endif //FRAGMENTATIONNODE_H
