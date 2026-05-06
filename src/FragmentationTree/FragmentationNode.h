//
// Created by gregj on 5/14/2025.
//

#ifndef FRAGMENTATIONNODE_H
#define FRAGMENTATIONNODE_H
#include <utility>

#include "../Chemicals/MolecularFormula/MolecularFormula.h"

struct FragmentationNode {
    FragmentationNode() = default;
    FragmentationNode(const int color, const double score,
        MolecularFormula formula)
        : color(color),
          score(score),
          formula(std::move(formula)) {
    }

    int color{};
    double score{};
    MolecularFormula formula;
    size_t amountOfDistinctColors = 0;
};

struct CompareFragmentationNodes {
    bool operator()(FragmentationNode const& s1, FragmentationNode const & s2) const {
        if (s1.amountOfDistinctColors < s2.amountOfDistinctColors) return false;
        if (s1.amountOfDistinctColors > s2.amountOfDistinctColors) return true;
        return s1.score >= s2.score;
    }
};
#endif //FRAGMENTATIONNODE_H
