//
// Created by gregj on 5/14/2025.
//

#ifndef FRAGMENTATIONNODE_H
#define FRAGMENTATIONNODE_H
#include <utility>

#include "../Chemicals/MolecularFormula/MolecularFormula.h"

struct FragmentationNode {
    FragmentationNode() = default;
    FragmentationNode(const int color, const double score, const double sub_tree_score,
        MolecularFormula formula)
        : color(color),
          score(score),
          subTreeScore(sub_tree_score),
          formula(std::move(formula)) {
    }

    int color{};
    double score{};
    double subTreeScore{};
    MolecularFormula formula;
};

struct CompareFragmentationNodes {
    bool operator()(FragmentationNode const& s1, FragmentationNode const & s2) const {
        if (s1.color != 0) return false;
        return s1.subTreeScore > s2.subTreeScore;
    }
};
#endif //FRAGMENTATIONNODE_H
