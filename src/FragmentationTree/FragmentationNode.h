//
// Created by gregj on 5/14/2025.
//

#ifndef FRAGMENTATIONNODE_H
#define FRAGMENTATIONNODE_H
#include <utility>

#include "../Chemicals/MolecularFormula/MolecularFormula.h"

struct FragmentationNode {
    FragmentationNode() = default;
    FragmentationNode(const int color, const int index, const double score, const double sub_tree_score,
        MolecularFormula formula)
        : color(color),
          index(index),
          score(score),
          subTreeScore(sub_tree_score),
          formula(std::move(formula)) {
    }

    int color{};
    int index{};
    double score{};
    double subTreeScore{};
    MolecularFormula formula;
    std::list<int> parentIndexes;
};

struct CompareFragmentationNodes {
    bool operator()(FragmentationNode const& s1, FragmentationNode const & s2) const {
        if (s1.color != 0) return false;
        if (s1.subTreeScore == s2.subTreeScore)
            return s1.index < s2.index;
        return s1.subTreeScore > s2.subTreeScore;
    }
};
#endif //FRAGMENTATIONNODE_H
