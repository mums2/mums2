//
// Created by gregj on 5/14/2025.
//

#ifndef FRAGMENTATIONNODE_H
#define FRAGMENTATIONNODE_H
#include "../Chemicals/MolecularFormula/MolecularFormula.h"

struct FragmentationNode {
    int color{};
    size_t index{};
    MolecularFormula formula;
};

#endif //FRAGMENTATIONNODE_H
