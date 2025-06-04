//
// Created by Gregory Johnson on 5/15/25.
//

#include "FragmentationTree/GreedyHeuristic.h"

std::string GreedyHeuristic::CalculateHeuristic(FragmentationTree& tree) {
    // In worse case scenario, return the first node, which will always
    // be a node of the parent index
    tree.SortFragmentationNodes();
    return tree.GetFragmentationNodes()[0].formula.GetMolecularFormula();
}
