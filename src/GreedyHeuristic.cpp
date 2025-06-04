//
// Created by Gregory Johnson on 5/15/25.
//

#include "FragmentationTree/GreedyHeuristic.h"

std::string GreedyHeuristic::CalculateHeuristic(FragmentationTree& tree) {
    // In worse case scenario, return the first node, which will always
    // be a node of the parent index
    tree.SortFragmentationNodes();
    // const std::vector<FragmentationNode>& nodes = tree.GetFragmentationNodes();
    //Get all nodes of color 0
    // for (const auto & node : nodes) {
    //     if (node.color != 0) continue;
    //     return node.formula.GetMolecularFormula();
    // }
    // Should be unreachable
    return tree.GetFragmentationNodes()[0].formula.GetMolecularFormula();
    // Rcpp::warning("Returning empty molecular formula. ");
    // return "";
}
