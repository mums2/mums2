//
// Created by Gregory Johnson on 5/15/25.
//

#include "FragmentationTree/GreedyHeuristic.h"

std::string GreedyHeuristic::CalculateHeuristic(FragmentationTree& tree) {
    // In worse case scenario, return the first node, which will always
    // be a node of the parent index
    tree.SortFragmentationNodes();
    const std::vector<FragmentationNode>& nodes = tree.GetFragmentationNodes();
    FragmentationNode result;
    result.index = -1;
    //Get all nodes of color 0
    for (const auto & node : nodes) {
        if (node.color != 0) continue;
        if (result.index == -1) {
            result = node;
            continue;
        }
        if (result.subTreeScore == node.subTreeScore) {
            if (result.index > node.index)
                return node.formula.GetMolecularFormula();
        }
        return result.formula.GetMolecularFormula();
    }
    // Should be unreachable
    Rcpp::warning("Returning empty molecular formula. ");
    return "";
}
