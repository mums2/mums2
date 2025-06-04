//
// Created by Gregory Johnson on 5/15/25.
//

#include "FragmentationTree/GreedyHeuristic.h"

std::string GreedyHeuristic::CalculateHeuristic(FragmentationTree& tree) {
    // In worse case scenario, return the first node, which will always
    // be a node of the parent index
    tree.SortFragmentationNodes();

    const std::vector<FragmentationNode>& nodes = tree.GetFragmentationNodes();
    const double maxScore = nodes[0].subTreeScore;
    size_t chosenIndex = 0;
    for (size_t i = 1; i < nodes.size(); ++i) {
        if (nodes[i].subTreeScore < maxScore) {
            break;
        }
        if (nodes[i].index < nodes[chosenIndex].index) chosenIndex = i;
    }
    Rcpp::Rcout << "Chosen Index: " << chosenIndex << std::endl;
    //Get all nodes of color 0
    // for (const auto & node : nodes) {
    //     if (node.color != 0) continue;
    //     return node.formula.GetMolecularFormula();
    // }
    // Should be unreachable
    return tree.GetFragmentationNodes()[chosenIndex].formula.GetMolecularFormula();
    // Rcpp::warning("Returning empty molecular formula. ");
    // return "";
}
