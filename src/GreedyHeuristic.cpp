//
// Created by Gregory Johnson on 5/15/25.
//

#include "FragmentationTree/GreedyHeuristic.h"

std::string GreedyHeuristic::CalculateHeuristic(FragmentationTree& tree) {
    tree.SortFragmentationNodes();
    const std::vector<FragmentationNode>& nodes = tree.GetFragmentationNodes();
    // If there is no subRoot, (which shouldn't be possible since there should not be circular sub-molecules)
    // then choose the highest score
    size_t candidateIndex = 0;
    for (const auto& node : nodes) {
        candidateIndex++;
        if (node.color != 0) continue;
        candidateIndex--;
        break;
    }
    if (candidateIndex >= nodes.size()) {
        Rcpp::warning("GreedyHeuristic: node index out of bounds, returning first valid result. ");
        return "";
    }
    const FragmentationNode& candidate = nodes[candidateIndex];
    return candidate.formula.GetMolecularFormula();
}
