//
// Created by Gregory Johnson on 5/15/25.
//

#include "FragmentationTree/GreedyHeuristic.h"

std::string GreedyHeuristic::CalculateHeuristic(FragmentationTree& tree) {
    const std::vector<FragmentationNode>& nodes = tree.GetFragmentationNodes();
    //Get all nodes of color 0
    size_t eligibleCandidateCounter = 0;
    while (nodes[eligibleCandidateCounter].color == 0) {
        eligibleCandidateCounter++;
    }
    std::vector<FragmentationNode> eligibleNodes(nodes.begin(), nodes.begin() + eligibleCandidateCounter);
    std::sort(eligibleNodes.begin(), eligibleNodes.end(), CompareFragmentationNodes());
    // If there is no subRoot, (which shouldn't be possible since there should not be circular sub-molecules)
    // then choose the highest score
    if (eligibleNodes.empty()) {
        Rcpp::warning("GreedyHeuristic: node index out of bounds, returning first valid result. ");
        return "";
    }
    const FragmentationNode& candidate = eligibleNodes[0];
    return candidate.formula.GetMolecularFormula();
}
