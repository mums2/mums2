//
// Created by Gregory Johnson on 5/15/25.
//

#include "FragmentationTree/GreedyHeuristic.h"

std::string GreedyHeuristic::CalculateHeuristic(FragmentationTree& tree) {
    tree.SortFragmentationNodes();
    const std::vector<FragmentationNode>& nodes = tree.GetFragmentationNodes();
    //Get all nodes of color 0
    size_t eligibleCandidateCounter = 0;
    while (nodes[eligibleCandidateCounter].color == 0) {
        eligibleCandidateCounter++;
    }
    // If there is no subRoot, (which shouldn't be possible since there should not be circular sub-molecules)
    // then choose the highest score
    const FragmentationNode& candidate = nodes[0];
    return candidate.formula.GetMolecularFormula();
}
