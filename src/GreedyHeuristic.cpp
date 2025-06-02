//
// Created by Gregory Johnson on 5/15/25.
//

#include "FragmentationTree/GreedyHeuristic.h"

std::string GreedyHeuristic::CalculateHeuristic(FragmentationTree& tree) {
    tree.SortFragmentationNodes();
    const std::vector<FragmentationNode>& nodes = tree.GetFragmentationNodes();
    //Get all nodes of color 0
    size_t eligibleCandidateCounter = 0;
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (nodes[i].color != 0) continue;
        return nodes[i].formula.GetMolecularFormula();
    }
}
