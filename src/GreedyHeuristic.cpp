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
    // Rcpp::Rcout << "Score: " << candidate.subTreeScore << std::endl;
    // Rcpp::Rcout << "Formula: " << candidate.formula.GetMolecularFormula() << std::endl;
    // Rcpp::Rcout << "Color: " << candidate.color << std::endl;

    // int maxUniqueColors = 0;
    // size_t currentMaxIndex = 0;
    //
    // for (const auto& root : roots) {
    //     std::vector<size_t> uniqueColors(colors, -1);
    //     int amountOfUniqueColors = 0;
    //     for (const auto& children : graph.GetEdges(root)) {
    //         if (uniqueColors[nodes[children].color] != -1) continue;
    //         uniqueColors[nodes[children].color] = 0;
    //         amountOfUniqueColors++;
    //     }
    //     if (amountOfUniqueColors < maxUniqueColors) continue;
    //     currentMaxIndex = root;
    //     maxUniqueColors = amountOfUniqueColors;
    // }

    //Rcpp::Rcout << "Max unique colors: " << maxUniqueColors << std::endl;
    //Rcpp::Rcout << "Best Formula: " << nodes[currentMaxIndex].formula.GetMolecularFormula() << " Mass: " << nodes[currentMaxIndex].formula.GetMass();

    // Print(visited, nodes);

}

void GreedyHeuristic::Print(const std::list<Vertex>& subtree, const std::vector<FragmentationNode>& nodes) const {
    DirectedAcyclicGraph graph;
    for (const auto& vertex : subtree) {
        graph.AddEdge(vertex.indexParentNode, vertex.indexChildNode);
    }
    graph.Print(nodes);
    const std::list<size_t> roots = graph.FindRoots();
    for (const auto& root : roots) {
        Rcpp::Rcout << "Root: " << root << std::endl;
    }
}

/*
 *
 * used[1,0] = true
 * used[3,0] = true
 * used[2,0] = true
 */
