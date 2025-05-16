//
// Created by Gregory Johnson on 5/15/25.
//

#include "FragmentationTree/GreedyHeuristic.h"

void GreedyHeuristic::CalculateHeuristic(FragmentationTree& tree) {
    tree.SortVertexList();
    const int colors = tree.GetUniqueColors();
    std::vector<int> usedColors(colors, -1);
    const std::vector<FragmentationNode>& nodes = tree.GetFragmentationNodes();
    const std::vector<Vertex>& vertexes = tree.GetVertexList();
    const DirectedAcyclicGraph& graph = tree.GetGraph();
    const size_t numVertices = vertexes.size();
    std::vector<std::vector<bool>> hasVisited(colors,
        std::vector<bool>(numVertices, false));
    std::list<Vertex> visited;
    double score = 0;
    for (const auto& vertex : vertexes) {
        // if (std::isnan(vertex.score)) continue;
        const FragmentationNode& node = nodes[vertex.indexChildNode];
        if (hasVisited[node.color][node.index]) continue;
        hasVisited[node.color][node.index] = true;
        visited.emplace_back(vertex);
        score += vertex.score;
        // Rcpp::Rcout << vertex.score << std::endl;
    }
    Rcpp::Rcout << "Score: " << score;
    Print(visited, nodes);
}

void GreedyHeuristic::Print(const std::list<Vertex>& subtree, const std::vector<FragmentationNode>& nodes) const {
    DirectedAcyclicGraph graph;
    for (const auto& vertex : subtree) {
        graph.AddEdge(vertex.indexParentNode, vertex.indexChildNode);
    }
    graph.Print(nodes);

}

/*
 *
 * used[1,0] = true
 * used[3,0] = true
 * used[2,0] = true
 */
