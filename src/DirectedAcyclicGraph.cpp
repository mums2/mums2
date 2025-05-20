//
// Created by gregj on 5/10/2025.
//

#include "DirectedAcyclicGraph/DirectedAcyclicGraph.h"
#include <Rcpp.h>
#include <fstream>
void DirectedAcyclicGraph::AddEdge(const size_t key, const size_t outGoingKey) {
    adjacencyList[key].push_back(outGoingKey);
    nodeParentCount[outGoingKey]++;
}

std::list<size_t> DirectedAcyclicGraph::GetEdges(const size_t key) const {
    if (adjacencyList.find(key) == adjacencyList.end()) return {};
    return adjacencyList.at(key);
}

void DirectedAcyclicGraph::Print(const std::vector<FragmentationNode>& nodes) const {
    std::string output;
    for (const auto& edge :  adjacencyList) {
        for (const auto& outGoingEdges : edge.second) {
            output += ("'" + nodes[edge.first].formula.GetMolecularFormula() + ":" +
                std::to_string(nodes[edge.first].color) + "' -> ");
            output += ("'" + nodes[outGoingEdges].formula.GetMolecularFormula() + ":" +
                std::to_string(nodes[outGoingEdges].color) + "'\n");
        }
        //Rcpp::Rcout << "\n";
    }
    std::ofstream outputFile;
    outputFile.open("F:\\mums2\\output.txt");
    outputFile << output;
    outputFile.close();
}

std::list<size_t> DirectedAcyclicGraph::FindRoots() const {
    std::list<size_t> roots;
    for (const auto& edge : adjacencyList) {
        if (nodeParentCount.find(edge.first) != nodeParentCount.end()) continue;
        // This is a root, it has no parent node
        roots.emplace_back(edge.first);
    }
    return roots;
}

