//
// Created by gregj on 5/10/2025.
//

#include "DirectedAcyclicGraph/DirectedAcyclicGraph.h"
#include <Rcpp.h>
void DirectedAcyclicGraph::AddEdge(const size_t key, const size_t outGoingKey) {
    adjacencyList[key].push_back(outGoingKey);
}

std::list<size_t> DirectedAcyclicGraph::GetEdges(const size_t key) const {
    if (adjacencyList.find(key) == adjacencyList.end()) return {};
    return adjacencyList.at(key);
}

void DirectedAcyclicGraph::Print(const std::vector<FragmentationNode>& nodes) const {
    for (const auto& edge :  adjacencyList) {
        Rcpp:: Rcout << edge.first << ":" << nodes[edge.first].color << "-> ";;
        for (const auto& outGoingEdges : edge.second) {
            Rcpp::Rcout << "( "<<outGoingEdges << ": " << nodes[outGoingEdges].color << "), ";
        }
        Rcpp::Rcout << "\n";
    }
}

