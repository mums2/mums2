//
// Created by gregj on 5/10/2025.
//

#ifndef DIRECTEDACYCLICGRAPH_H
#define DIRECTEDACYCLICGRAPH_H
#include <list>
#include <unordered_map>
#include <unordered_set>

#include "FragmentationNode.h"


class DirectedAcyclicGraph {
public:
    DirectedAcyclicGraph() = default;
    void AddEdge(size_t key, size_t outGoingKey);
    [[nodiscard]] std::list<size_t> GetEdges(size_t key) const;
    void Print(const std::vector<FragmentationNode>&) const;
    std::list<size_t> FindRoots() const;
private:
    std::unordered_map<size_t, size_t> nodeParentCount{};
    std::unordered_map<size_t,std::list<size_t>> adjacencyList{};
};


#endif //DIRECTEDACYCLICGRAPH_H
