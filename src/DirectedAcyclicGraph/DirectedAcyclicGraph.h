//
// Created by gregj on 5/10/2025.
//

#ifndef DIRECTEDACYCLICGRAPH_H
#define DIRECTEDACYCLICGRAPH_H
#include <list>
#include <unordered_map>
#include <unordered_set>


class DirectedAcyclicGraph {
public:
    DirectedAcyclicGraph() = default;
    void AddEdge(size_t key, size_t outGoingKey);
    [[nodiscard]] std::list<size_t> GetEdges(size_t key) const;
private:
    std::vector<std::unordered_set<size_t>> adjacencyListOfParents{};
    std::vector<std::unordered_set<size_t>> adjacencyListOfChildren{};
    std::unordered_map<size_t,std::list<size_t>> adjacencyList{};
};



#endif //DIRECTEDACYCLICGRAPH_H
