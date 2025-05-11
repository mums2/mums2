//
// Created by gregj on 5/10/2025.
//

#ifndef DIRECTEDACYCLICGRAPH_H
#define DIRECTEDACYCLICGRAPH_H
#include <list>
#include <unordered_map>



class DirectedAcyclicGraph {
public:
    DirectedAcyclicGraph() = default;
    void AddEdge(size_t key, size_t outGoingKey);
    std::list<size_t> GetEdges(size_t key);
private:
    std::unordered_map<size_t,std::list<size_t>> adjacencyList{};
};



#endif //DIRECTEDACYCLICGRAPH_H
