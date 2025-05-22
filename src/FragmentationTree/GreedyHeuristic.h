//
// Created by Gregory Johnson on 5/15/25.
//

#ifndef GREEDYHEURISTIC_H
#define GREEDYHEURISTIC_H
#include "FragmentationTree.h"


class GreedyHeuristic {
public:
    std::string CalculateHeuristic(FragmentationTree& tree);
    void Print(const std::list<Vertex>&, const std::vector<FragmentationNode>&) const;
};



#endif //GREEDYHEURISTIC_H
