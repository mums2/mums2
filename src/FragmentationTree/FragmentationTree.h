//
// Created by gregj on 5/10/2025.
//

#ifndef FRAGEMENTATIONTREE_H
#define FRAGEMENTATIONTREE_H
#include <string>
#include <unordered_map>
#include <vector>
#include <Rcpp.h>

#include "../DirectedAcyclicGraph/DirectedAcyclicGraph.h"
#include "../DirectedAcyclicGraph/FragmentationNode.h"


class FragmentationTree {
public:
    FragmentationTree() = default;
    void AddMolecularFormulasToGraph(const Rcpp::StringVector &molecularFormulas, const Rcpp::IntegerVector &color);
    void PrintGraph() const;

private:
    // Keys of the same color represent the same mz, int (isotope).
    DirectedAcyclicGraph graph;
    std::unordered_map<size_t, FragmentationNode> keyToMolecularFormulaMap;
};



#endif //FRAGEMENTATIONTREE_H
