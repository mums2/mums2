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


class FragmentationTree {
public:
    FragmentationTree() = default;
    void AddMolecularFormulasToGraph(const Rcpp::StringVector &molecularFormulas, const Rcpp::IntegerVector &color);

private:
    // Keys of the same color represent the same mz, int (isotope).
    DirectedAcyclicGraph graph;
    std::unordered_map<size_t, std::string> keyToMolecularFormulaMap;
    std::unordered_map<size_t, size_t> keyToColorMap;
};



#endif //FRAGEMENTATIONTREE_H
