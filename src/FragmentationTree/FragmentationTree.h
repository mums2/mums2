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
#include "../DirectedAcyclicGraph/Vertex.h"


class FragmentationTree {
public:
    FragmentationTree() = default;
    void AddMolecularFormulasToGraph(const Rcpp::StringVector &molecularFormulas,
        const Rcpp::IntegerVector &color, const Rcpp::IntegerVector& decompositionScores,
        double parentMass, int amountOfUniqueColors);
    void SortVertexList();
    const std::vector<FragmentationNode> &GetFragmentationNodes() const {return molecularNodeList;}
    const std::vector<Vertex> &GetVertexList() const {return vertexList;}
    int GetUniqueColors() const {return uniqueColors;}


private:
    // Keys of the same color represent the same mz, int (isotope).
    int uniqueColors;
    std::vector<FragmentationNode> molecularNodeList;
    std::vector<Vertex> vertexList;
};



#endif //FRAGEMENTATIONTREE_H
