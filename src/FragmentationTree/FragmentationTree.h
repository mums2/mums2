//
// Created by gregj on 5/10/2025.
//

#ifndef FRAGEMENTATIONTREE_H
#define FRAGEMENTATIONTREE_H
#include <mutex>
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
    FragmentationTree(const Rcpp::List& fragmentationData, double);
    void Initialize(const Rcpp::List& fragmentationData);
    void AddMolecularFormulasToGraph(const Rcpp::StringVector &molecularFormulas,
        const Rcpp::IntegerVector &color, const Rcpp::NumericVector& decompositionScores,
        const Rcpp::NumericVector& masses, double parentMass);
    const std::vector<FragmentationNode> &GetFragmentationNodes() const {return molecularNodeList;}
    const std::vector<Vertex> &GetVertexList() const {return vertexList;}
    void SortFragmentationNodes();
    void AddMolecularFormulaToGraph(int currentIndex);
    void CollectResultFromNode(const std::list<int>& parentIndexes, double subtreeScore, int index);


private:
    // Keys of the same color represent the same mz, int (isotope).
    std::vector<FragmentationNode> molecularNodeList;
    std::vector<Vertex> vertexList;
    std::mutex mutexLock;
    double parentMass;
    int size;
};



#endif //FRAGEMENTATIONTREE_H
