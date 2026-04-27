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
#include "FragmentationNode.h"


class FragmentationTree {
public:
    FragmentationTree() = default;
    FragmentationTree(const Rcpp::List& fragmentationData, double);
    FragmentationTree(const std::vector<FragmentationNode>&, double);
    const std::vector<FragmentationNode> &GetFragmentationNodes() const {return molecularNodeList;}
    void SortFragmentationNodes();
    void AddMolecularFormulaToGraph(int currentIndex);
    const std::vector<int>& GetColorZeroFormulas() {return colorZeroFormulas;}

private:
    void Initialize(const Rcpp::List& fragmentationData);
    void Initialize(const std::vector<FragmentationNode>& fragmentationData);
    void CollectResultFromNode(double subtreeScore, int index);
    // Keys of the same color represent the same mz, int (isotope).
    std::vector<FragmentationNode> molecularNodeList;
    std::vector<int> colorZeroFormulas;
    std::mutex mutexLock;
    double parentMass;
    int size;
};



#endif //FRAGEMENTATIONTREE_H
