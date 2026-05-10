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

#include "DetectNeutralLoses.h"
#include "FragmentationNode.h"
#include "../Decomposition/DecompResult.h"


class FragmentationTree {
public:
    FragmentationTree(const std::vector<DecompResult>&, double);
    int GetColorZeroCount() const {return colorZeroSize;}
    const std::vector<FragmentationNode> &GetFragmentationNodes() const {return molecularNodeList;}
    std::string GetBestFormula() const;
    void AddMolecularFormulaToGraph(int currentIndex);

private:
    void Initialize(const std::vector<DecompResult>& decompResults);
    void CollectResultFromNode(double subtreeScore, int index);
    // Keys of the same color represent the same mz, int (isotope).

    std::mutex mutexLock;
    double parentMass;
    int size;
    int colorZeroSize;
    std::vector<FragmentationNode> molecularNodeList;
    DetectNeutralLoses neutralLosesScorer;
};



#endif //FRAGEMENTATIONTREE_H
