//
// Created by gregj on 5/10/2025.
//

#include "FragmentationTree/FragmentationTree.h"
#include "Chemicals/MolecularFormula/MolecularFormula.h"
#include "DiversityMetrics/DiversityMetricFactory.h"

FragmentationTree::FragmentationTree(const Rcpp::List& fragmentationData, const double parentMass):
     parentMass(parentMass) {
        Initialize(fragmentationData);
    }

void FragmentationTree::Initialize(const Rcpp::List& fragmentationData) {
    const Rcpp::StringVector& molecularFormulas = fragmentationData["formula"];
    const Rcpp::IntegerVector& color = fragmentationData["color"];
    const Rcpp::NumericVector& decompositionScores = fragmentationData["score"];
    const Rcpp::NumericVector& masses = fragmentationData["mass"];
    size = molecularFormulas.size();
    molecularNodeList = std::vector<FragmentationNode>(size);
    for (int i = 0; i < size; i++) {
        molecularNodeList[i] = FragmentationNode(color[i], i,
            decompositionScores[i], 0, MolecularFormula(molecularFormulas[i], masses[i]));
    }
}

void FragmentationTree::SortFragmentationNodes() {
    std::stable_sort(molecularNodeList.begin(), molecularNodeList.end(), CompareFragmentationNodes());
}

void FragmentationTree::CollectResultFromNode(const std::list<int>& parentIndexes, const double subtreeScore,
    const int index) {
    mutexLock.lock();
    for (const auto& parentIndex : parentIndexes) {
        molecularNodeList[parentIndex].subTreeScore +=  subtreeScore;
    }
    molecularNodeList[index].subTreeScore += subtreeScore;
    mutexLock.unlock();
}


void FragmentationTree::AddMolecularFormulaToGraph(const int currentIndex) {
    const std::vector<FragmentationNode>& fragmentationNodes = molecularNodeList;
    const MolecularFormula& formula = fragmentationNodes[currentIndex].formula;
    const FragmentationNode& fragmentationNode =  fragmentationNodes[currentIndex];
    std::list<int> parentIndexes;
    double finalSubtreeScore = 0;
    for (int j = currentIndex + 1; j < size; j++) {
        const MolecularFormula& currentFormula = fragmentationNodes[j].formula;
        // Nodes with similar fragmentation colors should never be a subformula
        if (fragmentationNodes[j].color == fragmentationNode.color) continue;
        const int res = formula.CheckIfOtherIsSubFormula(currentFormula);
        if (res == 0) continue;
        if (res == 2) { // meaning parentFormula is a subformula of the child
            parentIndexes.emplace_back(j);
            continue;
        }
        const double lossMass = formula.GetLossMass(currentFormula);
        const double score = std::log(std::abs(1 - lossMass/parentMass)) +
            fragmentationNodes[j].score + fragmentationNodes[j].subTreeScore;
        finalSubtreeScore += score;
    }
    CollectResultFromNode(parentIndexes, finalSubtreeScore, currentIndex);
}

