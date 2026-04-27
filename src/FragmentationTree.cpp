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

FragmentationTree::FragmentationTree(const std::vector<FragmentationNode>& fragmentationData,
    const double parentMass):
     parentMass(parentMass) {
    Initialize(fragmentationData);
}

void FragmentationTree::Initialize(const Rcpp::List& fragmentationData) {
    const Rcpp::StringVector& molecularFormulas = fragmentationData["formula"];
    const Rcpp::IntegerVector& color = fragmentationData["color"];
    const Rcpp::NumericVector& decompositionScores = fragmentationData["score"];
    const Rcpp::NumericVector& masses = fragmentationData["mass"];
    size = molecularFormulas.size();
    colorZeroFormulas.reserve(size);
    molecularNodeList = std::vector<FragmentationNode>(size);
    for (int i = 0; i < size; i++) {
        molecularNodeList[i] = FragmentationNode(color[i], i,
            decompositionScores[i], 0, MolecularFormula(molecularFormulas[i], masses[i]));
        if (color[i] == 0) {
            colorZeroFormulas.emplace_back(i);
        }
    }
}

void FragmentationTree::Initialize(const std::vector<FragmentationNode>& fragmentationData) {
    size = static_cast<int>(fragmentationData.size());
    colorZeroFormulas.reserve(size);
    molecularNodeList = fragmentationData;
    for (int i = 0; i < size; i++) {
        if (fragmentationData[i].color == 0) {
            colorZeroFormulas.emplace_back(i);
        }
    }
}

void FragmentationTree::SortFragmentationNodes() {
    std::stable_sort(molecularNodeList.begin(), molecularNodeList.end(), CompareFragmentationNodes());
}

void FragmentationTree::CollectResultFromNode(const double subtreeScore,
    const int index) {
   // mutexLock.lock();
    // for (const auto& parentIndex : parentIndexes) {
    //     molecularNodeList[parentIndex].subTreeScore +=  subtreeScore;
    // }
    molecularNodeList[index].subTreeScore += subtreeScore;
 //   mutexLock.unlock();
}


void FragmentationTree::AddMolecularFormulaToGraph(const int currentIndex) {
    const std::vector<FragmentationNode>& fragmentationNodes = molecularNodeList;
    const MolecularFormula& formula = fragmentationNodes[currentIndex].formula;
    const FragmentationNode& fragmentationNode = fragmentationNodes[currentIndex];
    double finalSubtreeScore = 0;
    for (int j = 0; j < size; j++) {
        if (j == currentIndex) continue;
        const MolecularFormula& currentFormula = fragmentationNodes[j].formula;
        // Nodes with similar fragmentation colors should never be a subformula
        if (fragmentationNodes[j].color == fragmentationNode.color) continue;
        const bool res = formula.CheckIfOtherIsSubFormula(currentFormula);
        if (!res) continue;
        const double lossMass = formula.GetLossMass(currentFormula);
        const double score = std::log(std::abs(1 - lossMass/parentMass)) +
            fragmentationNodes[j].score;
        finalSubtreeScore += score;
    }
    CollectResultFromNode(finalSubtreeScore, currentIndex);
}

