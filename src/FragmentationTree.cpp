//
// Created by gregj on 5/10/2025.
//

#include "FragmentationTree/FragmentationTree.h"
#include "Chemicals/MolecularFormula/MolecularFormula.h"
#include "Decomposition/DecomposeMass.h"
#include "DiversityMetrics/DiversityMetricFactory.h"

FragmentationTree::FragmentationTree(const std::vector<DecompResult>& decompResults,
    const double parentMass):
     parentMass(parentMass) {
    Initialize(decompResults);
}

void FragmentationTree::Initialize(const std::vector<DecompResult>& decompResults) {


    size = 0;
    colorZeroSize = 0;
    for (const auto& decompResult : decompResults) {
        size += static_cast<int>(decompResult.formula.size());
    }
    molecularNodeList = std::vector<FragmentationNode>(size);
    int indexPosition = 0;
    int color = 0;
    for (const auto & fragData : decompResults) {
        const int currentSize = static_cast<int>(fragData.formula.size());
        const int indexToGetTo = (currentSize + indexPosition) - indexPosition;
        for (int i = 0; i < indexToGetTo; i++) {
            molecularNodeList[i + indexPosition] = FragmentationNode(color,
                       fragData.score[i], fragData.score[i],
                       MolecularFormula(fragData.formula[i], fragData.exactmass[i]));
            if (color == 0) {
                colorZeroSize++;
            }
        }
        indexPosition += currentSize;
        color++;
    }

}

std::string FragmentationTree::GetBestFormula() const {
    return std::max_element(molecularNodeList.cbegin(), molecularNodeList.cbegin() + colorZeroSize,
        CompareFragmentationNodes())->formula.GetMolecularFormula();
}

void FragmentationTree::CollectResultFromNode(const double subtreeScore,
    const int index) {
   mutexLock.lock();
    // for (const auto& parentIndex : parentIndexes) {
    //     molecularNodeList[parentIndex].subTreeScore +=  subtreeScore;
    // }
    molecularNodeList[index].subTreeScore += subtreeScore;
    mutexLock.unlock();
}


void FragmentationTree::AddMolecularFormulaToGraph(const int currentIndex) {
    const std::vector<FragmentationNode>& fragmentationNodes = molecularNodeList;
    const MolecularFormula& formula = fragmentationNodes[currentIndex].formula;
    double finalSubtreeScore = 0;
    for (int j = colorZeroSize; j < size; j++) {
        if (j == currentIndex) continue;
        const MolecularFormula& currentFormula = fragmentationNodes[j].formula;
        // Nodes with similar fragmentation colors should never be a subformula
        const bool res = formula.CheckIfOtherIsSubFormula(currentFormula);
        if (!res) continue;
        const double lossMass = formula.GetLossMass(currentFormula);
        const double score = std::log(std::abs(1 - lossMass/parentMass)) +
            fragmentationNodes[j].score;
        finalSubtreeScore += score;
    }
    CollectResultFromNode(finalSubtreeScore, currentIndex);
}

