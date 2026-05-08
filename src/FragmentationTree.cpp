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

    colorZeroSize = 0;
    size = 0;
    ranges = std::vector<int>(decompResults.size() + 1, 0);
    int count = 1;
    for (const auto& decompResult : decompResults) {
        size += static_cast<int>(decompResult.formula.size());
        ranges[count] = ranges[count - 1] + static_cast<int>(decompResult.formula.size());
        count++;
    }
    molecularNodeList = std::vector<FragmentationNode>(size);
    int indexPosition = 0;
    int color = 0;
    for (const auto & fragData : decompResults) {
        const int currentSize = static_cast<int>(fragData.formula.size());
        const int indexToGetTo = (currentSize + indexPosition) - indexPosition;
        for (int i = 0; i < indexToGetTo; i++) {
            molecularNodeList[i + indexPosition] = FragmentationNode(color, fragData.score[i],
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
    size_t currentBestIndex = 0;
    size_t topColorAmount = molecularNodeList[0].amountOfDistinctColors;
    double score = molecularNodeList[0].score;

    for (size_t i = 1; i < colorZeroSize; i++) {
        const FragmentationNode& node = molecularNodeList[i];
        if (node.amountOfDistinctColors > topColorAmount) {
            topColorAmount = node.amountOfDistinctColors;
            score = node.score;
            currentBestIndex = i;
            continue;
        }
        if (node.amountOfDistinctColors == topColorAmount && node.score > score) {
            score = node.score;
            currentBestIndex = i;
        }
    }
    const FragmentationNode& node = molecularNodeList[currentBestIndex];
    return node.formula.GetMolecularFormula();
    // const auto it =  std::max_element(molecularNodeList.cbegin(), molecularNodeList.cbegin() + colorZeroSize,
    //     CompareFragmentationNodes());
    // return it->formula.GetMolecularFormula();
}

void FragmentationTree::CollectResultFromNode(const double subtreeScore,
    const int index) {
   // mutexLock.lock();
    // for (const auto& parentIndex : parentIndexes) {
    //     molecularNodeList[parentIndex].subTreeScore +=  subtreeScore;
    // }
    molecularNodeList[index].score += subtreeScore;
    // mutexLock.unlock();
}

// 100 Range: 0, 100
// 20  Range: 100, 120
// 10  Range: 120, 130
// 5   Range: 130, 135
// TotalSize = 135
// TotalSize - (Summation(CurrentSize) from 1 to currentN
// Starting Index
// currentN = 1 -> TotalSize - CurrentSize(N1) -> 135 - 5 = 130
// currentN = 2 -> TotalSize - (CurrentSize(N1) + CurrentSize(N2)) -> 135 - (5 + 10) = 120
// currentN = 1 -> TotalSize) - (CurrentSize(N1) + CurrentSize(N2)) + CurrentSize(N3)) -> 135 - (5 + 10 + 20) = 100
// currentN = 1 -> TotalSize - (CurrentSize(N1) + CurrentSize(N2)) + CurrentSize(N3) + CurrentSize(N4))-> 135 - (5 + 10 + 20 + 100) = 0
// Ending Index
// (TotalSize - ((Summation(CurrentSize)) from 1 to currentN - 1)
// end1 = 1 -> TotalSize - CurrentSize(N0) -> 135 - 0 = 135
// end2 = 2 -> TotalSize - (CurrentSize(N1) -> 135 - (5) = 130
// end3 = 1 -> TotalSize) - (CurrentSize(N1) + CurrentSize(N2)) -> 135 - (5 + 10) = 120
// end4 = 1 -> TotalSize - (CurrentSize(N1) + CurrentSize(N2)) + CurrentSize(N3) -> 135 - (5 + 10 + 20) = 100
//
// for(size_t i = currentN; i < end1; i++)
// at 1 -> i = 130, end1 = 135
// at 2 -> i = 120, end1 = 130
// at 3 -> i = 100, end1 = 120
// at 4 -> i = 0, end1 = 100
// arr(0, 100, 120, 130, 135)
// arr[4 - 1] arr[4 - 0]
// arr[4 - 2] arr[4 - 1]
// arr[4 - 3] arr[4 - 2]
// arr[4 - 4] arr[4 - 3]
// May have to store the number of colors that explain the formulas
// And choose the formula that is most explained by the colors.
void FragmentationTree::AddMolecularFormulaToGraph(const int currentIndex,
    const size_t startingIndex, const size_t endingIndex, const DetectNeutralLoses& neutralLosesScorer) {
    const std::vector<FragmentationNode>& fragmentationNodes = molecularNodeList;
    const MolecularFormula& formula = fragmentationNodes[currentIndex].formula;
    const double fragmentHeteroCarbonRatio = formula.GetHeteroToCarbonRatio();
    double highestScore = 0;
    size_t chosenIndex = 0;
    bool hasScored = false;
    for (size_t j = startingIndex; j < endingIndex; j++) {
        const FragmentationNode& currentFragmentation = fragmentationNodes[j];
        const MolecularFormula& currentFormula = currentFragmentation.formula;
        // If it is a higher color, the mass will be higher than or equal to.
        // if (currentFragmentation.color > fragmentationNodes[currentIndex].color) continue;
        // Nodes with similar fragmentation colors should never be a subformula
        const bool res = formula.CheckIfOtherIsSubFormula(currentFormula);
        if (!res) continue;
        const double currentHeteroCarbonRatio = currentFormula.GetHeteroToCarbonRatio();
        double heteroCarbonRatioScore = 0;
        if (fragmentHeteroCarbonRatio < currentHeteroCarbonRatio) {
            heteroCarbonRatioScore = fragmentHeteroCarbonRatio - currentHeteroCarbonRatio;
        }
        const double neutralLoseScore = neutralLosesScorer.DetermineNeutralLoses(
            formula - currentFormula);

        const double lossMass = formula.GetLossMass(currentFormula);
        const double score = std::log(std::abs(1 - lossMass/parentMass)) +
            currentFragmentation.score + neutralLoseScore + heteroCarbonRatioScore;
        if (!hasScored || highestScore < score) {
            hasScored = true;
            highestScore = score;
            chosenIndex = j;
        }
    }
    mutexLock.lock();
    molecularNodeList[currentIndex].amountOfDistinctColors = molecularNodeList[chosenIndex].amountOfDistinctColors + 1;
    molecularNodeList[currentIndex].score += highestScore;
    // CollectResultFromNode(highestScore, currentIndex);
    mutexLock.unlock();
}
