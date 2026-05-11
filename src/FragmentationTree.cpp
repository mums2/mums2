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

    uniqueColors = decompResults.size();
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
                       fragData.score[i], MolecularFormula(fragData.formula[i],
                           fragData.exactmass[i]));
            if (color == 0) {
                colorZeroSize++;
            }
        }
        indexPosition += currentSize;
        color++;
    }

}

std::string FragmentationTree::GetBestFormula() const {
    int currentBestIndex = 0;
    int topColorAmount = molecularNodeList[0].amountOfDistinctColors;
    double score = molecularNodeList[0].score;

    for (int i = 1; i < colorZeroSize; i++) {
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

// May have to store the number of colors that explain the formulas
// And choose the formula that is most explained by the colors.
void FragmentationTree::AddMolecularFormulaToGraph(const int currentIndex) {
    const std::vector<FragmentationNode>& fragmentationNodes = molecularNodeList;
    const MolecularFormula& formula = fragmentationNodes[currentIndex].formula;
    const double fragmentHeteroCarbonRatio = formula.GetHeteroToCarbonRatio();

    std::vector<bool> hasScoredColors(uniqueColors, false);
    std::vector<double> highestColorScores(uniqueColors, 0);
    int distinctColors = 0;
    for (int j = colorZeroSize; j < size; j++) {
        if (j == currentIndex) continue;
        const FragmentationNode& currentFragmentation = fragmentationNodes[j];
        const MolecularFormula& currentFormula = currentFragmentation.formula;
        // Nodes with similar fragmentation colors should never be a subformula
        const bool res = formula.CheckIfOtherIsSubFormula(currentFormula);
        if (!res) continue;
        const double currentHeteroCarbonRatio = currentFormula.GetHeteroToCarbonRatio();
        double heteroCarbonRatioScore = 0;
        if (fragmentHeteroCarbonRatio > currentHeteroCarbonRatio) {
            heteroCarbonRatioScore = currentHeteroCarbonRatio - fragmentHeteroCarbonRatio;
        }
        const double neutralLoseScore = neutralLosesScorer.DetermineNeutralLoses(
            formula - currentFormula);

        const double lossMass = formula.GetLossMass(currentFormula);
        const double score = std::log(std::abs(1 - lossMass/parentMass)) +
            currentFragmentation.score + neutralLoseScore + heteroCarbonRatioScore;
        if (!hasScoredColors[currentFragmentation.color] || score > highestColorScores[currentFragmentation.color]) {
            highestColorScores[currentFragmentation.color] = score;
            hasScoredColors[currentFragmentation.color] = true;
        }
    }
    for (const auto& color : hasScoredColors) {
        if (color)
            distinctColors++;
    }
    // Score
    mutexLock.lock();
    molecularNodeList[currentIndex].amountOfDistinctColors = distinctColors;
    molecularNodeList[currentIndex].score += std::accumulate(highestColorScores.begin(),
        highestColorScores.end(), 0.0);
    mutexLock.unlock();
}
