//
// Created by gregj on 5/10/2025.
//

#include "FragmentationTree/FragmentationTree.h"



#include "Chemicals/MolecularFormula/MolecularFormula.h"
#include "DirectedAcyclicGraph/FragmentationNode.h"
#include "DiversityMetrics/DiversityMetricFactory.h"

void FragmentationTree::AddMolecularFormulasToGraph(const Rcpp::StringVector &molecularFormulas,
        const Rcpp::IntegerVector &color, const Rcpp::NumericVector& decompositionScores,
        const double parentMass, const int amountOfUniqueColors) {
    uniqueColors = amountOfUniqueColors;
    const size_t size = molecularFormulas.size();
    std::vector<FragmentationNode> fragmentationNodes(size);
    for (size_t i = 0; i < size; i++) {
        fragmentationNodes[i] = FragmentationNode{color[i], i,
            decompositionScores[i], 0, MolecularFormula(molecularFormulas[i])};
    }
    for (size_t i = 0; i < size; i++) {
        MolecularFormula& formula = fragmentationNodes[i].formula;
        for (size_t j = i + 1; j < size; j++) {
            MolecularFormula& currentFormula = fragmentationNodes[j].formula;
            // Nodes with similar fragmentation colors should never be a subformula
            if (fragmentationNodes[j].color == fragmentationNodes[i].color) continue;
            int res = formula.CheckIfOtherIsSubFormula(currentFormula);
            if (res == 0) continue;
            if (res == 2) { // meaning parentFormula is a subformula of the child
                fragmentationNodes[i].parentIndexes.emplace_back(j);
                continue;
            }
             const double lossMass = formula.GetLossMass(currentFormula);
             const double score = std::log(std::abs(1 - lossMass/parentMass)) +
                 fragmentationNodes[j].score + fragmentationNodes[j].subTreeScore;
            fragmentationNodes[i].subTreeScore += score;
        }
        for (const auto& parentIndex:  fragmentationNodes[i].parentIndexes) {
            fragmentationNodes[parentIndex].subTreeScore += (fragmentationNodes[i].subTreeScore +
                fragmentationNodes[parentIndex].score);
        }
   }
    molecularNodeList = fragmentationNodes;
}

void FragmentationTree::SortFragmentationNodes() {
    std::sort(molecularNodeList.begin(), molecularNodeList.end(), CompareFragmentationNodes());
}

