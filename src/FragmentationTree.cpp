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
    vertexList.reserve(size * size);
    for (size_t i = 0; i < size; i++) {
        fragmentationNodes[i] = FragmentationNode{true, color[i], i,
            decompositionScores[i], 0, MolecularFormula(molecularFormulas[i])};
        fragmentationNodes[i].self = fragmentationNodes[i].formula.GetMolecularFormula();
    }
    for (size_t i = 0; i < size; i++) {
        MolecularFormula& formula = fragmentationNodes[i].formula;
        for (size_t j = i + 1; j < size; j++) {

            MolecularFormula& currentFormula = fragmentationNodes[j].formula;
            if (formula.CheckIfOtherIsSubFormula(currentFormula)) {
                // We only want the subtrees that are roots and have the highest score
                // const MolecularFormula loss = MolecularFormula(formula - currentFormula);
                const double lossMass = formula.GetLossMass(currentFormula);
                const double score = std::log(std::abs(1 - lossMass/parentMass)) +
                    fragmentationNodes[j].score;
                fragmentationNodes[i].subTreeScore += (score + fragmentationNodes[j].subTreeScore);
                fragmentationNodes[j].isSubtreeRoot = false;
                fragmentationNodes[j].parentFormula.emplace_back(formula.GetMolecularFormula());
                continue;
            }
            if (currentFormula.CheckIfOtherIsSubFormula(formula)) {
                // const MolecularFormula loss = MolecularFormula(formula - currentFormula);
                const double lossMass = formula.GetLossMass(currentFormula);
                const double score = std::log(std::abs(1 - lossMass/parentMass)) +
                    fragmentationNodes[i].score;
                fragmentationNodes[j].subTreeScore += (score + fragmentationNodes[i].subTreeScore);
                fragmentationNodes[i].isSubtreeRoot = false;
                fragmentationNodes[i].parentFormula.emplace_back(currentFormula.GetMolecularFormula());
            }
        }
   }
    molecularNodeList = fragmentationNodes;
}

void FragmentationTree::SortFragmentationNodes() {
    std::sort(molecularNodeList.begin(), molecularNodeList.end(), CompareFragmentationNodes());
}

