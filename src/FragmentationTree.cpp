//
// Created by gregj on 5/10/2025.
//

#include "FragmentationTree/FragmentationTree.h"

#include "Chemicals/MolecularFormula/MolecularFormula.h"

void FragmentationTree::AddMolecularFormulasToGraph(const std::vector<std::string>& allIsotopeDecompositions,
                                                    const int color) {
    const size_t size = allIsotopeDecompositions.size();
    for (size_t i = 0; i < size; i++) {
        const std::string& molecularFormula = allIsotopeDecompositions[i];
        MolecularFormula formula(molecularFormula);
        size_t currentIndex = i + startingIndex;
        keyToColorMap[currentIndex] = color;
        keyToMolecularFormulaMap[currentIndex] = molecularFormula;
        for (size_t j = i + 1; j < size; j++) {
            const std::string& other = allIsotopeDecompositions[j];
            MolecularFormula currentFormula(allIsotopeDecompositions[j]);
            const size_t outGoingIndex = j + startingIndex;
            if (formula.CheckIfOtherIsSubFormula(currentFormula)) {
                graph.AddEdge(currentIndex, outGoingIndex);
                continue;
            }
            if (currentFormula.CheckIfOtherIsSubFormula(formula)) {
                graph.AddEdge(outGoingIndex, currentIndex);
            }

        }
    }
    startingIndex += size;
}
