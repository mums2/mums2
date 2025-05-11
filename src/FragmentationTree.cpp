//
// Created by gregj on 5/10/2025.
//

#include "FragmentationTree/FragmentationTree.h"



#include "Chemicals/MolecularFormula/MolecularFormula.h"
#include "DiversityMetrics/DiversityMetricFactory.h"

void FragmentationTree::AddMolecularFormulasToGraph(const Rcpp::StringVector& molecularFormulas,
    const Rcpp::IntegerVector& color) {

    const size_t size = molecularFormulas.size();
    std::vector<MolecularFormula> formulas(size);
    for (size_t i = 0; i < size; i++) {
        formulas[i] = MolecularFormula(molecularFormulas[i]);
    }
    for (size_t i = 0; i < size; i++) {
        const Rcpp::String& molecularFormula = molecularFormulas[i];
        keyToColorMap[i] = color[i];
        keyToMolecularFormulaMap[i] = molecularFormula;
        MolecularFormula& formula = formulas[i];
        for (size_t j = i + 1; j < size; j++) {
            MolecularFormula& currentFormula = formulas[j];
            // const Rcpp::String other = molecularFormulas[j];
            if (formula.CheckIfOtherIsSubFormula(currentFormula)) {
                graph.AddEdge(i, j);
                //Rcpp::Rcout << molecularFormula.get_cstring() << " -> " << other.get_cstring() << std::endl;
                continue;
            }
            if (currentFormula.CheckIfOtherIsSubFormula(formula)) {
                graph.AddEdge(j, i);
                //Rcpp::Rcout << other.get_cstring() << " -> " << molecularFormula.get_cstring() << std::endl;
            }

        }
   }
}
