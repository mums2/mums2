//
// Created by gregj on 5/10/2025.
//

#include "FragmentationTree/FragmentationTree.h"



#include "Chemicals/MolecularFormula/MolecularFormula.h"
#include "DirectedAcyclicGraph/FragmentationNode.h"
#include "DiversityMetrics/DiversityMetricFactory.h"

void FragmentationTree::AddMolecularFormulasToGraph(const Rcpp::StringVector& molecularFormulas,
    const Rcpp::IntegerVector& color) {

    const size_t size = molecularFormulas.size();
    std::vector<FragmentationNode> fragmentationNodes(size);
    for (size_t i = 0; i < size; i++) {
        fragmentationNodes[i].formula = MolecularFormula(molecularFormulas[i]);
    }
    for (size_t i = 0; i < size; i++) {
        MolecularFormula& formula = fragmentationNodes[i].formula;
        fragmentationNodes[i].color = color[i];
        keyToMolecularFormulaMap[i] = fragmentationNodes[i];
        for (size_t j = i + 1; j < size; j++) {
            MolecularFormula& currentFormula = fragmentationNodes[j].formula;
            if (formula.CheckIfOtherIsSubFormula(currentFormula)) {
                graph.AddEdge(i, j);
                continue;
            }
            if (currentFormula.CheckIfOtherIsSubFormula(formula)) {
                graph.AddEdge(j, i);
            }

        }
   }
}

void FragmentationTree::PrintGraph() const {
    for (const auto& node : keyToMolecularFormulaMap) {
        const auto& edges = graph.GetEdges(node.first);
        Rcpp::Rcout << node.first << ":" << node.second.color << " -> ";
        for (const auto& edge : edges) {
            Rcpp::Rcout << "(" << edge  << ": " << keyToMolecularFormulaMap.at(edge).color << "), ";
        }
        Rcpp::Rcout << "\n";
    }
}
