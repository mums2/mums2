//
// Created by gregj on 5/10/2025.
//

#include "FragmentationTree/FragmentationTree.h"



#include "Chemicals/MolecularFormula/MolecularFormula.h"
#include "DirectedAcyclicGraph/FragmentationNode.h"
#include "DiversityMetrics/DiversityMetricFactory.h"

void FragmentationTree::AddMolecularFormulasToGraph(const Rcpp::StringVector &molecularFormulas,
        const Rcpp::IntegerVector &color, const Rcpp::IntegerVector& decompositionScores,
        const double parentMass, const int amountOfUniqueColors) {
    uniqueColors = amountOfUniqueColors;
    const size_t size = molecularFormulas.size();
    std::vector<FragmentationNode> fragmentationNodes(size);
    vertexList.reserve(size * size);
    for (size_t i = 0; i < size; i++) {
        fragmentationNodes[i] = FragmentationNode{color[i], i,
            MolecularFormula(molecularFormulas[i])};
    }
    for (size_t i = 0; i < size; i++) {
        MolecularFormula& formula = fragmentationNodes[i].formula;
        for (size_t j = i + 1; j < size; j++) {
            MolecularFormula& currentFormula = fragmentationNodes[j].formula;
            if (formula.CheckIfOtherIsSubFormula(currentFormula)) {
                const MolecularFormula loss = MolecularFormula(formula - currentFormula);
                const double score = (1 - loss.GetMass()/parentMass) + decompositionScores[j];
                graph.AddEdge(i, j);
                vertexList.emplace_back(Vertex{i ,j, score});
                continue;
            }
            if (currentFormula.CheckIfOtherIsSubFormula(formula)) {
                const MolecularFormula loss = MolecularFormula(formula - currentFormula);
                const double score = (1 - loss.GetMass()/parentMass) + decompositionScores[i];
                graph.AddEdge(j, i);
                vertexList.emplace_back(Vertex{j ,i, score});
            }

        }
   }
    molecularNodeList = fragmentationNodes;
}

void FragmentationTree::SortVertexList() {
    std::sort(vertexList.begin(), vertexList.end(), CompareVertexes());
}

