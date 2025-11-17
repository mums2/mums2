//
// Created by Gregory Johnson on 11/13/25.
//

#include "AnnotationStructure/AnnotationController.h"

#include "Chemicals/MolecularFormula/MolecularFormulaSimilarity.h"




AnnotationController::AnnotationController(const std::vector<AnnotationNodeData> &annotations, const double minScoreThreshold,
                                           const double chemicalMinScore, const double precursorThreshold, const size_t minPeaks): annotations(annotations), minScoreThreshold(minScoreThreshold),
chemicalMinScore(chemicalMinScore), precursorThreshold(precursorThreshold), minPeaks(minPeaks) {}

bool AnnotationController::AddNodes(const std::vector<AnnotationNodeData> &nodes) {
    annotations.insert(annotations.end(), nodes.begin(), nodes.end());
    return true;
}

AnnotationNodeData AnnotationController::GetNode(const int index) {
    return annotations[index];
}

std::vector<AnnotationNodeData> AnnotationController::GetNodes(const std::vector<int>& index) const {
    std::vector<AnnotationNodeData> nodes(index.size());
    for (size_t i = 0; i < index.size(); ++i) {
        nodes[i] = annotations[index[i]];
    }
    return nodes;
}

std::queue<AnnotatedNode> AnnotationController::AnnotateFeature(const Feature &feature,
    const ScoringFactory& factory) {
    std::queue<AnnotatedNode> result;
    for (const auto& node : annotations) {
        if ((std::abs(feature.mz - node.precursorMz)) * 1e6 / feature.mz > precursorThreshold) continue;
        const double chemicalSimilarity = MolecularFormulaSimilarity::ComputeSimilarity(Rcpp::wrap(feature.formula),
            Rcpp::wrap(node.chemicalFormula));
        if (chemicalSimilarity < chemicalMinScore) continue;
        const double score = factory.CalculateScore(feature.spectra, node.spectra, minPeaks);
        if (score < minScoreThreshold) continue;
        AnnotatedNode annotatedNode;
        annotatedNode.feature = feature;
        annotatedNode.node = node;
        annotatedNode.formulaSimilarity = chemicalSimilarity;
        annotatedNode.score = score;
        result.push(annotatedNode);
    }
    return result;
}

