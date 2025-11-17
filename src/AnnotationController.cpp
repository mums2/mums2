//
// Created by Gregory Johnson on 11/13/25.
//

#include "AnnotationStructure/AnnotationController.h"
#include "Chemicals/MolecularFormula/MolecularFormulaSimilarity.h"
#include "CustomProgressBar/CliProgressBar.h"


AnnotationController::AnnotationController(const std::vector<AnnotationNodeData> &annotations):
annotations(annotations) {}

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

std::queue<AnnotatedNode> AnnotationController::AnnotateFeature(const std::vector<Feature> &features,
const ScoringFactory& factory, const double minScoreThreshold, const double chemicalMinScore,
const double precursorThreshold, const size_t minPeaks) const {

    std::queue<AnnotatedNode> result;
    CliProgressBar progressBar;
    int currentProgress = 0;
    const auto maxProgress = static_cast<float>(features.size());
    for (const auto& feature: features) {
        for (const auto& node : annotations) {
            double val = (std::abs(feature.mz - node.precursorMz)) * 1e6 / feature.mz ;
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
        progressBar.update(static_cast<float>(currentProgress++)/maxProgress);
    }
    progressBar.end_display();
    return result;
}

