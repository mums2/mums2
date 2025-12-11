//
// Created by Gregory Johnson on 11/13/25.
//

#include "AnnotationStructure/AnnotationController.h"
#include "Chemicals/MolecularFormula/MolecularFormulaSimilarity.h"
#include "CustomProgressBar/CliProgressBar.h"


AnnotationController::AnnotationController(const std::vector<AnnotationNode> &annotations):
annotations(annotations) {}

bool AnnotationController::AddNodes(const std::vector<AnnotationNode> &nodes) {
    annotations.insert(annotations.end(), nodes.begin(), nodes.end());
    return true;
}

bool AnnotationController::AddNodes(const AnnotationController &node) {
    annotations.insert(annotations.end(), node.annotations.begin(), node.annotations.end());
    return true;
}

AnnotationNode AnnotationController::GetNode(const int index) {
    return annotations[index];
}

std::vector<AnnotationNode> AnnotationController::GetNodes(const std::vector<int>& index) const {
    std::vector<AnnotationNode> nodes(index.size());
    for (size_t i = 0; i < index.size(); ++i) {
        nodes[i] = annotations[index[i]];
    }
    return nodes;
}
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
std::queue<AnnotatedNode> AnnotationController::AnnotateFeature(const std::vector<Feature> &features,
const ScoringFactory& factory, const double minScoreThreshold, const double chemicalMinScore,
const double precursorThreshold, const size_t minPeaks, const int threadCount) const {

    std::queue<AnnotatedNode> result;
    CliProgressBar progressBar;
    int currentProgress = 0;
    const int size = static_cast<int>(annotations.size());
    const size_t featureSize = features.size();
    const auto maxProgress = static_cast<float>(featureSize);
    std::mutex mutex;
    for (size_t i = 0; i < featureSize; ++i) {
        const Feature& feature = features[i];
        RcppThread::parallelFor(0, size, [this, &feature, &factory, &result,
            &minScoreThreshold, &chemicalMinScore, &precursorThreshold, &minPeaks, &mutex, &i](int j) {
            const AnnotationNode& node = annotations[j];
            if ((std::abs(feature.mz - node.precursorMz)) * 1e6 / feature.mz <= precursorThreshold) {
                const double chemicalSimilarity = MolecularFormulaSimilarity::ComputeSimilarity(feature.formula,
               node.chemicalFormula);
                if (chemicalSimilarity >= chemicalMinScore) {
                    const double score = factory.CalculateScore(feature.spectra, node.spectra, minPeaks);
                    if (score >= minScoreThreshold) {
                        AnnotatedNode annotatedNode;
                        annotatedNode.featureIndex = i;
                        annotatedNode.annotationNodeIndex = j;
                        annotatedNode.formulaSimilarity = chemicalSimilarity;
                        annotatedNode.score = score;
                        mutex.lock();
                        result.push(annotatedNode);
                        mutex.unlock();
                    }
                }
            }
        }, threadCount);
        progressBar.update(static_cast<float>(currentProgress++)/maxProgress);
    }
    progressBar.end_display();
    return result;
}

