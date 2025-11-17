//
// Created by Gregory Johnson on 11/13/25.
//

#ifndef MUMS2_ANNOTATIONCONTROLLER_H
#define MUMS2_ANNOTATIONCONTROLLER_H
#include "AnnotatedNode.h"
#include "AnnotationNode.h"
#include "../ScoringMethods/ScoringFactory.h"


class AnnotationController {
public:
    AnnotationController() = default;
    explicit AnnotationController(const std::vector<AnnotationNodeData>& annotations,
        double minScoreThreshold, double chemicalMinScore, double precursorThreshold, size_t minPeaks);
    bool AddNodes(const std::vector<AnnotationNodeData>& nodes);
    AnnotationNodeData GetNode(int index);
    std::vector<AnnotationNodeData> GetNodes(const std::vector<int>& index) const;
    std::queue<AnnotatedNode> AnnotateFeature(const Feature& feature, const ScoringFactory& factory);

private:
    std::vector<AnnotationNodeData> annotations;
    double minScoreThreshold, chemicalMinScore, precursorThreshold;
    size_t minPeaks;
};


#endif //MUMS2_ANNOTATIONCONTROLLER_H