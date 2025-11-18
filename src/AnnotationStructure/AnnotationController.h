//
// Created by Gregory Johnson on 11/13/25.
//

#ifndef MUMS2_ANNOTATIONCONTROLLER_H
#define MUMS2_ANNOTATIONCONTROLLER_H
#include <queue>
#include "AnnotatedNode.h"
#include "AnnotationNode.h"
#include "../ScoringMethods/ScoringFactory.h"


class AnnotationController {
public:
    AnnotationController() = default;
    ~AnnotationController() = default;
    explicit AnnotationController(const std::vector<AnnotationNode>& annotations);
    bool AddNodes(const std::vector<AnnotationNode>& nodes);
    AnnotationNode GetNode(int index);
    size_t NodeCount() const {return annotations.size();}
    std::vector<AnnotationNode> GetNodes(const std::vector<int>& index) const;
    std::queue<AnnotatedNode> AnnotateFeature(const std::vector<Feature>& features, const ScoringFactory& factory,
        double minScoreThreshold, double chemicalMinScore, double precursorThreshold, size_t minPeaks,
        int threadCount) const;

private:
    std::vector<AnnotationNode> annotations;
};


#endif //MUMS2_ANNOTATIONCONTROLLER_H