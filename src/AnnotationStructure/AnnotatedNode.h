//
// Created by Gregory Johnson on 11/13/25.
//

#ifndef MUMS2_ANNOTATEDNODE_H
#define MUMS2_ANNOTATEDNODE_H
#include <utility>

#include "Feature.h"
#include "AnnotationNode.h"

struct AnnotatedNode {
    AnnotatedNode() = default;
    size_t featureIndex = 0;
    size_t annotationNodeIndex = 0;
    double score = 0;
    double formulaSimilarity = 0;

};
#endif //MUMS2_ANNOTATEDNODE_H