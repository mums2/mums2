//
// Created by Gregory Johnson on 11/13/25.
//

#ifndef MUMS2_ANNOTATION_H
#define MUMS2_ANNOTATION_H
#include "AnnotatedNode.h"
#include <Rcpp.h>
#include <queue>


class Annotation {
public:
    explicit Annotation(const std::queue<AnnotatedNode>&);
    Rcpp::DataFrame CreateAnnotationDataFrame();
private:
    std::queue<AnnotatedNode> annotatedNodes;
};


#endif //MUMS2_ANNOTATION_H