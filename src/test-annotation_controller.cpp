//
// Created by Gregory Johnson on 11/13/25.
//
#include <testthat.h>
#include "AnnotationStructure/AnnotationController.h"
#include "AnnotationStructure/AnnotatedNode.h"
#include "AnnotationStructure/AnnotationNode.h"
#include "AnnotationStructure/Feature.h"
#include "ScoringMethods/ScoringFactory.h"
#include <queue>



// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("Test Annotation Controller") {

    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("The AddNode and AddNodes Function returns true") {
        AnnotationController anno;
        AnnotationNodeData node;
        std::vector<AnnotationNodeData> nodes;
        nodes.emplace_back(node);
        expect_true(anno.AddNodes(nodes));
    }

    test_that("GetNode Function returns node at the proper index") {
        AnnotationController anno;
        AnnotationNodeData node;
        std::vector<AnnotationNodeData> nodes;
        nodes.emplace_back(node);
        anno.AddNodes(nodes);
        auto node2 = anno.GetNode(0);
        expect_true(node2.name == node.name);
    }

    test_that("AnnotateFeature returns a list annotated feature") {
        AnnotationController anno;
        AnnotationNodeData node;
        std::vector<AnnotationNodeData> nodes;
        nodes.emplace_back(node);
        Feature feature;
        anno.AddNodes(nodes);
        Rcpp::List a = Rcpp::List::create();
        ScoringFactory score(a);
        std::queue<AnnotatedNode> annotations = anno.AnnotateFeature(feature, score);
        expect_true(annotations.empty());
    }

}
