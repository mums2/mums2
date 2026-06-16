//
// Created by Gregory Johnson on 11/19/25.
//
#include <Rcpp.h>
#include <string>
#include <vector>

#include "AnnotationStructure/AnnotationController.h"
#include "Spectra/ReadSpectra.h"


// [[Rcpp::export]]
Rcpp::List ReadMgf(const std::string& path, const bool interactive) {
    return ReadSpectra::ReadMGF(path, interactive);
}

// [[Rcpp::export]]
SEXP ReadMsp(const std::string& path, const bool interactive) {
    const std::vector<AnnotationNode> annotationData = ReadSpectra::ReadMSP(path, interactive);
    auto* controller = new AnnotationController(annotationData);
    return Rcpp::XPtr<AnnotationController>(controller);
}
