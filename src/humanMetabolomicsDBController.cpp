//
// Created by Gregory Johnson on 11/19/25.
//
#include "AnnotationStructure/AnnotationController.h"
#include "HMDB/HumanMetabolomicsDB.h"
#include "CustomProgressBar/CliProgressBar.h"
#include <Rcpp.h>

// [[Rcpp::export]]
SEXP CreateHumanMetabolomicsDB(const size_t nodeSize) {
    auto* hmdb = new HumanMetabolomicsDB(nodeSize);
    return Rcpp::XPtr<HumanMetabolomicsDB>(hmdb);
}

// [[Rcpp::export]]
void AddHumanMetabolomicNode(SEXP& hmdbPtr, const std::vector<std::string>& names,
    const std::vector<std::string>& values, const size_t index) {
    Rcpp::XPtr<HumanMetabolomicsDB> hmdbPointer(hmdbPtr);
    hmdbPointer.get()->AddHumanMetabolomicNode(HumanMetabolomicsDBNode(names,values), index);
}

// [[Rcpp::export]]
void AddSpectra(SEXP& hmdbPtr, const std::vector<std::string>& fileNames,
    const std::vector<std::string>& databaseNames) {
    Rcpp::XPtr<HumanMetabolomicsDB> hmdbPointer(hmdbPtr);
    CliProgressBar progressBar;
    float counter = 0;
    const float maxSize = static_cast<float>(fileNames.size());
    for (size_t i = 0; i < fileNames.size(); i++) {
        hmdbPointer.get()->AddSpectraFiles(fileNames[i], databaseNames[i]);
        counter++;
        progressBar.update(counter/maxSize);
    }
}

// [[Rcpp::export]]
void ProcessMs2Files(SEXP& hmdbPtr) {
    Rcpp::XPtr<HumanMetabolomicsDB> hmdbPointer(hmdbPtr);
    hmdbPointer.get()->ProcessSpectraFiles();
}

// [[Rcpp::export]]
SEXP CreateAnnotationController(SEXP& hmdbPtr) {
    Rcpp::XPtr<HumanMetabolomicsDB> hmdbPointer(hmdbPtr);
    AnnotationController* node = hmdbPointer.get()->ConstructDataBase();
    Rcpp::XPtr<AnnotationController> annotationPtr(node);
    return annotationPtr;
}