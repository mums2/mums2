//
// Created by Gregory Johnson on 11/19/25.
//
#include "AnnotationStructure/AnnotationController.h"
#include "HMDB/HumanMetabolomicsDB.h"
#include <Rcpp.h>

// [[Rcpp::export]]
SEXP CreateHumanMetabolomicsDB() {
    auto* hmdb = new HumanMetabolomicsDB();
    return Rcpp::XPtr<HumanMetabolomicsDB>(hmdb);
}

// [[Rcpp::export]]
void AddHumanMetabolomicNode(SEXP& hmdbPtr, const std::vector<std::string>& names,
    const std::vector<std::string>& values) {
    Rcpp::XPtr<HumanMetabolomicsDB> hmdbPointer(hmdbPtr);
    hmdbPointer.get()->AddHumanMetabolomicNode(HumanMetabolomicsDBNode(names,values));
}

// [[Rcpp::export]]
void PrintHMDBNames(const SEXP& hmdbPtr) {
    Rcpp::XPtr<HumanMetabolomicsDB> hmdbPointer(hmdbPtr);
    hmdbPointer.get()->PrintHumanMetabolomicsDB();
}
// [[Rcpp::export]]
void AddSpectra(SEXP& hmdbPtr, const std::vector<std::string>& fileNames,
    const std::vector<std::string>& databaseNames) {
    Rcpp::XPtr<HumanMetabolomicsDB> hmdbPointer(hmdbPtr);
    for (size_t i = 0; i < fileNames.size(); i++) {
        hmdbPointer.get()->AddSpectraFiles(fileNames[i], databaseNames[i]);
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