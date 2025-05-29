//
// Created by Gregory Johnson on 5/29/25.
//

#ifndef ANNOTATIONNODE_H
#define ANNOTATIONNODE_H
#include <Rcpp.h>

#include <utility>
#include "../../Distance/DataStructures/Spectra.h"

class AnnotationNode {
public:
    AnnotationNode():spectra() {}
    AnnotationNode(const Rcpp::NumericVector peakMz,
        const Rcpp::NumericVector peakIntensities,
        const Rcpp::String& formula,
        const double precursorMz, const int index):
        spectra(Spectra("", Rcpp::as<std::vector<double>>(peakMz),
                            Rcpp::as<std::vector<double>>(peakIntensities), precursorMz)),
                                   molecularFormula(formula),
                                   index(index),
                                   precursorMz(precursorMz) {}

    const Spectra& GetSpectra() const { return spectra; }
    int GetIndex() const { return index; }
    double GetPrecursorMz() const { return precursorMz; }
    const Rcpp::String& GetFormula() const { return molecularFormula; }
private:
    Spectra spectra;
    Rcpp::String molecularFormula;
    int index;
    double precursorMz;
};
#endif //ANNOTATIONNODE_H
