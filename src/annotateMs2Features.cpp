#include "Utils/Utils.h"
#include <Rcpp.h>
#include <algorithm>
#include "AnnotationStructure/Annotation.h"
#include "AnnotationStructure/AnnotationController.h"
#include "AnnotationStructure/Feature.h"
#include "Chemicals/MolecularFormula/MolecularFormulaSimilarity.h"
#include "CustomProgressBar/CliProgressBar.h"
#include "ScoringMethods/ScoringFactory.h"

// [[Rcpp::export]]
Rcpp::DataFrame AnnotateMs2Features(const Rcpp::DataFrame& queryList, const Rcpp::List querySpectra,
    const SEXP annotationController, const Rcpp::List& scoringParameters, const Rcpp::StringVector& formulas,
    const double precursorThreshold,const double minScoreThreshold, const double chemicalMinScore,
    const size_t minPeaks, const int threadCount) {
    const ScoringFactory factory(scoringParameters);
    const auto querySpectraSize = static_cast<size_t>(querySpectra.size());
    std::vector<Feature> queryFeatures(querySpectraSize);
    const Rcpp::NumericVector& queryMz = queryList["mz"];
    const Rcpp::NumericVector& queryRt = queryList["rt"];
    const Rcpp::StringVector& queryMs1Id = queryList["ms1_compound_id"];
    const Rcpp::StringVector& queryMs2Id = queryList["ms2_spectrum_id"];

    for (size_t i = 0; i < querySpectraSize; ++i) {
        Feature feature;
        const Rcpp::List& spectra = querySpectra[i];
        feature.formula = formulas[i];
        feature.ms1_id = queryMs1Id[i];
        feature.ms2_id = queryMs2Id[i];
        feature.mz = queryMz[i];
        feature.rt = queryRt[i];
        feature.spectra = Spectra("", spectra["mz"], spectra["intensity"], feature.mz);
        queryFeatures[i] = feature;
    }

    const Rcpp::XPtr<AnnotationController> ptr(annotationController);
    const std::queue<AnnotatedNode> results = ptr.get()->AnnotateFeature(queryFeatures, factory, minScoreThreshold,
        chemicalMinScore, precursorThreshold, minPeaks, threadCount);
    Annotation annotation(results);
    return annotation.CreateAnnotationDataFrame();
}

// [[Rcpp::export]]
size_t GetNodeCount(const SEXP& annotationController) {
    const Rcpp::XPtr<AnnotationController> ptr(annotationController);
    return ptr.get()->NodeCount();
}

// [[Rcpp::export]]
SEXP GetNode(const SEXP& annotationController, const int index) {
    const Rcpp::XPtr<AnnotationController> ptr(annotationController);
    const AnnotationNode data = ptr.get()->GetNode(index);
    return Rcpp::List::create(Rcpp::Named("info") =
                Rcpp::List::create(Rcpp::Named("keys") = data.keys,
                    Rcpp::Named("values") = data.values),
            Rcpp::Named("spec") = Rcpp::List::create(Rcpp::Named("mz") = data.spectra.mz,
                Rcpp::Named("intensity") = data.spectra.intensity));


}

// [[Rcpp::export]]
SEXP CombineReferenceDatabases(const SEXP& annotationController, const SEXP& otherAnnotationController) {
    const Rcpp::XPtr<AnnotationController> ptr(annotationController);
    const Rcpp::XPtr<AnnotationController> other(otherAnnotationController);
    const AnnotationController& controller = *ptr.get();
    const AnnotationController& otherController = *other.get();
    auto* newController = new AnnotationController();
    newController->AddNodes(controller);
    newController->AddNodes(otherController);
    return Rcpp::XPtr<AnnotationController>(newController);
}