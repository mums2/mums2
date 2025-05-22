#include "AnnotateMs2/AnnotateMs2.h"
#include "Utils/Utils.h"
#include <Rcpp.h>
#include "ScoringMethods/ScoringFactory.h"

// [[Rcpp::export]]
Rcpp::DataFrame AnnotateMs2Features(const std::vector<std::string>& variableId, const std::vector<std::string>& ms2Id, 
          const std::vector<float>& ms2Mz, const std::vector<float>& ms2Rt, const Rcpp::List& ms2Spectra,
          const Rcpp::List& reference, const Rcpp::List& parameters, const double precursorThreshold, double minScore,
          const size_t minPeaks) {

    const ScoringFactory factory(parameters);
    AnnotateMs2 annotateMs2(minPeaks);
    annotateMs2.createQueryList(variableId, ms2Id, ms2Mz, ms2Rt, ms2Spectra);
    annotateMs2.createRefList(reference);
    Rcpp::DataFrame matches = annotateMs2.getMatches(precursorThreshold, factory, minScore);
    
    return matches;

    
}