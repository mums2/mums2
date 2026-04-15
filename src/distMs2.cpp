#include "Distance/Distance.h"
#include <Rcpp.h>

//[[Rcpp::export]]
Rcpp::DataFrame  distMS2(const Rcpp::List spectraDataList, const Rcpp::List parameters,
    const double precursor_thresh, const double cutoff, const int minPeaks, const int numberOfThreads){

    Distance distance;
    distance.CreateSpectraList(spectraDataList);
    const ScoringFactory factory(parameters);
    distance.CalculateDistances(precursor_thresh, cutoff, factory, minPeaks, numberOfThreads);
    Rcpp::DataFrame dist = distance.ExtractMatrix();

    return dist;
}

//[[Rcpp::export]]
double CompareMs2(const std::vector<double>& mz1, const std::vector<double>& mz2,
    const std::vector<double>& int1, const std::vector<double>& int2, const double precMz1,
    const double precMz2, const Rcpp::List parameters) {
    Spectra spectra("", mz1, int1, precMz1);
    Spectra spectra2("", mz2, int2, precMz2);
    const ScoringFactory factory(parameters);
    return factory.CalculateScore(spectra, spectra2, 0);
}
