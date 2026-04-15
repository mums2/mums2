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
