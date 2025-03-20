#include "Distance/Distance.h"
#include <Rcpp.h>

//[[Rcpp::export]]
Rcpp::DataFrame  distMS2(const Rcpp::List spectraDataList, const Rcpp::List parameters, double precursor_thresh, double cutoff){

    Distance distance;
    distance.CreateSpectraList(spectraDataList);
    const ScoringFactory factory(parameters);
    distance.CalculateDistances(precursor_thresh, cutoff, factory);
    Rcpp::DataFrame dist = distance.ExtractMatrix();

    return dist;
}