//
// Created by gregj on 1/9/2025.
//

#include "MothurShannonCalculator.h"

double MothurShannonCalculator::Calculate(const Rcpp::List &vectors) const {
    // Maybe I send over another vector
    if (vectors.size() < 2)
        return 0;
    const Rcpp::NumericVector& abundance = vectors[0];
    //each position in the sequence abundance is correlated to the abundance
    // at each i'th position, it represents the amount of sequences
    const Rcpp::NumericVector& sequenceAbundance = vectors[1];
    for (size_t i = 0; i < abundance.size(); i++) {

    }
    return 0;
}
