#ifndef DISTANCE
#define DISTANCE

#include "DataStructures/Spectra.h"
#include "DataStructures/SparseValue.h"
#include <vector>
#include <iostream>
#include <Rcpp.h>
#include <cmath>
#include "../ScoringMethods/ScoringFactory.h"


class Distance final {
public:

    Distance() = default;
    ~Distance() = default;
    void CreateSpectraList(Rcpp::List data);
    void CalculateDistances(double prec_threshold, double cutoff, const ScoringFactory &scoreMethod, int minPeaks);
    Rcpp::DataFrame ExtractMatrix();

private:
    std::vector<Spectra> spectraList;
    std::queue<SparseValue> sparseMatrix;
};

#endif //DISTANCE