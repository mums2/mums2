//
// Created by Gregory Johnson on 12/20/24.
//
#include <iostream>
#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include "DiversityMetrics/Diversity.h"
// #include "../../../../Downloads/gperftools-2.15/src/gperftools/profiler.h"
#include "Rarefy/Rarefaction.h"
#include "DiversityMetrics/DiversityCalculator.h"
#include "DiversityMetrics/DiversityMetricFactory.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix CalculateDiversity(const Rcpp::NumericMatrix& abundances, std::string& diversityIndex) {
    std::transform(diversityIndex.begin(), diversityIndex.end(), diversityIndex.begin(), tolower);
    Diversity* diversity = DiversityMetricFactory::ChooseDiversityBasedOnIndex(diversityIndex);
    if(diversity == nullptr) {
        Rcpp::stop("Diversity Metric not found");
    }
    Rcpp::NumericMatrix results = diversity->CalculateDiversity(abundances, diversityIndex);
    return results;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RarefactionCalculation(const Rcpp::NumericMatrix& communityMatrix, const int64_t size,
    const int64_t threshold) {
    const int row = communityMatrix.nrow();
    const int col = communityMatrix.ncol();

    std::vector<int> indexToName(col);
    std::iota(indexToName.begin(), indexToName.end(), 0);
    Rarefaction rarefaction;
    Rcpp::NumericMatrix resultantMatrix(row, col);
    for(int i = 0; i < row; i++) {
        Rcpp::NumericVector community = communityMatrix(i, Rcpp::_);
        std::vector<int64_t> communityVector = Rcpp::as<std::vector<int64_t>>(community);
        // Rcpp::Rcout << communityVector << std::endl;
        // Rarefy
        // We are going to have to switch between the tranpose of the community matrix and
        // the orignal matrix.
        // Samples are represented by rows and columns represent species.
        // So ex..
        //          FD39  FD09 <- is backwards            species1     <- Correct
        // Species1  0      1                       FD39    1
        //                                          DF09    0
        const auto results = rarefaction.Rarefy(indexToName, communityVector, size, threshold);
        for(int j = 0; j < col; j++) {
            resultantMatrix(i, j) = results[j];
        }
    }
    return resultantMatrix;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix FasterAvgDist(const Rcpp::NumericMatrix& communityMatrix, const std::string& index,
    const int64_t size, const int64_t threshold, const int iterations = 1000) {

    for(int i = 0; i < iterations; i++) {

    }
    return {};
}


