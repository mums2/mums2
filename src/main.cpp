//
// Created by Gregory Johnson on 12/20/24.
//
#include <iostream>
#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <numeric>
// #include "../../../../Downloads/gperftools-2.15/src/gperftools/profiler.h"
#include "Rarefy/Rarefaction.h"
#include "DiversityMetrics/DiversityCalculator.h"
#include "DiversityMetrics/DiversityMetricFactory.h"

// [[Rcpp::export]]
std::vector<double> CalculateDiversity(const Rcpp::NumericMatrix& abundances, std::string& diversityIndex) {
    // Beta diversity requires two vectors
    // Alpha diversity requires one
    // Rarefaction requires one but also other parameters
    // Should I make a class to deal with the parametrization?
    std::transform(diversityIndex.begin(), diversityIndex.end(), diversityIndex.begin(), tolower);
    const DiversityCalculator* diversity = DiversityMetricFactory::ChooseDiversityMetricBasedOnName(diversityIndex);
    const int size = abundances.nrow();
    const Rcpp::CharacterVector samples = Rcpp::rownames(abundances);
    std::vector<double> results(size);
    for(int i = 0; i < size; i++) {
        Rcpp::NumericVector abundance = abundances(i, Rcpp::_);
        std::vector<double> diversities = Rcpp::as<std::vector<double>>(abundance);
        results[i] = diversity->Calculate({diversities});
    }
    return results;
}
// [[Rcpp::export]]
double Test() {
    Rcpp::NumericVector vec(5);
    auto sum = std::accumulate(vec.begin(), vec.end(), 0.0);
    return sum;
}
void DiversityTest() {
    Diversity* diversity = DiversityMetricFactory::ChooseDiversityBasedOnIndex("shannon");
    //diversity->CalculateDiversity(matrix, "shannon")

}


// [[Rcpp::export]]
std::vector<std::vector<int64_t>> RarefactionCalculation(const Rcpp::NumericMatrix& communityMatrix, const int64_t size,
    const int64_t threshold) {
    const int row = communityMatrix.nrow();
    const int col = communityMatrix.ncol();

    std::vector<int> indexToName(col);
    std::iota(indexToName.begin(), indexToName.end(), 0);
    Rarefaction rarefaction;
    std::vector<std::vector<int64_t>> communityList(row, std::vector<int64_t>());
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
        communityList[i] = rarefaction.Rarefy(indexToName, communityVector, size, threshold);
    }
    return communityList;
}

// [[Rcpp::export]]
std::vector<std::vector<int>> TestMat() {
    return std::vector<std::vector<int>>(10, std::vector<int>(10, 0));
}
