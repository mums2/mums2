//
// Created by Gregory Johnson on 12/20/24.
//
#include <iostream>
#include <Rcpp.h>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include "DataStructures/CommunityMatrix.h"
#include "DiversityMetrics/Diversity.h"
#include "Rarefy/Rarefaction.h"
#include "DiversityMetrics/DiversityMetricFactory.h"

Rcpp::NumericMatrix CalculateDiversity(const Rcpp::NumericMatrix& abundances, const std::string& diversityIndex) {
    std::string index = diversityIndex;
    std::transform(index.begin(), index.end(), index.begin(), tolower);
    Diversity* diversity = DiversityMetricFactory::ChooseDiversityBasedOnIndex(index);
    if(diversity == nullptr) {
        Rcpp::stop("Diversity Metric not found");
    }
    Rcpp::NumericMatrix results = diversity->CalculateDiversity(abundances, index);
    delete diversity;
    return results;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix CalculateDiversityCommunityObject(const SEXP& communityMatrix, const std::string& diversityIndex) {
    const Rcpp::XPtr<CommunityMatrix> matrix(communityMatrix);
    const Rcpp::NumericMatrix abundances = matrix.get()->GetCommunityMatrix();
    std::string index = diversityIndex;
    std::transform(index.begin(), index.end(), index.begin(), tolower);
    Diversity* diversity = DiversityMetricFactory::ChooseDiversityBasedOnIndex(index);
    if(diversity == nullptr) {
        Rcpp::stop("Diversity Metric not found");
    }
    Rcpp::NumericMatrix results = diversity->CalculateDiversity(abundances, index);
    delete diversity;
    return results;
}

// [[Rcpp::export]]
SEXP CreateCommunityMatrix(Rcpp::NumericMatrix communityMatrix) {
    auto* matrix = new CommunityMatrix(communityMatrix);
    matrix->InitializeMatrix();
    return Rcpp::XPtr<CommunityMatrix>(matrix);
}

// [[Rcpp::export]]
SEXP GetCommunityMatrix(SEXP communityMatrix) {
    const Rcpp::XPtr<CommunityMatrix> matrix(communityMatrix);
    return matrix.get()->GetCommunityMatrix();
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RarefactionCalculation(const SEXP& communityMatrix, const uint32_t size,
    const uint32_t threshold) {

    const Rcpp::XPtr<CommunityMatrix> matrix(communityMatrix);
    const int row = matrix.get()->GetRow();
    const int col = matrix.get()->GetColumn();
    const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);
    Rarefaction rarefaction;
    std::vector<std::vector<uint32_t>>& allIndexes = matrix.get()->GetAllIndexes();
    const std::vector<std::vector<uint32_t>>& eligibleIndexes = matrix.get()->GetColumnEligibleIndexes();
    std::vector<std::vector<uint32_t>> eligibleAbundances = matrix.get()->GetRowAbundances();
    const std::vector<uint32_t>& sums = matrix.get()->GetSums();
    Rcpp::NumericMatrix resultantMatrix(row, col);
    for(int i = 0; i < row; i++) {
        const std::vector<uint32_t> communityVector = matrix.get()->GetCommunityMatrixByRow(i);
        const auto results = rarefaction.Rarefy(communityVector, eligibleIndexes[i], allIndexes[i],
                                                 size, sums[i], threshold);

        for(const auto& index : eligibleIndexes[i]) {
            resultantMatrix(i, index) = results.at(index);
        }

    }
    Rcpp::rownames(resultantMatrix) = samples;
    return resultantMatrix;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix FasterAvgDist(const SEXP& communityMatrix, const std::string& index,
    const uint32_t size, const uint32_t threshold, const int iterations = 1000) {
    const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);
    const Rcpp::XPtr<CommunityMatrix> communityObject(communityMatrix);

    Rcpp::NumericMatrix diversity = CalculateDiversity(RarefactionCalculation(communityObject,
        size, threshold), index);
    for(int i = 1; i < iterations; i++) {
        Rcpp::NumericMatrix rarefyMatrix = RarefactionCalculation(communityObject,
            size, threshold);
        diversity += CalculateDiversity(rarefyMatrix, index);
    }
    diversity = diversity/iterations;
    Rcpp::colnames(diversity) = samples;
    Rcpp::rownames(diversity) = samples;

    return diversity;
}

// [[Rcpp::export]]
void SizeOfCommunityObject(const SEXP& communityMatrix) {
    const Rcpp::XPtr<CommunityMatrix> communityObject(communityMatrix);
    communityObject->SeeSizes();
}