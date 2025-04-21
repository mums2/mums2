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
#include "Spectra/ReadSpectra.h"

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
    const Rcpp::CharacterVector& rowNames = matrix.get()->GetRowNames();
    const Rcpp::CharacterVector& columnNames = matrix.get()->GetColumnNames();

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
    Rcpp::rownames(resultantMatrix) = rowNames;
    Rcpp::colnames(resultantMatrix) = columnNames;
    return resultantMatrix;
}

// [[Rcpp::export]]
Rcpp::DataFrame FasterAvgDist(const SEXP& communityMatrix, const std::string& index,
    const uint32_t size, const uint32_t threshold, const int iterations = 1000) {
    const Rcpp::XPtr<CommunityMatrix> communityObject(communityMatrix);
    const Rcpp::CharacterVector samples = communityObject.get()->GetSampleNames();

    Rcpp::NumericMatrix diversityMatrix = CalculateDiversity(RarefactionCalculation(communityObject,
        size, threshold), index);
    for(int i = 1; i < iterations; i++) {
        Rcpp::NumericMatrix rarefyMatrix = RarefactionCalculation(communityObject,
            size, threshold);
        diversityMatrix += CalculateDiversity(rarefyMatrix, index);
    }
    diversityMatrix = diversityMatrix/iterations;
    Rcpp::colnames(diversityMatrix) = samples;
    Rcpp::rownames(diversityMatrix) = samples;
    const int sampleSize = std::pow(samples.size(), 2);
    Rcpp::CharacterVector firstSample(sampleSize);
    Rcpp::CharacterVector otherSample(sampleSize);
    Rcpp::NumericVector diversityResults(sampleSize);
    int currentIndex = 0;
    for(int i = 0; i < samples.size(); i++) {
        for(int j = 0; j < samples.size(); j++) {
            firstSample[currentIndex] = samples[i];
            otherSample[currentIndex] = samples[j];
            diversityResults[currentIndex++] = diversityMatrix(i, j);
        }
    }

    return Rcpp::DataFrame::create(Rcpp::Named("firstSample") = firstSample,
        Rcpp::Named("otherSample") = otherSample,
        Rcpp::Named("diversity") = diversityResults);
}

// [[Rcpp::export]]
Rcpp::List ReadMgf(const std::string& path) {
    ReadSpectra spectra;
    return(spectra.ReadMGF(path));
}

// [[Rcpp::export]]
Rcpp::List ReadMsp(const std::string& path) {
    ReadSpectra spectra;
    return(spectra.ReadMSP(path));
}

// [[Rcpp::export]]
double DotProduct(Rcpp::NumericVector x, Rcpp::NumericVector y) {
    double result = 0;
    for (int i = 0; i < x.length(); i++) {
        result += x[i] * y[i];
    }
    return result;
}