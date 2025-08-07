//
// Created by gregj on 2/18/2025.
//

#include "DataStructures/CommunityMatrix.h"
#include <numeric>
CommunityMatrix::CommunityMatrix(const Rcpp::NumericMatrix &matrix) :communityMatrix(matrix) {}



void CommunityMatrix::InitializeMatrix() {

    rowNames = Rcpp::rownames(communityMatrix);
    colNames = Rcpp::colnames(communityMatrix);
    row = communityMatrix.nrow();
    col = communityMatrix.ncol();
    communityAbundances = std::vector<std::vector<uint64_t>>(row);
    eligibleRowIndexes = std::vector<std::vector<uint64_t>>(row);
    abundancesRanges = std::vector<std::vector<uint64_t>>(row);
    sums = std::vector<uint64_t>(row);
    auto cppMatrix = std::vector<double>(row * col);
    Rcpp::NumericMatrix resultantMatrix(row, col);
    for(int i = 0; i < row; i++) {
        Rcpp::NumericVector community = communityMatrix(i, Rcpp::_);
        std::vector<uint64_t> communityVector = Rcpp::as<std::vector<uint64_t>>(community);
        communityAbundances[i] = communityVector;
        const size_t size = communityVector.size();
        eligibleRowIndexes[i].reserve(size);
        abundancesRanges[i].reserve(size);
        sums[i] = std::accumulate(communityVector.begin(), communityVector.end(), 0LL);
        uint64_t previousAbundance = 0;
        for(size_t j = 0; j < size; j++) {
            uint64_t abundance = communityVector[j];
            if (abundance <= 0) continue;
            eligibleRowIndexes[i].emplace_back(j);
            previousAbundance = abundance + previousAbundance;
            abundancesRanges[i].emplace_back(previousAbundance);
        }
    }
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            cppMatrix[i * col + j] = static_cast<double>(communityAbundances[i][j]);
        }
    }
    cppMatrixOfAbundances = CppMatrix(cppMatrix, row, col);
}

std::vector<uint64_t> CommunityMatrix::GetCommunityMatrixByRow(const int row) const {
    Rcpp::NumericVector community = communityMatrix(row, Rcpp::_);
    const std::vector<uint64_t> communityVector = Rcpp::as<std::vector<uint64_t>>(community);
    return communityVector;
}

Rcpp::CharacterVector CommunityMatrix::GetSampleNames() {
    return Rcpp::rownames(communityMatrix);
}