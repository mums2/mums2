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
    communityAbundances = std::vector<std::vector<uint32_t>>(row);
    eligibleRowIndexes = std::vector<std::vector<uint32_t>>(row);
    allIndexes = std::vector<std::vector<uint32_t>>(row);
    sums = std::vector<uint32_t>(row);
    Rcpp::NumericMatrix resultantMatrix(row, col);
    for(int i = 0; i < row; i++) {
        Rcpp::NumericVector community = communityMatrix(i, Rcpp::_);
        std::vector<uint32_t> communityVector = Rcpp::as<std::vector<uint32_t>>(community);
        communityAbundances[i] = communityVector;
        const size_t size = communityVector.size();
        eligibleRowIndexes[i].reserve(size);
        sums[i] = std::accumulate(communityVector.begin(), communityVector.end(), 0LL);
        allIndexes[i] = std::vector<uint32_t>(sums[i]);
        size_t index = 0;
        uint32_t currentPosition = 0;
        for(size_t j = 0; j < size; j++) {
            uint32_t val = communityVector[j];
            if(val > 0) {
                eligibleRowIndexes[i].emplace_back(j);
                for (uint32_t k = 0; k < val; k++) {
                    allIndexes[i][index++] = currentPosition;
                }
                currentPosition++;
            }
        }
    }
}

std::vector<uint32_t> CommunityMatrix::GetCommunityMatrixByRow(const int row) const {
    Rcpp::NumericVector community = communityMatrix(row, Rcpp::_);
    const std::vector<uint32_t> communityVector = Rcpp::as<std::vector<uint32_t>>(community);
    return communityVector;
}

Rcpp::CharacterVector CommunityMatrix::GetSampleNames() {
    return Rcpp::rownames(communityMatrix);
}