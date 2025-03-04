//
// Created by gregj on 2/18/2025.
//

#include "DataStructures/CommunityMatrix.h"
#include <numeric>
CommunityMatrix::CommunityMatrix(const Rcpp::NumericMatrix &matrix) :communityMatrix(matrix) {}



void CommunityMatrix::InitializeMatrix() {

    row = communityMatrix.nrow();
    col = communityMatrix.ncol();
    rowAbundance = std::vector<std::vector<int64_t>>(row);
    eligibleRowIndexes = std::vector<std::vector<int64_t>>(row);
    allIndexes = std::vector<std::vector<int64_t>>(row);
    sums = std::vector<int64_t>(row);
    Rcpp::NumericMatrix resultantMatrix(row, col);
    for(int i = 0; i < row; i++) {
        Rcpp::NumericVector community = communityMatrix(i, Rcpp::_);
        std::vector<int64_t> communityVector = Rcpp::as<std::vector<int64_t>>(community);
        const size_t size = communityVector.size();
        eligibleRowIndexes[i].reserve(size);
        rowAbundance[i].reserve(size);
        sums[i] = std::accumulate(communityVector.begin(), communityVector.end(), 0LL);
        allIndexes[i] = std::vector<int64_t>(sums[i]);
        size_t index = 0;
        int64_t currentPosition = 0;
        for(size_t j = 0; j < size; j++) {
            int64_t val = communityVector[j];
            if(val > 0) {
                eligibleRowIndexes[i].emplace_back(j);
                rowAbundance[i].emplace_back(val);
                for (int64_t k = 0; k < val; k++) {
                    allIndexes[i][index++] = currentPosition;
                }
                currentPosition++;
            }
        }
    }
}

std::vector<int64_t> CommunityMatrix::GetCommunityMatrixByRow(const int row) const {
    Rcpp::NumericVector community = communityMatrix(row, Rcpp::_);
    const std::vector<int64_t> communityVector = Rcpp::as<std::vector<int64_t>>(community);
    return communityVector;
}
