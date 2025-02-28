//
// Created by gregj on 2/18/2025.
//

#include "DataStructures/CommunityMatrix.h"
#include <numeric>
CommunityMatrix::CommunityMatrix(const Rcpp::NumericMatrix &matrix) :communityMatrix(matrix) {}



void CommunityMatrix::InitializeMatrix() {
    int64_t max = 0;
    size_t index = 0;
    row = communityMatrix.nrow();
    col = communityMatrix.ncol();
    rowAbundance = std::vector<std::vector<int64_t>>(row);
    eligibleRowIndexes = std::vector<std::vector<int64_t>>(row);
    abundanceRanges = std::vector<std::vector<int64_t>>(row);
    allIndexes = std::vector<std::vector<int64_t>>(row);
    sums = std::vector<int64_t>(row);
    const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);
    std::vector<int> indexToName(col);
    Rcpp::NumericMatrix resultantMatrix(row, col);
    // std::vector<int64_t> abundanceRanges(1,0);
    for(int i = 0; i < row; i++) {
        Rcpp::NumericVector community = communityMatrix(i, Rcpp::_);
        std::vector<int64_t> communityVector = Rcpp::as<std::vector<int64_t>>(community);
        const size_t size = communityVector.size();
        abundanceRanges[i].reserve(size + 1);
        abundanceRanges[i].emplace_back(0);
        eligibleRowIndexes[i].reserve(size);
        rowAbundance[i].reserve(size);
        sums[i] = std::accumulate(communityVector.begin(), communityVector.end(), 0LL);
        allIndexes[i] = std::vector<int64_t>(sums[i]);
        std::iota(allIndexes[i].begin(), allIndexes[i].end(), 0);
        int count = 1;
        for(size_t j = 0; j < size; j++) {
            int64_t val = communityVector[j];
            if(val > 0) {
                eligibleRowIndexes[i].emplace_back(j);
                rowAbundance[i].emplace_back(val);
                abundanceRanges[i].emplace_back(abundanceRanges[i][count++ - 1] + val);
            }
        }
        if (sums[i] > max) {
            max = sums[i];
            index = i;
        }

    }
}

std::vector<int64_t> CommunityMatrix::GetCommunityMatrixByRow(const int row) const {
    Rcpp::NumericVector community = communityMatrix(row, Rcpp::_);
    const std::vector<int64_t> communityVector = Rcpp::as<std::vector<int64_t>>(community);
    return communityVector;
}
