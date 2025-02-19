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
    sums = std::vector<int64_t>(row);
    const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);
    std::vector<int> indexToName(col);
    Rcpp::NumericMatrix resultantMatrix(row, col);
    for(int i = 0; i < row; i++) {
        Rcpp::NumericVector community = communityMatrix(i, Rcpp::_);
        std::vector<int64_t> communityVector = Rcpp::as<std::vector<int64_t>>(community);
        size_t communityVectorSize = communityVector.size();
        eligibleRowIndexes[i].reserve(communityVectorSize);
        rowAbundance[i].reserve(communityVectorSize);
        sums[i] = std::accumulate(communityVector.begin(), communityVector.end(), 0LL);
        for(size_t j = 0; j < communityVectorSize; j++) {
            int64_t val = communityVector[j];
            if(val > 0) {
                eligibleRowIndexes[i].emplace_back(j);
                rowAbundance[i].emplace_back(val);
            }
        }
    }

}

const std::vector<int64_t>& CommunityMatrix::GetCommunityMatrixByRow(int row) const {
    Rcpp::NumericVector community = communityMatrix(row, Rcpp::_);
    const std::vector<int64_t>& communityVector = Rcpp::as<std::vector<int64_t>>(community);
    return communityVector;
}
