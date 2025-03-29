//
// Created by gregj on 2/18/2025.
//

#include "DataStructures/CommunityMatrix.h"
#include <numeric>
CommunityMatrix::CommunityMatrix(const Rcpp::NumericMatrix &matrix) :communityMatrix(matrix) {}



void CommunityMatrix::InitializeMatrix() {

    row = communityMatrix.nrow();
    col = communityMatrix.ncol();
    rowAbundance = std::vector<std::vector<uint>>(row);
    eligibleRowIndexes = std::vector<std::vector<uint>>(row);
    allIndexes = std::vector<std::vector<uint>>(row);
    sums = std::vector<uint>(row);
    Rcpp::NumericMatrix resultantMatrix(row, col);
    for(int i = 0; i < row; i++) {
        Rcpp::NumericVector community = communityMatrix(i, Rcpp::_);
        std::vector<uint> communityVector = Rcpp::as<std::vector<uint>>(community);
        const size_t size = communityVector.size();
        eligibleRowIndexes[i].reserve(size);
        rowAbundance[i].reserve(size);
        sums[i] = std::accumulate(communityVector.begin(), communityVector.end(), 0LL);
        allIndexes[i] = std::vector<uint>(sums[i]);
        size_t index = 0;
        uint currentPosition = 0;
        for(size_t j = 0; j < size; j++) {
            uint val = communityVector[j];
            if(val > 0) {
                eligibleRowIndexes[i].emplace_back(j);
                rowAbundance[i].emplace_back(val);
                for (uint k = 0; k < val; k++) {
                    allIndexes[i][index++] = currentPosition;
                }
                currentPosition++;
            }
        }
    }
}

std::vector<uint> CommunityMatrix::GetCommunityMatrixByRow(const int row) const {
    Rcpp::NumericVector community = communityMatrix(row, Rcpp::_);
    const std::vector<uint> communityVector = Rcpp::as<std::vector<uint>>(community);
    return communityVector;
}
