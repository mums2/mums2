//
// Created by gregj on 2/18/2025.
//

#ifndef COMMUNITYMATRIX_H
#define COMMUNITYMATRIX_H
#include <cstdint>
#include <Rcpp.h>
#include <vector>

class CommunityMatrix {
public:
    explicit CommunityMatrix(const Rcpp::NumericMatrix& matrix);
    ~CommunityMatrix() = default;
    void InitializeMatrix();
    const Rcpp::NumericMatrix& GetCommunityMatrix() const {return communityMatrix;}
    const std::vector<int64_t>& GetCommunityMatrixByRow(int row) const;
    std::vector<std::vector<int64_t>> GetRowAbundances() const {return rowAbundance;}
    const std::vector<std::vector<int64_t>>& GetColumnEligibleIndexes() const {return eligibleRowIndexes;}
    const std::vector<int64_t>& GetSums() const {return sums;}
    const int& GetRow() const {return row;}
    const int& GetColumn() const {return col;}

private:
    int row;
    int col;
    Rcpp::NumericMatrix communityMatrix;
    std::vector<int64_t> sums;
    std::vector<std::vector<int64_t>> rowAbundance; // or abundances. The index represents the row number
    std::vector<std::vector<int64_t>> eligibleRowIndexes; // The values in the row that do not have 0 abundance



};



#endif //COMMUNITYMATRIX_H
