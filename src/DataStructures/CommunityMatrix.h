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
    std::vector<uint> GetCommunityMatrixByRow(int row) const;
    std::vector<std::vector<uint>> GetRowAbundances() const {return rowAbundance;}
    const std::vector<std::vector<uint>>& GetColumnEligibleIndexes() const {return eligibleRowIndexes;}
    std::vector<std::vector<uint>>& GetAllIndexes() {return allIndexes;}

    const std::vector<uint>& GetSums() const {return sums;}
    const int& GetRow() const {return row;}
    const int& GetColumn() const {return col;}

private:
    int row = 0;
    int col = 0;
    Rcpp::NumericMatrix communityMatrix;
    std::vector<uint> sums;
    std::vector<std::vector<uint>> rowAbundance; // or abundances. The index represents the row number
    std::vector<std::vector<uint>> eligibleRowIndexes; // The values in the row that do not have 0 abundance
    std::vector<std::vector<uint>> allIndexes;



};



#endif //COMMUNITYMATRIX_H
