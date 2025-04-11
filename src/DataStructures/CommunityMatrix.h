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
    std::vector<uint32_t> GetCommunityMatrixByRow(int row) const;


    std::vector<std::vector<uint32_t>> GetRowAbundances() const {return rowAbundance;}
    const std::vector<std::vector<uint32_t>>& GetColumnEligibleIndexes() const {return eligibleRowIndexes;}
    std::vector<std::vector<uint32_t>>& GetAllIndexes() {return allIndexes;}

    Rcpp::CharacterVector GetSampleNames();
    const std::vector<uint32_t>& GetSums() const {return sums;}
    const int& GetRow() const {return row;}
    const int& GetColumn() const {return col;}
    const Rcpp::CharacterVector& GetRowNames() const {return rowNames;}
    const Rcpp::CharacterVector& GetColumnNames() const {return colNames;}


private:
    int row = 0;
    int col = 0;
    Rcpp::CharacterVector rowNames;
    Rcpp::CharacterVector colNames;
    Rcpp::NumericMatrix communityMatrix;
    std::vector<uint32_t> sums;
    std::vector<std::vector<uint32_t>> rowAbundance; // or abundances. The index represents the row number
    std::vector<std::vector<uint32_t>> eligibleRowIndexes; // The values in the row that do not have 0 abundance
    std::vector<std::vector<uint32_t>> allIndexes;



};



#endif //COMMUNITYMATRIX_H
