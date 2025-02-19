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
    explicit CommunityMatrix(const NumericMatrix& matrix);
    ~CommunityMatrix();
    void InitializeMatrix();

private:
    NumericMatrix communityMatrix{};
    std::vector<std::vector<int64_t>> rowAbundance; // or abundances. The index represents the row number
    std::vector<std::vector<int64_t>> eligibleRowIndexes; // The values in the row that do not have 0 abundance

};



#endif //COMMUNITYMATRIX_H
