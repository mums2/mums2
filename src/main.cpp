//
// Created by Gregory Johnson on 12/20/24.
//
#include <iostream>
#include <Rcpp.h>
#include <cstdint>
#include <numeric>
#include <algorithm>

#include "DataStructures/CommunityMatrix.h"
#include "DiversityMetrics/Diversity.h"
#include "Rarefy/Rarefaction.h"
#include "DiversityMetrics/DiversityMetricFactory.h"
#include "Random/RandomizationMethods.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix CalculateDiversity(const Rcpp::NumericMatrix& abundances, const std::string& diversityIndex) {
    std::string index = diversityIndex;
    std::transform(index.begin(), index.end(), index.begin(), tolower);
    Diversity* diversity = DiversityMetricFactory::ChooseDiversityBasedOnIndex(index);
    if(diversity == nullptr) {
        Rcpp::stop("Diversity Metric not found");
    }
    Rcpp::NumericMatrix results = diversity->CalculateDiversity(abundances, index);
    delete diversity;
    return results;
}


/* TODO
    We are going to pay the price of generating a vector of a size of 1billion.
    Do that when we create the community matrix.
    This vector will have a size of 1 -> the sum of the LARGEST row.
    Therefore, for rows that are smaller than it, it may miss the values a few times.
    So when we past the values to a function, we have to create an offset map. So we know
    Which value was chosen. This allows us to not have to copy the matrix over (saving time and space)
*/


// [[Rcpp::export]]
Rcpp::NumericMatrix RarefactionCalculation(const Rcpp::NumericMatrix& communityMatrix, const int64_t size,
    const int64_t threshold) {
    const int row = communityMatrix.nrow();
    const int col = communityMatrix.ncol();
    const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);
    std::vector<int> indexToName(col);
    std::iota(indexToName.begin(), indexToName.end(), 0);
    Rarefaction rarefaction;
    Rcpp::NumericMatrix resultantMatrix(row, col);
    for(int i = 0; i < row; i++) {
        Rcpp::NumericVector community = communityMatrix(i, Rcpp::_);
        std::vector<int64_t> communityVector = Rcpp::as<std::vector<int64_t>>(community);
        size_t communityVectorSize = communityVector.size();
        std::vector<int64_t> eligibleIndexes;
        std::vector<int64_t> eligibleAbundances;
        std::vector<int64_t> abundanceRanges(1,0);
       // abundanceRanges.reserve(communityVectorSize);
        eligibleAbundances.reserve(communityVectorSize);
        eligibleIndexes.reserve(communityVectorSize);
        for(size_t j = 0; j < communityVectorSize; j++) {
            int64_t val = communityVector[j];
            if(val > 0) {
                eligibleIndexes.emplace_back(j);
                eligibleAbundances.emplace_back(val);
               // abundanceRanges.emplace_back(abundanceRanges[count++ - 1] + val);
            }
        }
        // Rcpp::Rcout << communityVector << std::endl;
        // Rarefy
        // We are going to have to switch between the transpose of the community matrix and
        // the original matrix.
        // Samples are represented by rows and columns represent species.
        // So ex..
        //          FD39  FD09 <- is backwards            species1     <- Correct
        // Species1  0      1                       FD39    1
        //                                          DF09    0

        // const auto results = rarefaction.Rarefy(indexToName, communityVector,
        //
        //                                         eligibleIndexes, eligibleAbundances, abundanceRanges, size, threshold);
        // This will return a vector of the size of eligibile Indexes, so when we assemble it
        // We can use that to our advantage
        const int64_t sum = std::accumulate(eligibleAbundances.begin(), eligibleAbundances.end(), 0LL);
        const auto results = rarefaction.Rarefy2(communityVector, eligibleIndexes, eligibleAbundances,
                                                 size, sum, threshold);
        for(const auto& index : eligibleIndexes) {
            resultantMatrix(i, index) = results[index];
        }
        // for(int j = 0; j < col; j++) {
        //     resultantMatrix(i, j) = results[j];
        // }
    }
    Rcpp::rownames(resultantMatrix) = samples;
    return resultantMatrix;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix FasterAvgDist(const Rcpp::NumericMatrix& communityMatrix, const std::string& index,
    const int64_t size, const int64_t threshold, const int iterations = 1000) {
    const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);

    Rcpp::NumericMatrix diversity = CalculateDiversity(RarefactionCalculation(communityMatrix,
        size, threshold), index);
    for(int i = 1; i < iterations; i++) {
        Rcpp::NumericMatrix rarefyMatrix = RarefactionCalculation(communityMatrix, size, threshold);
        diversity += CalculateDiversity(rarefyMatrix, index);
    }
    diversity = diversity/iterations;
    Rcpp::colnames(diversity) = samples;
    Rcpp::rownames(diversity) = samples;

    return diversity;
}

struct CountIndexPair {
    int64_t index;
    int64_t abundance;
    bool operator<(const CountIndexPair& other) const {
        return abundance < other.abundance;
    }
};

// [[Rcpp::export]]
SEXP CreateCommunityMatrix(Rcpp::NumericMatrix communityMatrix) {
    auto* matrix = new CommunityMatrix(communityMatrix);
    matrix->InitializeMatrix();
    return Rcpp::XPtr<CommunityMatrix>(matrix);
}

// [[Rcpp::export]]
SEXP GetCommunityMatrix(SEXP communityMatrix) {
    const Rcpp::XPtr<CommunityMatrix> matrix(communityMatrix);
    return matrix.get()->GetCommunityMatrix();
}





// [[Rcpp::export]]
Rcpp::NumericMatrix RarefactionCalculation2(const SEXP& communityMatrix, const int64_t size,
    const int64_t threshold) {
    const Rcpp::XPtr<CommunityMatrix> matrix(communityMatrix);
    const int row = matrix.get()->GetRow();
    const int col = matrix.get()->GetColumn();
    const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);
    Rarefaction rarefaction;
    const std::vector<std::vector<int64_t>>& eligibleIndexes = matrix.get()->GetColumnEligibleIndexes();
    std::vector<std::vector<int64_t>> eligibleAbundances = matrix.get()->GetRowAbundances();
    const std::vector<int64_t>& sums = matrix.get()->GetSums();
    Rcpp::NumericMatrix resultantMatrix(row, col);
    for(int i = 0; i < row; i++) {
        const std::vector<int64_t>& communityVector = matrix.get()->GetCommunityMatrixByRow(i);
        const auto results = rarefaction.Rarefy2(communityVector, eligibleIndexes[i], eligibleAbundances[i],
                                                 size, sums[i], threshold);

        for(const auto& index : eligibleIndexes[i]) {
            resultantMatrix(i, index) = results.at(index);
        }

    }
    Rcpp::rownames(resultantMatrix) = samples;
    return resultantMatrix;
}
