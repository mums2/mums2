//
// Created by Gregory Johnson on 12/20/24.
//
#include <iostream>
#include <Rcpp.h>
#include <cstdint>
#include <numeric>
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
        // This will return a vector of the size of eligibileIndexes, so when we assemble it
        // We can use that to our advantage
        const auto results = rarefaction.Rarefy2(communityVector, eligibleIndexes, eligibleAbundances,
            size, threshold);
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
std::vector<size_t> UpdatedValue(const std::vector<int64_t> &weightRanges,
    const int64_t sizeToPull, const int64_t sum) {
    std::set<CountIndexPair> vals;
    for(int i = 0; i < weightRanges.size(); i++) {
        vals.insert(CountIndexPair{i,weightRanges[i]});
    }
    std::vector<size_t> indexes(sizeToPull);
    for(int i = 0; i < sizeToPull; i++) {
        const CountIndexPair randomNum{-1,
            static_cast<int64_t>(R::runif(0, static_cast<double>(sum - i)))};
        // const auto randomNum = static_cast<int64_t>(R::runif(0, static_cast<double>(sum - i)));
        const auto val = vals.upper_bound(randomNum);
        const CountIndexPair updatedPair{val->index , val->abundance - 1};
        // std::cout << val->index << " " << val->abundance << std::endl;
        indexes[i] = val->index;
        vals.erase(val);
        vals.insert(updatedPair);
    }
    return indexes;
}

