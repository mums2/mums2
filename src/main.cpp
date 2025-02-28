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

// [[Rcpp::export]]
Rcpp::NumericMatrix CalculateDiversityCommunityObject(const SEXP& communityMatrix, const std::string& diversityIndex) {
    const Rcpp::XPtr<CommunityMatrix> matrix(communityMatrix);
    const Rcpp::NumericMatrix abundances = matrix.get()->GetCommunityMatrix();
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
        // const auto results = rarefaction.Rarefy2(communityVector, eligibleIndexes, TODO,
        //                                          TODO, size, sum, threshold);
        // for(const auto& index : eligibleIndexes) {
        //     resultantMatrix(i, index) = results[index];
        // }
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
void ShuffleVectorNoConversion(Rcpp::NumericVector& vec) {
    const int indexSize = static_cast<int>(vec.size());
    vec = Rcpp::sample(vec, indexSize);
}
// [[Rcpp::export]]
void ShuffleVectorConversion(const std::vector<int64_t>& vec) {
    const Rcpp::NumericVector indexes = Rcpp::wrap(vec);
    const int indexSize = static_cast<int>(vec.size());
    Rcpp::NumericVector shuffledValues = Rcpp::sample(indexes, indexSize);
}


// [[Rcpp::export]]
std::vector<int64_t> ShuffleWithRandomNumbers(std::vector<int64_t>& vec) {
    const int size = static_cast<int>(vec.size());
    std::vector<int64_t> shuffledVector(size, 0);
    for(int i = 0; i < vec.size(); i++) {
        const auto value = static_cast<int64_t>(R::runif(i, size));
        shuffledVector[i] = vec[value];
        const auto currentVal = vec[i];
        vec[i] = vec[value];
        vec[value] = currentVal;
    }
    return shuffledVector;
}



void ShuffleVector(std::vector<int64_t>& vec) {
    const Rcpp::NumericVector indexes = Rcpp::wrap(vec);
    const int indexSize = static_cast<int>(vec.size());
    Rcpp::NumericVector shuffledValues = Rcpp::sample(indexes, indexSize);
    vec = Rcpp::as<std::vector<int64_t>>(shuffledValues);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RarefactionCalculation2(const SEXP& communityMatrix, const int64_t size,
    const int64_t threshold) {

    const auto startTime = std::chrono::steady_clock::now();
    const Rcpp::XPtr<CommunityMatrix> matrix(communityMatrix);
    const int row = matrix.get()->GetRow();
    const int col = matrix.get()->GetColumn();

    const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);
    Rarefaction rarefaction;
    std::vector<std::vector<int64_t>>& allIndexes = matrix.get()->GetAllIndexes();
    //ShuffleVectorNoConversion(allIndexes);
    const std::vector<std::vector<int64_t>>& abundanceRanges  = matrix.get()->GetAbundanceRanges();
    const std::vector<std::vector<int64_t>>& eligibleIndexes = matrix.get()->GetColumnEligibleIndexes();
    std::vector<std::vector<int64_t>> eligibleAbundances = matrix.get()->GetRowAbundances();
    const std::vector<int64_t>& sums = matrix.get()->GetSums();
    const auto endTime = std::chrono::steady_clock::now();
    const auto time = std::chrono::duration_cast<std::chrono::milliseconds>(endTime-startTime).count();
    Rcpp::Rcout << "Initialization Time: " << time << std::endl;
    Rcpp::NumericMatrix resultantMatrix(row, col);
    for(int i = 0; i < row; i++) {
        const std::vector<int64_t> communityVector = matrix.get()->GetCommunityMatrixByRow(i);
        const auto results = rarefaction.Rarefy2(communityVector, eligibleIndexes[i], allIndexes[i],
                                                 abundanceRanges[i], size, sums[i], threshold);

        for(const auto& index : eligibleIndexes[i]) {
            resultantMatrix(i, index) = results.at(index);
        }

    }
    Rcpp::rownames(resultantMatrix) = samples;
    return resultantMatrix;
}

int randWrapper(const int n) { return floor(R::unif_rand()*n); }
// [[Rcpp::export]]
void randomShuffle2(Rcpp::NumericVector& a) {
    // clone a into b to leave a alone
    int n = a.size();
    int j;

    // Fisher-Yates Shuffle Algorithm
    for (int i = 0; i < n - 1; i++) {
        j = i + randWrapper(n - i);
        std::swap(a[i], a[j]);
    }
}

// [[Rcpp::export]]
Rcpp::NumericMatrix RarefactionCalculationFisherYates(const SEXP& communityMatrix, const int64_t size,
    const int64_t threshold) {
    //const auto startTime = std::chrono::steady_clock::now();
    const Rcpp::XPtr<CommunityMatrix> matrix(communityMatrix);
    const int row = matrix.get()->GetRow();
    const int col = matrix.get()->GetColumn();

    const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);
    Rarefaction rarefaction;
    std::vector<std::vector<int64_t>>& allIndexes = matrix.get()->GetAllIndexes();
    // const auto shuffleStartTime = std::chrono::steady_clock::now();
    // //randomShuffle2(allIndexes);
    const std::vector<std::vector<int64_t>>& abundanceRanges  = matrix.get()->GetAbundanceRanges();
    const std::vector<std::vector<int64_t>>& eligibleIndexes = matrix.get()->GetColumnEligibleIndexes();
    std::vector<std::vector<int64_t>> eligibleAbundances = matrix.get()->GetRowAbundances();
    const std::vector<int64_t>& sums = matrix.get()->GetSums();
    Rcpp::NumericMatrix resultantMatrix(row, col);
    // const auto endTime = std::chrono::steady_clock::now();
    // const auto time = std::chrono::duration_cast<std::chrono::milliseconds>(endTime-startTime).count();
    // Rcpp::Rcout << "Initialization Time: " << time << " Milliseconds " << std::endl;
    for(int i = 0; i < row; i++) {
        const std::vector<int64_t> communityVector = matrix.get()->GetCommunityMatrixByRow(i);
        const auto results = rarefaction.Rarefy2(communityVector, eligibleIndexes[i], allIndexes[i],
                                                 abundanceRanges[i], size, sums[i], threshold);

        for(const auto& index : eligibleIndexes[i]) {
            resultantMatrix(i, index) = results.at(index);
        }

    }
    Rcpp::rownames(resultantMatrix) = samples;
    return resultantMatrix;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix FasterAvgDist2(const SEXP& communityMatrix, const std::string& index,
    const int64_t size, const int64_t threshold, const int iterations = 1000) {
    const Rcpp::CharacterVector samples = Rcpp::rownames(communityMatrix);
    const Rcpp::XPtr<CommunityMatrix> communityObject(communityMatrix);

    Rcpp::NumericMatrix diversity = CalculateDiversity(RarefactionCalculationFisherYates(communityObject,
        size, threshold), index);
    for(int i = 1; i < iterations; i++) {
        Rcpp::NumericMatrix rarefyMatrix = RarefactionCalculationFisherYates(communityObject,
            size, threshold);
        diversity += CalculateDiversity(rarefyMatrix, index);
    }
    diversity = diversity/iterations;
    Rcpp::colnames(diversity) = samples;
    Rcpp::rownames(diversity) = samples;

    return diversity;
}



