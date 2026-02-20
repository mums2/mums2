//
// Created by Gregory Johnson on 11/19/25.
//
#include <RcppThread.h>
#include <Rcpp.h>
#include <string>
#include <vector>
#include <mutex>


#include "CustomProgressBar/CliProgressBar.h"
#include "DataStructures/CommunityMatrix.h"
#include "DataStructures/CppMatrix.h"
#include "DiversityMetrics/Diversity.h"
#include "DiversityMetrics/DiversityMetricFactory.h"
#include "Math/ParallelRandomNumberSitmo.h"
#include "Rarefy/Rarefaction.h"


CppMatrix CalculateDiversity(const CppMatrix& abundances, const std::string& diversityIndex) {
    std::string index = diversityIndex;
    std::transform(index.begin(), index.end(), index.begin(), tolower);
    Diversity* diversity = DiversityMetricFactory::ChooseDiversityBasedOnIndex(index);
    if(diversity == nullptr) {
        Rcpp::stop("Diversity Metric not found");
    }
    CppMatrix results = diversity->CalculateDiversity(abundances, index);
    delete diversity;
    return results;
}

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
Rcpp::NumericMatrix RarefactionCalculation(const SEXP& communityMatrix, const uint64_t size,
    const uint64_t threshold, const int numberOfThreads, const int seed = 123) {

    const Rcpp::XPtr<CommunityMatrix> matrix(communityMatrix);
    const int row = matrix.get()->GetRow();
    const int col = matrix.get()->GetColumn();
    const Rcpp::CharacterVector& rowNames = matrix.get()->GetRowNames();
    const Rcpp::CharacterVector& columnNames = matrix.get()->GetColumnNames();
    std::vector<std::string> names = Rcpp::as<std::vector<std::string> >(rowNames);
    const std::vector<std::vector<uint64_t>>& abundanceRanges = matrix.get()->GetAbundanceRanges();
    const std::vector<std::vector<uint64_t>>& communityAbundances = matrix.get()->GetCommunityAbundances();
    const std::vector<std::vector<uint64_t>>& eligibleIndexes = matrix.get()->GetColumnEligibleIndexes();
    const std::vector<uint64_t>& sums = matrix.get()->GetSums();
    Rcpp::NumericMatrix resultantMatrix(row, col);
    std::vector<ParallelRandomNumberSitmo> rngEngines(row);
    for (int i = 0; i < row; ++i) {
        rngEngines[i] = ParallelRandomNumberSitmo(seed + i);
    }

    std::mutex mutex;
    RcppThread::parallelFor(0, row, [&communityAbundances, &eligibleIndexes,
    &abundanceRanges, &resultantMatrix, &sums, &rngEngines, &mutex, &size, &threshold](int i) {
        const std::vector<uint64_t> result = Rarefaction::Rarefy(communityAbundances[i], eligibleIndexes[i],
          abundanceRanges[i],rngEngines[i], size, sums[i], threshold);
        mutex.lock();
        for(const auto& index : eligibleIndexes[i]) {
          resultantMatrix(i, index) = result.at(index);
        }
        mutex.unlock();
    }, numberOfThreads);
    Rcpp::rownames(resultantMatrix) = rowNames;
    Rcpp::colnames(resultantMatrix) = columnNames;
    return resultantMatrix;
}

CppMatrix RarefactionCalculationParallelized(const std::vector<std::vector<uint64_t>>& communityAbundances,
    const std::vector<std::vector<uint64_t>>& eligibleIndexes,
    const std::vector<std::vector<uint64_t>>& abundanceRanges,
    const std::vector<uint64_t>& sums,
    ParallelRandomNumberSitmo& rngEngine,
    const int rows,
    const int columns,
    const uint64_t size,
    const uint64_t threshold) {

    std::vector<double> resultantMatrix(rows * columns);
    for (int i = 0; i < rows; i++) {
        const std::vector<uint64_t> data = Rarefaction::Rarefy(communityAbundances[i], eligibleIndexes[i],
        abundanceRanges[i],rngEngine, size, sums[i], threshold);
        size_t counter = 0;
        for (int j = columns * i; j < columns * i + columns; j++) {
            resultantMatrix[j] = static_cast<double>(data[counter++]);
        }
    }
    return CppMatrix(resultantMatrix, rows, columns);
}
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix FasterAvgDist(const SEXP& communityMatrix, const std::string& index,
    const uint64_t size, const uint64_t threshold, const bool subsample, const int numberOfThreads,
    const int iterations = 1000, const int seed = 123) {
    CliProgressBar p;
    const Rcpp::XPtr<CommunityMatrix> communityObject(communityMatrix);
    const Rcpp::CharacterVector samples = communityObject.get()->GetSampleNames();
    const size_t matrixRowSize = samples.size();
    int row = samples.size();
    if (index == "simpson" || index == "shannon")
        row = 1;
    CppMatrix diversityMatrix(std::vector<double>(row * samples.size(), 0), row, samples.size());
    const std::vector<std::vector<uint64_t>>& abundanceRanges = communityObject.get()->GetAbundanceRanges();
    const std::vector<std::vector<uint64_t>>& communityAbundances = communityObject.get()->GetCommunityAbundances();
    const std::vector<std::vector<uint64_t>>& eligibleIndexes = communityObject.get()->GetColumnEligibleIndexes();
    const std::vector<uint64_t>& sums = communityObject.get()->GetSums();
    std::vector<ParallelRandomNumberSitmo> rngEngines(iterations);
    for (int i = 0; i < iterations; ++i) {
        rngEngines[i] = ParallelRandomNumberSitmo(seed + i);
    }
    std::mutex mutex;
    int columnSize = communityObject.get()->GetColumn();
    int currentProgress = 0;
    const CppMatrix& matrix = communityObject.get()->GetCppMatrixOfAbundances();
    RcppThread::parallelFor(0, iterations, [&communityAbundances, &eligibleIndexes, &abundanceRanges,
        &diversityMatrix, &rngEngines, &matrix, &sums, &matrixRowSize, &columnSize, &size, &threshold,
        &subsample, &index, &iterations, &mutex, &p, &currentProgress](int i) {
        CppMatrix rarefyMatrix;

        if (subsample) {
            rarefyMatrix = RarefactionCalculationParallelized(communityAbundances, eligibleIndexes,
                 abundanceRanges, sums, rngEngines[i], matrixRowSize, columnSize, size, threshold);
        }
        CppMatrix diversityCppMatrix;
        if (subsample)
            diversityCppMatrix = CalculateDiversity(rarefyMatrix, index);
        else
            diversityCppMatrix = CalculateDiversity(matrix, index);
        mutex.lock();
        diversityMatrix += diversityCppMatrix;
        p.update(static_cast<float>(currentProgress++)/static_cast<float>(iterations));
        mutex.unlock();
    }, numberOfThreads);
    diversityMatrix/=iterations;
    p.end_display();

    Rcpp::NumericMatrix resultantMatrix = diversityMatrix.ToRcppMatrix();
    Rcpp::colnames(resultantMatrix) = samples;
    if(diversityMatrix.GetRowSize() <= 1) return resultantMatrix; // alpha diversity
    Rcpp::rownames(resultantMatrix) = samples;
    return resultantMatrix;
}
