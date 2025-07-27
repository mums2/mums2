//
// Created by Gregory Johnson on 12/20/24.
//
#include <iostream>
#include <Rcpp.h>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include <RcppThread.h>
#include <mutex>
#include "Chemicals/MolecularFormula/MolecularFormula.h"
#include "CustomProgressBar/CliProgressBar.h"
#include "CustomProgressBar/ETAProgressBar.h"
#include "DataStructures/CommunityMatrix.h"
#include "DataStructures/CppMatrix.h"
#include "DiversityMetrics/Diversity.h"
#include "Rarefy/Rarefaction.h"
#include "DiversityMetrics/DiversityMetricFactory.h"
#include "FragmentationTree/FragmentationTree.h"
#include "FragmentationTree/GreedyHeuristic.h"
#include "Math/ParallelRandomNumberSitmo.h"
#include "Math/VectorMath.h"
#include "Spectra/ReadSpectra.h"

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


// Rcpp::NumericMatrix CalculateDiversityCommunityObject(const SEXP& communityMatrix, const std::string& diversityIndex) {
//     const Rcpp::XPtr<CommunityMatrix> matrix(communityMatrix);
//     const Rcpp::NumericMatrix abundances = matrix.get()->GetCommunityMatrix();
//     std::string index = diversityIndex;
//     std::transform(index.begin(), index.end(), index.begin(), tolower);
//     Diversity* diversity = DiversityMetricFactory::ChooseDiversityBasedOnIndex(index);
//     if(diversity == nullptr) {
//         Rcpp::stop("Diversity Metric not found");
//     }
//     Rcpp::NumericMatrix results = diversity->CalculateDiversity(abundances, index);
//     delete diversity;
//     return results;
// }

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

CppMatrix RarefactionCalculation2(const std::vector<std::vector<uint64_t>>& communityAbundances,
    const std::vector<std::vector<uint64_t>>& eligibleIndexes,
    const std::vector<std::vector<uint64_t>>& abundanceRanges,
    const std::vector<uint64_t>& sums,
    std::vector<ParallelRandomNumberSitmo>& rngEngines,
    const int rows,
    const int columns,
    const uint64_t size,
    const uint64_t threshold) {

    std::vector<double> resultantMatrix(rows * columns);
    for (int i = 0; i < rows; i++) {
        const std::vector<uint64_t> data = Rarefaction::Rarefy(communityAbundances[i], eligibleIndexes[i],
        abundanceRanges[i],rngEngines[i], size, sums[i], threshold);
        size_t counter = 0;
        for (int j = columns * i; j < columns * i + columns; j++) {
            resultantMatrix[j] = static_cast<double>(data[counter++]);
        }
    }
    return CppMatrix(resultantMatrix, rows, columns);
}

Rcpp::NumericMatrix TwoDimVectorToNumericMatrix(const std::vector<std::vector<uint64_t>>& matrix){
    Rcpp::NumericMatrix result(matrix.size(), matrix[0].size());
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix[i].size(); ++j) {
            result(i, j) = matrix[i][j];
        }
    }
    return result;
}

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
    // for(int i = 0; i < iterations; i++) {
    //     Rcpp::NumericMatrix rarefyMatrix = communityObject.get()->GetCommunityMatrix();
    //     if (subsample) {
    //         rarefyMatrix = RarefactionCalculation(communityObject,
    //             size, threshold, numberOfThreads, seed);
    //     }
    //     diversityMatrix += CalculateDiversity(rarefyMatrix, index);
    //     p.update(static_cast<float>(i)/static_cast<float>(iterations));
    // }
    // Rcpp::NumericMatrix rarefyMatrix = communityObject.get()->GetCommunityMatrix();
    const std::vector<std::vector<uint64_t>>& abundanceRanges = communityObject.get()->GetAbundanceRanges();
    const std::vector<std::vector<uint64_t>>& communityAbundances = communityObject.get()->GetCommunityAbundances();
    const std::vector<std::vector<uint64_t>>& eligibleIndexes = communityObject.get()->GetColumnEligibleIndexes();
    const std::vector<uint64_t>& sums = communityObject.get()->GetSums();
    std::vector<ParallelRandomNumberSitmo> rngEngines(matrixRowSize);
    for (int i = 0; i < matrixRowSize; ++i) {
        rngEngines[i] = ParallelRandomNumberSitmo(seed + i);
    }
    std::mutex mutex;
    int columnSize = communityObject.get()->GetColumn();
    int currentProgress = 0;
    const CppMatrix& matrix = communityObject.get()->GetCppMatrixOfAbundances();
    RcppThread::parallelFor(0, iterations, [&communityAbundances, &eligibleIndexes, &abundanceRanges,
        &diversityMatrix, &rngEngines, &matrix, &sums, &matrixRowSize, &columnSize, &size, &threshold, &subsample, &index,
        &iterations, &mutex, &p, &currentProgress](int i) {
        CppMatrix rarefyMatrix;
        if (subsample) {
            rarefyMatrix = RarefactionCalculation2(communityAbundances, eligibleIndexes,
                 abundanceRanges, sums, rngEngines, matrixRowSize, columnSize, size, threshold);
        }
        mutex.lock();
        if (subsample)
            diversityMatrix += CalculateDiversity(rarefyMatrix, index);
        // else
        //     diversityMatrix += CalculateDiversity(matrix, index);
        p.update(static_cast<float>(currentProgress++)/static_cast<float>(iterations));
        mutex.unlock();
    }, numberOfThreads,numberOfThreads);
    diversityMatrix/=iterations;
    p.end_display();

    Rcpp::NumericMatrix resultantMatrix = diversityMatrix.ToRcppMatrix();
    Rcpp::colnames(resultantMatrix) = samples;
    if(diversityMatrix.GetRowSize() <= 1) return resultantMatrix; // alpha diversity
    Rcpp::rownames(resultantMatrix) = samples;
    return resultantMatrix;
}

// [[Rcpp::export]]
Rcpp::List ReadMgf(const std::string& path) {
    ReadSpectra spectra;
    return spectra.ReadMGF(path);
}

// [[Rcpp::export]]
Rcpp::List ReadMsp(const std::string& path) {
    ReadSpectra spectra;
    return spectra.ReadMSP(path);
}

// [[Rcpp::export]]
Rcpp::NumericVector CompareMS2Ms1(const Rcpp::NumericVector& mz2, const Rcpp::NumericVector& mz1,
    const Rcpp::NumericVector& rt2, const Rcpp::NumericVector& rt1, const double mzThreshold,
    const double rtThreshold) {
    const auto currentSize = static_cast<size_t>(mz1.size());
    Rcpp::NumericVector resultsIndexes(currentSize, -1); // -1 means no match
    for (size_t i = 0; i < currentSize; i++) {
        double currentMz1 = mz1[i];
        double currentRt1 = rt1[i];
        Rcpp::NumericVector mzError = Rcpp::abs(currentMz1 - mz2) * 1e6 / currentMz1;
        Rcpp::NumericVector rtError = Rcpp::abs(currentRt1 - rt2);
        double bestDotProduct = 0;
        for (int j = 0; j < mzError.size(); j++) { // Pick score with the closest dotProduct value
            if (mzError[j] > mzThreshold || rtError[j] > rtThreshold) continue; // Over the threshold
            
            // Otherwise
            // Check if the similarity score (the dot product) is closer than the last one
            // If so replace
            double dotProduct = VectorMath::CosineScore({mz1[i], rt1[i]}, {mz2[j], rt2[j]});
            if (dotProduct < bestDotProduct) continue;
            resultsIndexes[i] = j + 1; // To match with R indexes add 1
            bestDotProduct = dotProduct;
            if (bestDotProduct >= 1) break; // If they are equal, we should break the loop and move on
        }
        
    }
    return resultsIndexes;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::export]]
std::string ComputeFragmentationTree(const Rcpp::List& molecularFormulas,
    const double parentMass, const int numberOfThreads) {
    const int size = molecularFormulas.size();
    FragmentationTree tree(molecularFormulas, parentMass);
    RcppThread::parallelFor(0, size, [&tree](int i) {
        tree.AddMolecularFormulaToGraph(i);
    }, numberOfThreads);
    return GreedyHeuristic::CalculateHeuristic(tree);
}

// [[Rcpp::export]]
SEXP CreateProgressBarObject() {
    auto* progressBar = new CliProgressBar();
    return Rcpp::XPtr<CliProgressBar>(progressBar);
}

// [[Rcpp::export]]
void IncrementProgressBar(SEXP& progressBar, const float progress) {
    const Rcpp::XPtr<CliProgressBar> cliProgressBar(progressBar);
    cliProgressBar.get()->update(progress);
}

// [[Rcpp::export]]
void DestroyProgressBar(SEXP& progressBar) {
    const Rcpp::XPtr<CliProgressBar> cliProgressBar(progressBar);
    cliProgressBar.get()->end_display();
}

// [[Rcpp::export]]
Rcpp::NumericMatrix Test() {
    std::vector<std::vector<double>> mat(3);
    mat[0] = {1,2,3};
    mat[1] = {4,5,6};
    mat[2] = {7,8,9};
    std::vector<std::vector<double>> mat2(3);
    mat2[0] = {1,2.2,3};
    mat2[1] = {4,5.2,6.1};
    mat2[2] = {7,8,9};
    const CppMatrix matrix(mat);
    const CppMatrix matrix2(mat2);
    const auto res = matrix + matrix2;
    std::vector<std::vector<double>> expectedResult(3);
    expectedResult[0] = {2,4,6};
    expectedResult[1] = {8,10,12};
    expectedResult[2] = {14,16,18};
    const CppMatrix expectedMatrix(expectedResult);
    return matrix2.ToRcppMatrix();
}

// [[Rcpp::export]]
void TestNumericMatrix(Rcpp::NumericVector& vector) {
    vector.attr("dim") = Rcpp::Dimension(100, 100);
    Rcpp::NumericMatrix m = Rcpp::as<Rcpp::NumericMatrix>(vector);
    Rcpp::NumericMatrix m2 = Rcpp::as<Rcpp::NumericMatrix>(vector);
    m + m2;

}

// [[Rcpp::export]]

void TestCppMatrix() {
    std::vector<double> mat(std::vector<double>(10000, 1));
    std::vector<double> mat2(std::vector<double>(10000, 1));
    const CppMatrix matrix(mat, 100, 100);
    const CppMatrix matrix2(mat2, 100, 100);
    const auto res = matrix + matrix2;
}

// [[Rcpp::export]]
void TestTraditional() {
    std::vector<std::vector<double>> mat(100, std::vector<double>(100, 1));
    std::vector<std::vector<double>> mat2(100, std::vector<double>(100, 1));
    std::vector<std::vector<double>> result(100, std::vector<double>(100, 0));
    for (int i = 0 ; i < mat.size() ; i++) {
        for (int j = 0 ; j < mat[0].size() ; j++) {
            result[i][j] = mat[i][j] + mat2[i][j];
        }
    }
}