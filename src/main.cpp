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
#include "DiversityMetrics/Diversity.h"
#include "Rarefy/Rarefaction.h"
#include "DiversityMetrics/DiversityMetricFactory.h"
#include "FragmentationTree/FragmentationTree.h"
#include "FragmentationTree/GreedyHeuristic.h"
#include "Math/ParallelRandomNumberSitmo.h"
#include "Math/VectorMath.h"
#include "Spectra/ReadSpectra.h"

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

// [[Rcpp::depends(sitmo)]]
// [[Rcpp::export]]
Rcpp::NumericMatrix RarefactionCalculation(const SEXP& communityMatrix, const uint32_t size,
    const uint32_t threshold, const int seed = 123) {

    const Rcpp::XPtr<CommunityMatrix> matrix(communityMatrix);
    const int row = matrix.get()->GetRow();
    const int col = matrix.get()->GetColumn();
    const Rcpp::CharacterVector& rowNames = matrix.get()->GetRowNames();
    const Rcpp::CharacterVector& columnNames = matrix.get()->GetColumnNames();
    std::vector<std::string> names = Rcpp::as<std::vector<std::string> >(rowNames);
    const std::vector<std::vector<uint32_t>>& abundanceRanges = matrix.get()->GetAbundanceRanges();
    const std::vector<std::vector<uint32_t>>& communityAbundances = matrix.get()->GetCommunityAbundances();
    const std::vector<std::vector<uint32_t>>& eligibleIndexes = matrix.get()->GetColumnEligibleIndexes();
    const std::vector<uint32_t>& sums = matrix.get()->GetSums();
    Rcpp::NumericMatrix resultantMatrix(row, col);
    std::vector<sitmo::prng> rngEngines(row);
    for (int i = 0; i < row; ++i) {
        rngEngines[i] = ParallelRandomNumberSitmo::CreateRandomNumberGenerateSitmo(seed + i);
    }

    std::mutex mutex;
    RcppThread::parallelFor(0, row, [&communityAbundances, &eligibleIndexes,
    &abundanceRanges, &resultantMatrix, &sums, &rngEngines, &mutex, &size, &threshold](int i) {
        const std::vector<uint32_t> result = Rarefaction::Rarefy(communityAbundances[i], eligibleIndexes[i],
          abundanceRanges[i],rngEngines[i], size, sums[i], threshold);
        mutex.lock();
        for(const auto& index : eligibleIndexes[i]) {
          resultantMatrix(i, index) = result.at(index);
        }
        mutex.unlock();
    }, 16);
    Rcpp::rownames(resultantMatrix) = rowNames;
    Rcpp::colnames(resultantMatrix) = columnNames;
    return resultantMatrix;
}



// [[Rcpp::export]]
Rcpp::NumericMatrix FasterAvgDist(const SEXP& communityMatrix, const std::string& index,
    const uint32_t size, const uint32_t threshold, const bool subsample,
    const int iterations = 1000, const int seed = 123) {
    CliProgressBar p;
    const Rcpp::XPtr<CommunityMatrix> communityObject(communityMatrix);
    const Rcpp::CharacterVector samples = communityObject.get()->GetSampleNames();
    int row = samples.size();
    if (index == "simpson" || index == "shannon")
        row = 1;
    Rcpp::NumericMatrix diversityMatrix(row, samples.size());
    for(int i = 0; i < iterations; i++) {
        Rcpp::NumericMatrix rarefyMatrix = communityObject.get()->GetCommunityMatrix();
        if (subsample) {
            rarefyMatrix = RarefactionCalculation(communityObject,
                size, threshold);
        }
        diversityMatrix += CalculateDiversity(rarefyMatrix, index);
        p.update(static_cast<float>(i)/static_cast<float>(iterations));
    }
    diversityMatrix = diversityMatrix/iterations;
    Rcpp::colnames(diversityMatrix) = samples;
    if(diversityMatrix.rows() <= 1) return diversityMatrix; // alpha diversity
    Rcpp::rownames(diversityMatrix) = samples;
    p.end_display();
    return diversityMatrix;
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



// [[Rcpp::depends(sitmo)]]
#include <sitmo.h>
// [[Rcpp::export]]
Rcpp::NumericVector runif_sitmo(unsigned int n, double min = 0.0, double max = 1.0, uint32_t seed = 1) {
    Rcpp::NumericVector o(n);

    // Create a prng engine
    sitmo::prng eng(seed);
    // Obtain the range between max and min
    double dis = max - min;

    for(int i = 0; i < n; ++i) {
        // Sample from the RNG and divide it by the maximum value possible (can also use SITMO_RAND_MAX, which is 4294967295)
        // Apply appropriate scale (MAX-MIN)
        o[i] = min + (static_cast<double>(eng()) / (sitmo::prng::max())) * (dis);
    }

    return o;
}

// [[Rcpp::export]]
int Test(int seed) {
    return 0;
}