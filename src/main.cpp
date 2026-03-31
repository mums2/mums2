//
// Created by Gregory Johnson on 12/20/24.
//

#include <Rcpp.h>
#include <RcppThread.h>
#include <string>
#include "FragmentationTree/FragmentationTree.h"
#include "FragmentationTree/GreedyHeuristic.h"
#include "Math/VectorMath.h"

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
    FragmentationTree tree(molecularFormulas, parentMass);
    const std::vector<int>& availableIndexes = tree.GetColorZeroFormulas();
    const int size = static_cast<int>(availableIndexes.size());
    RcppThread::parallelFor(0, size, [&tree, &availableIndexes](int i) {
        tree.AddMolecularFormulaToGraph(availableIndexes[i]);
    }, numberOfThreads);
    return GreedyHeuristic::CalculateHeuristic(tree);
}

#include "../../../../Downloads/gperftools-2.15/src/gperftools/profiler.h"

// [[Rcpp::export]]
SEXP start_profiler(const SEXP& str) {
    ProfilerStart(Rcpp::as<const char*>(str));
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP stop_profiler() {
    ProfilerStop();
    return R_NilValue;
}