//
// Created by Gregory Johnson on 12/20/24.
//
#include <iostream>
#include <Rcpp.h>
#include <cstdint>
#include <numeric>
#include <algorithm>

#include "Chemicals/MolecularFormula/MolecularFormula.h"
#include "DataStructures/CommunityMatrix.h"
#include "DiversityMetrics/Diversity.h"
#include "Rarefy/Rarefaction.h"
#include "DiversityMetrics/DiversityMetricFactory.h"
#include "FragmentationTree/FragmentationTree.h"
#include "FragmentationTree/GreedyHeuristic.h"
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

// [[Rcpp::export]]
Rcpp::NumericMatrix RarefactionCalculation(const SEXP& communityMatrix, const uint32_t size,
    const uint32_t threshold) {

    const Rcpp::XPtr<CommunityMatrix> matrix(communityMatrix);
    const int row = matrix.get()->GetRow();
    const int col = matrix.get()->GetColumn();
    const Rcpp::CharacterVector& rowNames = matrix.get()->GetRowNames();
    const Rcpp::CharacterVector& columnNames = matrix.get()->GetColumnNames();

    Rarefaction rarefaction;
    std::vector<std::vector<uint32_t>>& allIndexes = matrix.get()->GetAllIndexes();
    const std::vector<std::vector<uint32_t>>& eligibleIndexes = matrix.get()->GetColumnEligibleIndexes();
    std::vector<std::vector<uint32_t>> eligibleAbundances = matrix.get()->GetRowAbundances();
    const std::vector<uint32_t>& sums = matrix.get()->GetSums();
    Rcpp::NumericMatrix resultantMatrix(row, col);
    for(int i = 0; i < row; i++) {
        const std::vector<uint32_t> communityVector = matrix.get()->GetCommunityMatrixByRow(i);
        const auto results = rarefaction.Rarefy(communityVector, eligibleIndexes[i], allIndexes[i],
                                                 size, sums[i], threshold);

        for(const auto& index : eligibleIndexes[i]) {
            resultantMatrix(i, index) = results.at(index);
        }

    }
    Rcpp::rownames(resultantMatrix) = rowNames;
    Rcpp::colnames(resultantMatrix) = columnNames;
    return resultantMatrix;
}

// [[Rcpp::export]]
Rcpp::DataFrame FasterAvgDist(const SEXP& communityMatrix, const std::string& index,
    const uint32_t size, const uint32_t threshold, const int iterations = 1000) {
    const Rcpp::XPtr<CommunityMatrix> communityObject(communityMatrix);
    const Rcpp::CharacterVector samples = communityObject.get()->GetSampleNames();

    Rcpp::NumericMatrix diversityMatrix = CalculateDiversity(RarefactionCalculation(communityObject,
        size, threshold), index);
    for(int i = 1; i < iterations; i++) {
        Rcpp::NumericMatrix rarefyMatrix = RarefactionCalculation(communityObject,
            size, threshold);
        diversityMatrix += CalculateDiversity(rarefyMatrix, index);
    }
    diversityMatrix = diversityMatrix/iterations;
    Rcpp::colnames(diversityMatrix) = samples;
    if(diversityMatrix.rows() <= 1) return diversityMatrix; // alpha diversity
    Rcpp::rownames(diversityMatrix) = samples;
    const int sampleSize = std::pow(samples.size(), 2);
    Rcpp::CharacterVector firstSample(sampleSize);
    Rcpp::CharacterVector otherSample(sampleSize);
    Rcpp::NumericVector diversityResults(sampleSize);
    int currentIndex = 0;
    for(int i = 0; i < samples.size(); i++) {
        for(int j = 0; j < samples.size(); j++) {
            firstSample[currentIndex] = samples[i];
            otherSample[currentIndex] = samples[j];
            diversityResults[currentIndex++] = diversityMatrix(i, j);
        }
    }

    return Rcpp::DataFrame::create(Rcpp::Named("firstSample") = firstSample,
        Rcpp::Named("otherSample") = otherSample,
        Rcpp::Named("diversity") = diversityResults);
}

// [[Rcpp::export]]
Rcpp::List ReadMgf(const std::string& path) {
    ReadSpectra spectra;
    return(spectra.ReadMGF(path));
}

// [[Rcpp::export]]
Rcpp::List ReadMsp(const std::string& path) {
    ReadSpectra spectra;
    return(spectra.ReadMSP(path));
}

double DotProduct(Rcpp::NumericVector x, Rcpp::NumericVector y) {
    double dotValue = Rcpp::sum(x * y);
    double magnitudeOne = std::sqrt(Rcpp::sum(Rcpp::pow(x, 2)));
    double magnitudeTwo = std::sqrt(Rcpp::sum(Rcpp::pow(y, 2)));
    return dotValue / (magnitudeOne * magnitudeTwo);
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
        Rcpp::NumericVector mzError = Rcpp::abs(currentMz1 - mz2);
        Rcpp::NumericVector rtError = Rcpp::abs(currentRt1 - rt2);
        double bestDotProduct = 0;
        for (int j = 0; j < mzError.size(); j++) { // Pick score with the closest dotProduct value
            if (mzError[j] > mzThreshold || rtError[j] > rtThreshold) continue; // Over the threshold
            
            // Otherwise
            // Check if the similarity score (the dot product) is closer than the last one
            // If so replace
            double dotProduct = DotProduct({mz1[i], rt1[i]}, {mz2[j], rt2[j]});
            if (dotProduct < bestDotProduct) continue;
            resultsIndexes[i] = j + 1; // To match with R indexes add 1
            bestDotProduct = dotProduct;
            if (bestDotProduct >= 1) break; // If they are equal, we should break the loop and move on
        }
        
    }
    return resultsIndexes;
}

// [[Rcpp::export]]
void GetMolecularFormula(const std::string& formula) {
    const MolecularFormula molecularFormula(formula);
}

// [[Rcpp::export]]
std::string SubtractMolecularFormula(const std::string& formula, const std::string& otherFormula) {
    const MolecularFormula molecularFormula(formula);
    const MolecularFormula otherMolecularFormula(otherFormula);
    return molecularFormula - otherMolecularFormula;
}

// [[Rcpp::export]]
bool CheckIfSubFormula(const std::string& formula, const std::string& otherFormula) {
    const MolecularFormula molecularFormula(formula);
    const MolecularFormula otherMolecularFormula(otherFormula);
    return molecularFormula.CheckIfOtherIsSubFormula(otherMolecularFormula);
}

// [[Rcpp::export]]
void FragmentationTreeTest(const Rcpp::List& molecularFormulas,
    const double parentMass, const int amountOfColors) {
    FragmentationTree tree;
    tree.AddMolecularFormulasToGraph(molecularFormulas["formula"], molecularFormulas["color"],
        molecularFormulas["scores"], parentMass, amountOfColors);
    GreedyHeuristic greedy;
   // greedy.CalculateHeuristic(tree);
}


// [[Rcpp::export]]
void test(const Rcpp::CharacterVector& vec) {
    std::vector<std::string> vecTest = Rcpp::as<std::vector<std::string>>(vec);
    Rcpp::Rcout << vec[1] << std::endl;
}
// [[Rcpp::export]]
void test2(const std::vector<std::string>& vec) {
    Rcpp::Rcout << vec[1] << std::endl;
}

// [[Rcpp::export]]
void test3(const Rcpp::String& vec) {
    Rcpp::Rcout << std::strlen(vec.get_cstring());
}


#include "../../../../Downloads/gperftools-2.15/src/gperftools/profiler.h"
#include <Rcpp.h>

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